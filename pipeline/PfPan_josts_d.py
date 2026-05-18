#!/usr/bin/env python3
"""
Fast parallel Jost's D analysis using zarr format with biological consequence prediction
"""

import os
import sys
import argparse
import subprocess
from tqdm import tqdm
import time
import pandas as pd
import numpy as np
import allel
import zarr
from multiprocessing import Pool, cpu_count
import warnings
warnings.filterwarnings('ignore')

print('scikit-allel', allel.__version__)

def load_metadata(metadata_file, sample_names):
    """Load sample metadata and align with VCF samples"""
    print("Loading metadata...")
    metadata = pd.read_csv(metadata_file)
    metadata['sample'] = metadata['sample'].str.strip()
    
    # Align with VCF samples like your code
    df_samples = pd.DataFrame({'sample': sample_names})
    df_samples = df_samples.merge(metadata, on='sample', how='left')
    
    print(f"Loaded metadata for {len(df_samples)} VCF samples")
    print(f"Samples with region info: {df_samples['Region'].notna().sum()}")
    print(f"Regions: {df_samples['Region'].value_counts().to_dict()}")
    
    return df_samples

def load_watchlist(watchlist_file):
    """Load watchlist genes"""
    print("Loading watchlist genes...")
    with open(watchlist_file, 'r') as f:
        watchlist = [line.strip().replace('gene:', '') for line in f if line.strip()]
    print(f"Loaded {len(watchlist)} watchlist genes")
    return watchlist

def load_gff_comprehensive(gff_file):
    """Load and parse GFF file with comprehensive feature extraction"""
    print("Loading comprehensive GFF annotations...")
    
    features = {
        'genes': [],
        'exons': [],
        'cds': [],
        'utrs': [],
        'repeats': []
    }
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
                
            seqid, source, feature_type, start, end, score, strand, phase, attributes = parts
            
            # Extract attributes
            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attr_dict[key] = value
            
            feature_data = {
                'chr': seqid,
                'start': int(start),
                'end': int(end),
                'strand': strand,
                'phase': phase if phase != '.' else None,
                'attributes': attr_dict
            }
            
            if feature_type == 'gene':
                gene_id = attr_dict.get('ID', attr_dict.get('gene_id', ''))
                gene_name = attr_dict.get('Name', attr_dict.get('gene_name', gene_id))
                
                feature_data.update({
                    'gene_id': gene_id,
                    'gene_name': gene_name,
                    'biotype': attr_dict.get('biotype', attr_dict.get('gene_biotype', 'protein_coding'))
                })
                features['genes'].append(feature_data)
            
            elif feature_type in ['exon', 'CDS']:
                parent_id = attr_dict.get('Parent', attr_dict.get('gene_id', ''))
                feature_data.update({
                    'parent_id': parent_id,
                    'exon_number': attr_dict.get('exon_number', attr_dict.get('rank', '')),
                    'feature_type': feature_type
                })
                
                if feature_type == 'exon':
                    features['exons'].append(feature_data)
                elif feature_type == 'CDS':
                    features['cds'].append(feature_data)
            
            elif feature_type in ['five_prime_UTR', '5UTR', 'three_prime_UTR', '3UTR', 'UTR']:
                parent_id = attr_dict.get('Parent', attr_dict.get('gene_id', ''))
                feature_data.update({
                    'parent_id': parent_id,
                    'utr_type': feature_type
                })
                features['utrs'].append(feature_data)
            
            elif 'repeat' in feature_type.lower() or feature_type in ['tandem_repeat', 'microsatellite']:
                repeat_info = {
                    'chr': seqid,
                    'start': int(start),
                    'end': int(end),
                    'repeat_type': feature_type,
                    'repeat_unit': attr_dict.get('repeat_unit', ''),
                    'repeat_count': attr_dict.get('repeat_count', ''),
                    'motif': attr_dict.get('motif', attr_dict.get('sequence', ''))
                }
                features['repeats'].append(repeat_info)
    
    # Convert to DataFrames
    feature_dfs = {}
    for feature_type, feature_list in features.items():
        if feature_list:
            feature_dfs[feature_type] = pd.DataFrame(feature_list)
            print(f"Loaded {len(feature_list)} {feature_type}")
        else:
            feature_dfs[feature_type] = pd.DataFrame()
    
    return feature_dfs

def predict_variant_consequence(variant_row, feature_dfs):
    """Predict biological consequence of a variant"""
    chr_name = variant_row['CHR']
    pos = variant_row['POS']
    ref = str(variant_row['REF'])
    alt = str(variant_row['ALT'])
    svlen = variant_row.get('SVLEN', len(alt) - len(ref))
    
    consequences = []
    details = {}
    
    # Basic variant classification
    if len(ref) == 1 and len(alt) == 1:
        var_type = 'SNV'
    elif len(ref) > len(alt):
        var_type = 'deletion'
    elif len(ref) < len(alt):
        var_type = 'insertion'
    else:
        var_type = 'complex'
    
    details['variant_type'] = var_type
    details['length_change'] = len(alt) - len(ref)
    
    # Check if it's a large structural variant
    if abs(svlen) >= 50:
        details['is_structural_variant'] = True
        if abs(svlen) >= 1000:
            consequences.append('large_structural_variant')
        else:
            consequences.append('structural_variant')
    else:
        details['is_structural_variant'] = False
    
    # Gene-level analysis
    genes_df = feature_dfs.get('genes', pd.DataFrame())
    if not genes_df.empty:
        chr_genes = genes_df[genes_df['chr'] == chr_name]
        
        # Find overlapping genes
        overlapping_genes = chr_genes[
            (chr_genes['start'] <= pos) & (chr_genes['end'] >= pos)
        ]
        
        if not overlapping_genes.empty:
            gene_info = overlapping_genes.iloc[0]
            details['gene_id'] = gene_info['gene_id']
            details['gene_name'] = gene_info['gene_name']
            details['gene_biotype'] = gene_info.get('biotype', 'protein_coding')
            
            # CDS analysis for protein-coding genes
            if gene_info.get('biotype', 'protein_coding') == 'protein_coding':
                cds_consequences = analyze_cds_consequence(
                    pos, ref, alt, gene_info, feature_dfs.get('cds', pd.DataFrame())
                )
                consequences.extend(cds_consequences['consequences'])
                details.update(cds_consequences['details'])
            
            # Exon analysis
            exon_consequences = analyze_exon_consequence(
                pos, ref, alt, gene_info, feature_dfs.get('exons', pd.DataFrame())
            )
            consequences.extend(exon_consequences['consequences'])
            details.update(exon_consequences['details'])
            
        else:
            # Intergenic variant - check for regulatory regions
            consequences.append('intergenic')
            
            # Check for proximity to genes (regulatory potential)
            if not chr_genes.empty:
                distances = np.minimum(
                    np.abs(chr_genes['start'] - pos),
                    np.abs(chr_genes['end'] - pos)
                )
                min_distance = distances.min()
                
                if min_distance <= 2000:
                    consequences.append('regulatory_region')
                    details['nearest_gene_distance'] = int(min_distance)
                    nearest_gene = chr_genes.loc[distances.idxmin()]
                    details['nearest_gene'] = nearest_gene['gene_name']
    
    # Repeat region analysis
    repeats_df = feature_dfs.get('repeats', pd.DataFrame())
    if not repeats_df.empty:
        repeat_analysis = analyze_repeat_consequence(pos, ref, alt, repeats_df, chr_name)
        if repeat_analysis:
            consequences.extend(repeat_analysis['consequences'])
            details.update(repeat_analysis['details'])
    
    # UTR analysis
    utrs_df = feature_dfs.get('utrs', pd.DataFrame())
    if not utrs_df.empty:
        utr_analysis = analyze_utr_consequence(pos, chr_name, utrs_df)
        if utr_analysis:
            consequences.extend(utr_analysis['consequences'])
            details.update(utr_analysis['details'])
    
    # Overall severity assessment
    severity = assess_variant_severity(consequences, details)
    details['predicted_severity'] = severity
    
    return {
        'consequences': ';'.join(consequences) if consequences else 'unknown',
        'details': details
    }

def analyze_cds_consequence(pos, ref, alt, gene_info, cds_df):
    """Analyze consequence within coding sequence"""
    if cds_df.empty:
        return {'consequences': [], 'details': {}}
    
    consequences = []
    details = {}
    
    # Find overlapping CDS
    gene_cds = cds_df[cds_df['parent_id'].str.contains(gene_info['gene_id'], na=False)]
    overlapping_cds = gene_cds[
        (gene_cds['start'] <= pos) & (gene_cds['end'] >= pos)
    ]
    
    if not overlapping_cds.empty:
        cds_info = overlapping_cds.iloc[0]
        details['in_cds'] = True
        
        # Calculate position within CDS
        if gene_info['strand'] == '+':
            cds_pos = pos - cds_info['start'] + 1
        else:
            cds_pos = cds_info['end'] - pos + 1
        
        details['cds_position'] = cds_pos
        
        # Determine codon position
        codon_pos = cds_pos % 3
        if codon_pos == 0:
            codon_pos = 3
        details['codon_position'] = codon_pos
        
        # Length change analysis for indels
        length_change = len(alt) - len(ref)
        
        if length_change == 0:
            # SNV in CDS
            if codon_pos == 3:  # Third codon position - often synonymous
                consequences.append('coding_sequence_variant')
                details['likely_synonymous'] = True
            else:
                consequences.append('missense_variant')
                details['likely_synonymous'] = False
        
        elif length_change % 3 == 0:
            # In-frame indel
            if length_change > 0:
                consequences.append('inframe_insertion')
                details['inframe_change'] = True
                details['amino_acid_change'] = length_change // 3
            else:
                consequences.append('inframe_deletion')
                details['inframe_change'] = True
                details['amino_acid_change'] = abs(length_change) // 3
        
        else:
            # Frameshift
            consequences.append('frameshift_variant')
            details['inframe_change'] = False
            details['frameshift'] = True
            
            # Predict premature stop
            if length_change < 0:
                consequences.append('protein_truncating_variant')
    
    else:
        details['in_cds'] = False
    
    return {'consequences': consequences, 'details': details}

def analyze_exon_consequence(pos, ref, alt, gene_info, exons_df):
    """Analyze consequence within exons"""
    if exons_df.empty:
        return {'consequences': [], 'details': {}}
    
    consequences = []
    details = {}
    
    # Find overlapping exons
    gene_exons = exons_df[exons_df['parent_id'].str.contains(gene_info['gene_id'], na=False)]
    overlapping_exons = gene_exons[
        (gene_exons['start'] <= pos) & (gene_exons['end'] >= pos)
    ]
    
    if not overlapping_exons.empty:
        exon_info = overlapping_exons.iloc[0]
        details['in_exon'] = True
        details['exon_number'] = exon_info.get('exon_number', 'unknown')
        
        # Check if near splice sites
        distance_to_start = abs(pos - exon_info['start'])
        distance_to_end = abs(pos - exon_info['end'])
        
        if distance_to_start <= 2 or distance_to_end <= 2:
            consequences.append('splice_region_variant')
            details['near_splice_site'] = True
            details['splice_distance'] = min(distance_to_start, distance_to_end)
        
        if distance_to_start == 1 or distance_to_end == 1:
            consequences.append('splice_site_variant')
            details['affects_splice_site'] = True
    
    else:
        details['in_exon'] = False
        
        # Check if in intron
        if not gene_exons.empty:
            # Find flanking exons
            upstream_exons = gene_exons[gene_exons['end'] < pos]
            downstream_exons = gene_exons[gene_exons['start'] > pos]
            
            if not upstream_exons.empty and not downstream_exons.empty:
                consequences.append('intron_variant')
                details['in_intron'] = True
                
                # Check proximity to splice sites
                nearest_upstream = upstream_exons.loc[upstream_exons['end'].idxmax()]
                nearest_downstream = downstream_exons.loc[downstream_exons['start'].idxmin()]
                
                dist_to_upstream_splice = pos - nearest_upstream['end']
                dist_to_downstream_splice = nearest_downstream['start'] - pos
                
                if dist_to_upstream_splice <= 10 or dist_to_downstream_splice <= 10:
                    consequences.append('splice_region_variant')
                    details['near_splice_site'] = True
    
    return {'consequences': consequences, 'details': details}

def analyze_repeat_consequence(pos, ref, alt, repeats_df, chr_name):
    """Analyze if variant affects tandem repeats (VNTRs, microsatellites)"""
    if repeats_df.empty:
        return None
    
    chr_repeats = repeats_df[repeats_df['chr'] == chr_name]
    overlapping_repeats = chr_repeats[
        (chr_repeats['start'] <= pos) & (chr_repeats['end'] >= pos)
    ]
    
    if not overlapping_repeats.empty:
        repeat_info = overlapping_repeats.iloc[0]
        
        consequences = []
        details = {}
        
        repeat_type = repeat_info.get('repeat_type', 'tandem_repeat')
        repeat_unit = repeat_info.get('repeat_unit', '')
        motif = repeat_info.get('motif', repeat_unit)
        
        details['in_tandem_repeat'] = True
        details['repeat_type'] = repeat_type
        details['repeat_motif'] = motif
        details['repeat_unit'] = repeat_unit
        
        # Classify repeat type
        if repeat_type == 'microsatellite' or (motif and len(motif) <= 6):
            consequences.append('microsatellite_variant')
            details['is_microsatellite'] = True
            
            # Check if it's a simple repeat expansion/contraction
            length_change = len(alt) - len(ref)
            if motif and length_change % len(motif) == 0:
                if length_change > 0:
                    consequences.append('tandem_repeat_expansion')
                    details['repeat_expansion'] = True
                    details['repeat_units_added'] = length_change // len(motif)
                elif length_change < 0:
                    consequences.append('tandem_repeat_contraction')
                    details['repeat_contraction'] = True
                    details['repeat_units_lost'] = abs(length_change) // len(motif)
        
        elif 'VNTR' in repeat_type.upper() or repeat_info.get('repeat_count', ''):
            consequences.append('vntr_variant')
            details['is_vntr'] = True
        
        else:
            consequences.append('repeat_region_variant')
        
        return {'consequences': consequences, 'details': details}
    
    return None

def analyze_utr_consequence(pos, chr_name, utrs_df):
    """Analyze consequence in UTR regions"""
    if utrs_df.empty:
        return None
    
    chr_utrs = utrs_df[utrs_df['chr'] == chr_name]
    overlapping_utrs = chr_utrs[
        (chr_utrs['start'] <= pos) & (chr_utrs['end'] >= pos)
    ]
    
    if not overlapping_utrs.empty:
        utr_info = overlapping_utrs.iloc[0]
        utr_type = utr_info.get('utr_type', 'UTR')
        
        consequences = []
        details = {}
        
        if 'five' in utr_type or '5' in utr_type:
            consequences.append('5_prime_UTR_variant')
            details['in_5_prime_utr'] = True
        elif 'three' in utr_type or '3' in utr_type:
            consequences.append('3_prime_UTR_variant')
            details['in_3_prime_utr'] = True
        else:
            consequences.append('UTR_variant')
            details['in_utr'] = True
        
        details['utr_type'] = utr_type
        
        return {'consequences': consequences, 'details': details}
    
    return None

def assess_variant_severity(consequences, details):
    """Assess overall variant severity based on consequences"""
    high_impact = [
        'frameshift_variant', 'protein_truncating_variant', 'splice_site_variant',
        'large_structural_variant'
    ]
    
    moderate_impact = [
        'missense_variant', 'inframe_insertion', 'inframe_deletion',
        'structural_variant', 'splice_region_variant'
    ]
    
    low_impact = [
        'coding_sequence_variant', 'synonymous_variant', '5_prime_UTR_variant',
        '3_prime_UTR_variant', 'intron_variant'
    ]
    
    modifier = [
        'intergenic', 'regulatory_region', 'repeat_region_variant',
        'microsatellite_variant', 'vntr_variant'
    ]
    
    if any(cons in consequences for cons in high_impact):
        return 'HIGH'
    elif any(cons in consequences for cons in moderate_impact):
        return 'MODERATE'
    elif any(cons in consequences for cons in low_impact):
        return 'LOW'
    elif any(cons in consequences for cons in modifier):
        return 'MODIFIER'
    else:
        return 'UNKNOWN'

def calculate_josts_d(allele_counts_pop1, allele_counts_pop2):
    """Calculate Jost's D between two populations"""
    # Get allele frequencies
    af1 = allele_counts_pop1.to_frequencies()
    af2 = allele_counts_pop2.to_frequencies()
    
    # Handle cases with different number of alleles
    if af1.shape[1] < 2 or af2.shape[1] < 2:
        return np.full(len(af1), np.nan)
    
    # For biallelic sites, use reference allele frequency
    p1 = af1[:, 0]  # Reference allele frequency in pop1
    p2 = af2[:, 0]  # Reference allele frequency in pop2
    
    # Handle invalid frequencies
    valid = ~np.isnan(p1) & ~np.isnan(p2) & (p1 >= 0) & (p1 <= 1) & (p2 >= 0) & (p2 <= 1)
    
    if not np.any(valid):
        return np.full(len(af1), np.nan)
    
    # Calculate diversity measures for biallelic sites
    # Within-population diversity: Hs = 2*p*(1-p)
    hs1 = 2 * p1 * (1 - p1)
    hs2 = 2 * p2 * (1 - p2)
    hs_mean = (hs1 + hs2) / 2
    
    # Total diversity
    p_total = (p1 + p2) / 2
    ht = 2 * p_total * (1 - p_total)
    
    # Jost's D = (Ht - Hs) / (1 - Hs)
    josts_d = np.where(
        (hs_mean >= 1) | (hs_mean == 0), 
        0, 
        np.maximum(0, (ht - hs_mean) / (1 - hs_mean))
    )
    
    # Set invalid positions to NaN
    josts_d[~valid] = np.nan
    
    return josts_d

def process_chromosome_chunk(args):
    """Process a chunk of variants for Jost's D calculation"""
    chunk_info, callset_path, df_samples = args
    chr_name, start_idx, end_idx = chunk_info
    
    try:
        # Open zarr callset
        callset = zarr.open(callset_path, mode='r')
        
        # Get data for this chunk
        chunk_slice = slice(start_idx, end_idx)
        
        gt_chunk = allel.GenotypeChunkedArray(callset['calldata']['GT'][chunk_slice])
        pos_chunk = callset['variants/POS'][chunk_slice]
        ref_chunk = callset['variants/REF'][chunk_slice]
        alt_chunk = callset['variants/ALT'][chunk_slice]
        chrom_chunk = callset['variants/CHROM'][chunk_slice]
        
        # Get SVLEN if available
        svlen_chunk = None
        if 'variants/SVLEN' in callset:
            svlen_chunk = callset['variants/SVLEN'][chunk_slice]
        elif 'variants/INFO_SVLEN' in callset:
            svlen_chunk = callset['variants/INFO_SVLEN'][chunk_slice]
        else:
            # Calculate length difference for indels
            ref_lengths = np.array([len(str(r)) for r in ref_chunk])
            alt_lengths = np.array([len(str(a[0]) if hasattr(a, '__len__') and len(a) > 0 else str(a)) 
                                   for a in alt_chunk])
            svlen_chunk = alt_lengths - ref_lengths
        
        print(f"  Processing chunk {chr_name}:{start_idx}-{end_idx} ({len(pos_chunk)} variants)")
        
        # Get unique regions
        regions = df_samples['Region'].dropna().unique()
        chunk_results = []
        
        for focal_region in regions:
            # Get sample indices for focal vs others
            focal_mask = (df_samples['Region'] == focal_region).values
            other_mask = (df_samples['Region'] != focal_region) & df_samples['Region'].notna().values
            
            focal_indices = np.where(focal_mask)[0]
            other_indices = np.where(other_mask)[0]
            
            if len(focal_indices) == 0 or len(other_indices) == 0:
                continue
            
            # Calculate allele counts
            gt_focal = gt_chunk.take(focal_indices, axis=1)
            gt_others = gt_chunk.take(other_indices, axis=1)
            
            ac_focal = gt_focal.count_alleles()
            ac_others = gt_others.count_alleles()
            
            # Calculate Jost's D
            josts_d = calculate_josts_d(ac_focal, ac_others)
            
            # Calculate genotype statistics for focal population
            n_total_focal = len(focal_indices)
            n_ref_focal = np.sum(gt_focal.is_hom_ref(), axis=1)
            n_alt_focal = np.sum(gt_focal.is_hom_alt() | gt_focal.is_het(), axis=1)
            n_missing_focal = np.sum(gt_focal.is_missing(), axis=1)
            
            # Calculate genotype statistics for other populations combined
            n_total_others = len(other_indices)
            n_ref_others = np.sum(gt_others.is_hom_ref(), axis=1)
            n_alt_others = np.sum(gt_others.is_hom_alt() | gt_others.is_het(), axis=1)
            n_missing_others = np.sum(gt_others.is_missing(), axis=1)
            
            # Handle ALT alleles
            if alt_chunk.ndim > 1:
                alt_alleles = []
                for i in range(len(alt_chunk)):
                    alt_row = alt_chunk[i]
                    if hasattr(alt_row, '__len__') and len(alt_row) > 0:
                        # Take first non-empty alt allele
                        alt_val = alt_row[0]
                        if isinstance(alt_val, bytes):
                            alt_val = alt_val.decode()
                        if alt_val and alt_val != '.':
                            alt_alleles.append(alt_val)
                        else:
                            alt_alleles.append('.')
                    else:
                        alt_alleles.append(str(alt_row) if alt_row else '.')
            else:
                alt_alleles = [a.decode() if isinstance(a, bytes) else str(a) for a in alt_chunk]
            
            # Handle REF alleles
            ref_alleles = [r.decode() if isinstance(r, bytes) else str(r) for r in ref_chunk]
            
            # Handle SVLEN
            if svlen_chunk is not None:
                svlen_values = svlen_chunk
            else:
                svlen_values = np.zeros(len(pos_chunk))
            
            # Create results for this chunk
            chunk_result = pd.DataFrame({
                'CHR': [chr_name] * len(pos_chunk) if isinstance(chr_name, str) else 
                       [c.decode() if isinstance(c, bytes) else str(c) for c in chrom_chunk],
                'POS': pos_chunk,
                'REF': ref_alleles,
                'ALT': alt_alleles,
                'SVLEN': svlen_values,
                'focal_region': focal_region,
                'josts_d': josts_d,
                
                # Focal population stats
                'n_total_focal': n_total_focal,
                'n_ref_focal': n_ref_focal,
                'n_alt_focal': n_alt_focal,
                'n_missing_focal': n_missing_focal,
                'pct_ref_focal': np.round(n_ref_focal / n_total_focal * 100, 2),
                'pct_alt_focal': np.round(n_alt_focal / n_total_focal * 100, 2),
                'pct_missing_focal': np.round(n_missing_focal / n_total_focal * 100, 2),
                
                # Others population stats
                'n_total_others': n_total_others,
                'n_ref_others': n_ref_others,
                'n_alt_others': n_alt_others,
                'n_missing_others': n_missing_others,
                'pct_ref_others': np.round(n_ref_others / n_total_others * 100, 2),
                'pct_alt_others': np.round(n_alt_others / n_total_others * 100, 2),
                'pct_missing_others': np.round(n_missing_others / n_total_others * 100, 2),
                
                # Allele frequencies
                'af_focal': ac_focal.to_frequencies()[:, 1] if ac_focal.shape[1] > 1 else np.zeros(len(pos_chunk)),
                'af_others': ac_others.to_frequencies()[:, 1] if ac_others.shape[1] > 1 else np.zeros(len(pos_chunk))
            })
            
            chunk_results.append(chunk_result)
        
        if chunk_results:
            return pd.concat(chunk_results, ignore_index=True)
        else:
            return pd.DataFrame()
    
    except Exception as e:
        print(f"Error processing chunk {chr_name}:{start_idx}-{end_idx}: {e}")
        import traceback
        traceback.print_exc()
        return pd.DataFrame()

def create_chunks(callset, chunk_size=10000):
    """Create chunks for parallel processing"""
    n_variants = callset['variants/POS'].shape[0]
    
    # Get chromosome information
    chroms = callset['variants/CHROM'][:]
    positions = callset['variants/POS'][:]
    
    # Create chunks by chromosome
    chunks = []
    unique_chroms = np.unique(chroms)
    
    for chrom in unique_chroms:
        chrom_mask = chroms == chrom
        chrom_indices = np.where(chrom_mask)[0]
        
        if len(chrom_indices) == 0:
            continue
        
        # Split chromosome into chunks
        start_idx = chrom_indices[0]
        end_idx = chrom_indices[-1] + 1
        
        for chunk_start in range(start_idx, end_idx, chunk_size):
            chunk_end = min(chunk_start + chunk_size, end_idx)
            chunks.append((chrom.decode() if isinstance(chrom, bytes) else chrom, 
                          chunk_start, chunk_end))
    
    print(f"Created {len(chunks)} chunks across {len(unique_chroms)} chromosomes")
    return chunks

def annotate_with_biological_consequences(variants_df, feature_dfs, upstream_bp=2000, downstream_bp=2000):
    """Annotate variants with comprehensive biological consequences"""
    print(f"Annotating {len(variants_df)} variants with biological consequences...")
    
    annotated_variants = []
    
    for _, variant in tqdm(variants_df.iterrows(), total=len(variants_df), desc="Predicting consequences"):
        # Get basic gene annotation
        chr_name = variant['CHR']
        pos = variant['POS']
        
        # Find overlapping genes
        genes_df = feature_dfs.get('genes', pd.DataFrame())
        if not genes_df.empty:
            chr_genes = genes_df[genes_df['chr'] == chr_name]
            gene_overlaps = chr_genes[
                (chr_genes['start'] <= pos) & (chr_genes['end'] >= pos)
            ]
            
            if not gene_overlaps.empty:
                gene_info = gene_overlaps.iloc[0]
                variant['gene_id'] = gene_info['gene_id']
                variant['gene_name'] = gene_info['gene_name'] 
                variant['gene_biotype'] = gene_info.get('biotype', 'protein_coding')
                variant['region_type'] = 'gene'
            else:
                # Check upstream/downstream
                if not chr_genes.empty:
                    distances = np.minimum(
                        np.abs(chr_genes['start'] - pos),
                        np.abs(chr_genes['end'] - pos)
                    )
                    min_distance = distances.min()
                    
                    if min_distance <= upstream_bp:
                        nearest_gene = chr_genes.loc[distances.idxmin()]
                        variant['gene_id'] = nearest_gene['gene_id']
                        variant['gene_name'] = nearest_gene['gene_name']
                        variant['gene_biotype'] = nearest_gene.get('biotype', 'protein_coding')
                        
                        if pos < nearest_gene['start']:
                            variant['region_type'] = 'upstream'
                        elif pos > nearest_gene['end']:
                            variant['region_type'] = 'downstream'
                        else:
                            variant['region_type'] = 'intergenic'
                    else:
                        variant['gene_id'] = 'intergenic'
                        variant['gene_name'] = 'intergenic'
                        variant['gene_biotype'] = 'intergenic'
                        variant['region_type'] = 'intergenic'
                else:
                    variant['gene_id'] = 'intergenic'
                    variant['gene_name'] = 'intergenic'
                    variant['gene_biotype'] = 'intergenic'
                    variant['region_type'] = 'intergenic'
        
        # Predict biological consequence
        consequence_prediction = predict_variant_consequence(variant, feature_dfs)
        
        # Add consequence columns
        variant['predicted_consequences'] = consequence_prediction['consequences']
        
        # Add detailed consequence information
        details = consequence_prediction['details']
        variant['variant_type'] = details.get('variant_type', 'unknown')
        variant['predicted_severity'] = details.get('predicted_severity', 'UNKNOWN')
        variant['length_change'] = details.get('length_change', 0)
        variant['is_structural_variant'] = details.get('is_structural_variant', False)
        variant['in_cds'] = details.get('in_cds', False)
        variant['in_exon'] = details.get('in_exon', False)
        variant['frameshift'] = details.get('frameshift', False)
        variant['inframe_change'] = details.get('inframe_change', False)
        variant['affects_splice_site'] = details.get('affects_splice_site', False)
        variant['in_tandem_repeat'] = details.get('in_tandem_repeat', False)
        variant['is_microsatellite'] = details.get('is_microsatellite', False)
        variant['is_vntr'] = details.get('is_vntr', False)
        variant['repeat_type'] = details.get('repeat_type', '')
        variant['repeat_motif'] = details.get('repeat_motif', '')
        variant['repeat_expansion'] = details.get('repeat_expansion', False)
        variant['repeat_contraction'] = details.get('repeat_contraction', False)
        variant['amino_acid_change'] = details.get('amino_acid_change', 0)
        variant['codon_position'] = details.get('codon_position', 0)
        variant['likely_synonymous'] = details.get('likely_synonymous', False)
        
        annotated_variants.append(variant)
    
    return pd.DataFrame(annotated_variants)

def filter_and_annotate_variants(combined_results, feature_dfs, watchlist, args):
    """Filter variants and create annotated tables with biological consequences"""
    
    # Calculate thresholds
    print("Calculating Jost's D thresholds...")
    thresholds = combined_results.groupby('focal_region')['josts_d'].quantile(args.percentile/100).reset_index()
    thresholds.columns = ['focal_region', 'jd_threshold']
    
    combined_results = combined_results.merge(thresholds, on='focal_region')
    combined_results['high_josts_d'] = (combined_results['josts_d'] >= combined_results['jd_threshold']) & \
                                      combined_results['josts_d'].notna()
    
    # Basic quality filters
    print("Applying quality filters...")
    filtered_results = combined_results[
        (combined_results['pct_missing_focal'] < 50) &
        (combined_results['n_alt_focal'] > 5) &
        (combined_results['josts_d'].notna())
    ].copy()
    
    print(f"After quality filters: {len(filtered_results)} variants")
    
    # Get high Jost's D variants
    high_jd = filtered_results[filtered_results['high_josts_d'] == True].copy()
    print(f"High Jost's D variants: {len(high_jd)}")
    
    # Get structural variants (>=50bp)
    sv_variants = filtered_results[np.abs(filtered_results['SVLEN']) >= 50].copy()
    print(f"Structural variants (>=50bp): {len(sv_variants)}")
    
    # High Jost's D SVs
    high_jd_sv = high_jd[np.abs(high_jd['SVLEN']) >= 50].copy()
    print(f"High Jost's D structural variants: {len(high_jd_sv)}")
    
    # High MAF variants (increased threshold)
    high_maf = filtered_results[filtered_results['pct_alt_focal'] >= args.maf_threshold].copy()
    print(f"High MAF variants (>={args.maf_threshold}%): {len(high_maf)}")
    
    # High MAF SVs
    high_maf_sv = high_maf[np.abs(high_maf['SVLEN']) >= 50].copy()
    print(f"High MAF structural variants: {len(high_maf_sv)}")
    
    # Annotate variants with genes and biological consequences (limit to reasonable sizes)
    results_dict = {}
    
    # Annotate high Jost's D variants
    if not high_jd.empty and len(high_jd) <= args.max_annotate:
        print("Annotating high Jost's D variants with biological consequences...")
        high_jd_annotated = annotate_with_biological_consequences(high_jd, feature_dfs)
        high_jd_annotated['is_watchlist'] = high_jd_annotated.apply(
            lambda row: any(w in str(row['gene_name']) or w in str(row['gene_id']) 
                           for w in watchlist), axis=1
        )
        results_dict['high_josts_d_annotated'] = high_jd_annotated
        results_dict['watchlist_high_jd'] = high_jd_annotated[high_jd_annotated['is_watchlist'] == True]
    else:
        print(f"Skipping detailed annotation for high Jost's D variants ({len(high_jd)} > {args.max_annotate})")
        results_dict['high_josts_d_annotated'] = high_jd
        results_dict['watchlist_high_jd'] = pd.DataFrame()
    
    # Annotate high Jost's D SVs
    if not high_jd_sv.empty and len(high_jd_sv) <= args.max_annotate:
        print("Annotating high Jost's D structural variants...")
        high_jd_sv_annotated = annotate_with_biological_consequences(high_jd_sv, feature_dfs)
        high_jd_sv_annotated['is_watchlist'] = high_jd_sv_annotated.apply(
            lambda row: any(w in str(row['gene_name']) or w in str(row['gene_id']) 
                           for w in watchlist), axis=1
        )
        results_dict['high_josts_d_sv_annotated'] = high_jd_sv_annotated
        results_dict['watchlist_high_jd_sv'] = high_jd_sv_annotated[high_jd_sv_annotated['is_watchlist'] == True]
    else:
        print(f"Skipping annotation for high Jost's D SVs ({len(high_jd_sv)} variants)")
        results_dict['high_josts_d_sv_annotated'] = high_jd_sv
        results_dict['watchlist_high_jd_sv'] = pd.DataFrame()
    
    # Annotate high MAF variants  
    if not high_maf.empty and len(high_maf) <= args.max_annotate:
        print("Annotating high MAF variants...")
        high_maf_annotated = annotate_with_biological_consequences(high_maf, feature_dfs)
        high_maf_annotated['is_watchlist'] = high_maf_annotated.apply(
            lambda row: any(w in str(row['gene_name']) or w in str(row['gene_id'])
                           for w in watchlist), axis=1
        )
        results_dict['high_maf_annotated'] = high_maf_annotated
        results_dict['watchlist_high_maf'] = high_maf_annotated[high_maf_annotated['is_watchlist'] == True]
    else:
        print(f"Skipping annotation for high MAF variants ({len(high_maf)} > {args.max_annotate})")
        results_dict['high_maf_annotated'] = high_maf
        results_dict['watchlist_high_maf'] = pd.DataFrame()
    
    # Annotate high MAF SVs
    if not high_maf_sv.empty and len(high_maf_sv) <= args.max_annotate:
        print("Annotating high MAF structural variants...")
        high_maf_sv_annotated = annotate_with_biological_consequences(high_maf_sv, feature_dfs)
        high_maf_sv_annotated['is_watchlist'] = high_maf_sv_annotated.apply(
            lambda row: any(w in str(row['gene_name']) or w in str(row['gene_id'])
                           for w in watchlist), axis=1
        )
        results_dict['high_maf_sv_annotated'] = high_maf_sv_annotated
        results_dict['watchlist_high_maf_sv'] = high_maf_sv_annotated[high_maf_sv_annotated['is_watchlist'] == True]
    else:
        print(f"Skipping annotation for high MAF SVs ({len(high_maf_sv)} variants)")
        results_dict['high_maf_sv_annotated'] = high_maf_sv
        results_dict['watchlist_high_maf_sv'] = pd.DataFrame()
    
    # Add basic tables
    results_dict['all_variants'] = combined_results
    results_dict['filtered_variants'] = filtered_results
    results_dict['structural_variants'] = sv_variants
    
    return results_dict

def main(args):
    print("="*60)
    print("JOST'S D ANALYSIS WITH BIOLOGICAL CONSEQUENCE PREDICTION")
    print("="*60)
    
    # Convert VCF to zarr if needed
    if not os.path.exists(args.zarr):
        print(f"Converting VCF to zarr: {args.vcf} -> {args.zarr}")
        allel.vcf_to_zarr(args.vcf, args.zarr)
    
    # Open zarr callset
    print("Opening zarr callset...")
    callset = zarr.open(args.zarr, mode='r')
    
    # Get sample names (like your code)
    sample_names = callset['samples'][:]
    if isinstance(sample_names[0], bytes):
        sample_names = [s.decode() for s in sample_names]
    
    print(f"Found {len(sample_names)} samples in zarr file")
    
    # Check available fields
    print("Available fields in callset:")
    for key in sorted(callset.keys()):
        print(f"  {key}: {callset[key].shape if hasattr(callset[key], 'shape') else 'group'}")
    
    # Load metadata and align with samples
    df_samples = load_metadata(args.metadata, sample_names)
    
    # Load other data
    watchlist = load_watchlist(args.watchlist)
    feature_dfs = load_gff_comprehensive(args.gff)
    
    # Filter samples with region information
    df_samples_with_region = df_samples.dropna(subset=['Region'])
    if len(df_samples_with_region) == 0:
        print("ERROR: No samples have region information!")
        return
    
    print(f"Using {len(df_samples_with_region)} samples with region information")
    
    # Create chunks for processing
    chunks = create_chunks(callset, chunk_size=args.chunk_size)
    
    # Prepare arguments for parallel processing
    pool_args = [(chunk, args.zarr, df_samples) for chunk in chunks]
    
    # Process chunks in parallel
    print(f"Processing {len(chunks)} chunks with {args.threads} threads...")
    with Pool(args.threads) as pool:
        chunk_results = []
        for result in tqdm(pool.imap(process_chromosome_chunk, pool_args), 
                          total=len(pool_args), desc="Processing chunks"):
            if not result.empty:
                chunk_results.append(result)
    
    if not chunk_results:
        print("No results generated!")
        return
    
    # Combine all results
    print("Combining chunk results...")
    combined_results = pd.concat(chunk_results, ignore_index=True)
    print(f"Combined results: {len(combined_results)} variant-population combinations")
    
    # Filter and annotate variants with biological consequences
    print("Filtering and annotating variants with biological consequences...")
    results_dict = filter_and_annotate_variants(combined_results, feature_dfs, watchlist, args)
    
    # Save results
    print("\nSaving results...")
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Save all result tables
    for name, df in results_dict.items():
        if not df.empty:
            filename = f"{name}.csv"
            filepath = os.path.join(args.output_dir, filename)
            df.to_csv(filepath, index=False)
            print(f"Saved {len(df)} rows to {filename}")
    
    # Create summary of biological consequences
    print("\nCreating biological consequence summary...")
    consequence_summary = {}
    
    for name, df in results_dict.items():
        if 'annotated' in name and 'predicted_consequences' in df.columns:
            # Count consequence types
            all_consequences = []
            for cons_str in df['predicted_consequences'].dropna():
                all_consequences.extend(cons_str.split(';'))
            
            cons_counts = pd.Series(all_consequences).value_counts()
            consequence_summary[name] = cons_counts
            
            # Summary by severity
            if 'predicted_severity' in df.columns:
                severity_counts = df['predicted_severity'].value_counts()
                print(f"\n{name} - Severity distribution:")
                for severity, count in severity_counts.items():
                    print(f"  {severity}: {count}")
            
            # Functional variant counts
            if len(df) > 0:
                functional_counts = {
                    'Frameshift variants': df['frameshift'].sum() if 'frameshift' in df.columns else 0,
                    'In-frame changes': df['inframe_change'].sum() if 'inframe_change' in df.columns else 0,
                    'Splice site variants': df['affects_splice_site'].sum() if 'affects_splice_site' in df.columns else 0,
                    'Microsatellites': df['is_microsatellite'].sum() if 'is_microsatellite' in df.columns else 0,
                    'VNTRs': df['is_vntr'].sum() if 'is_vntr' in df.columns else 0,
                    'Repeat expansions': df['repeat_expansion'].sum() if 'repeat_expansion' in df.columns else 0,
                }
                
                print(f"\n{name} - Functional variant types:")
                for var_type, count in functional_counts.items():
                    if count > 0:
                        print(f"  {var_type}: {count}")
    
    # Save consequence summary
    consequence_summary_df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in consequence_summary.items()]))
    consequence_summary_df.to_csv(os.path.join(args.output_dir, 'consequence_summary.csv'))
    
    # Print summary
    print("\n" + "="*60)
    print("ANALYSIS SUMMARY")
    print("="*60)
    print(f"Total variant-population combinations: {len(results_dict['all_variants'])}")
    print(f"After quality filters: {len(results_dict['filtered_variants'])}")
    print(f"High Jost's D variants ({args.percentile}th percentile): {len(results_dict['high_josts_d_annotated'])}")
    print(f"Structural variants (>=50bp): {len(results_dict['structural_variants'])}")
    print(f"High Jost's D structural variants: {len(results_dict['high_josts_d_sv_annotated'])}")
    print(f"High MAF variants (>={args.maf_threshold}%): {len(results_dict['high_maf_annotated'])}")
    print(f"High MAF structural variants: {len(results_dict['high_maf_sv_annotated'])}")
    print(f"Watchlist high Jost's D: {len(results_dict['watchlist_high_jd'])}")
    print(f"Watchlist high Jost's D SVs: {len(results_dict['watchlist_high_jd_sv'])}")
    print(f"Watchlist high MAF: {len(results_dict['watchlist_high_maf'])}")
    print(f"Watchlist high MAF SVs: {len(results_dict['watchlist_high_maf_sv'])}")
    
    # Print top results for watchlist with biological consequences
    def print_top_results_with_consequences(df, title, sort_col):
        if not df.empty and 'predicted_consequences' in df.columns:
            print(f"\n{title}:")
            display_cols = ['CHR', 'POS', 'focal_region', 'gene_name', 'josts_d', 'pct_alt_focal', 
                           'SVLEN', 'variant_type', 'predicted_consequences', 'predicted_severity']
            available_cols = [col for col in display_cols if col in df.columns]
            top_results = df.nlargest(min(5, len(df)), sort_col)[available_cols]
            print(top_results.to_string(index=False, max_colwidth=30))
    
    # Display top results with consequences
    print_top_results_with_consequences(results_dict['watchlist_high_jd'], 
                                       "Top 5 watchlist variants by Jost's D (with consequences)", 
                                       'josts_d')
    
    print_top_results_with_consequences(results_dict['watchlist_high_jd_sv'], 
                                       "Top 5 watchlist structural variants by Jost's D", 
                                       'josts_d')
    
    print(f"\nResults saved to: {args.output_dir}")
    print("\nKey files generated:")
    print("  - *_annotated.csv: Variants with biological consequence predictions")
    print("  - consequence_summary.csv: Summary of predicted consequences")
    print("  - watchlist_*: Drug resistance gene variants")
    
    print("\nBiological consequence columns added:")
    print("  - predicted_consequences: Functional impact prediction")
    print("  - predicted_severity: HIGH/MODERATE/LOW/MODIFIER")
    print("  - variant_type: SNV/insertion/deletion/complex")
    print("  - frameshift: Causes protein frameshift")
    print("  - inframe_change: In-frame insertion/deletion")
    print("  - affects_splice_site: Near splice junctions")
    print("  - is_microsatellite/is_vntr: Tandem repeat variants")
    print("  - repeat_expansion/contraction: Repeat length changes")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Josts D analysis with biological consequence prediction')
    parser.add_argument('--vcf', required=True, help='VCF file path')
    parser.add_argument('--zarr', required=True, help='Zarr file path (will be created if not exists)')
    parser.add_argument('--metadata', required=True, help='Metadata CSV file')
    parser.add_argument('--gff', required=True, help='GFF file')
    parser.add_argument('--watchlist', required=True, help='Watchlist genes file')
    parser.add_argument('--output-dir', default='josts_d_biological_results', help='Output directory')
    parser.add_argument('--threads', type=int, default=cpu_count(), help='Number of threads')
    parser.add_argument('--percentile', type=float, default=99, help='Percentile threshold for high Josts D')
    parser.add_argument('--chunk-size', type=int, default=10000, help='Variants per chunk')
    parser.add_argument('--maf-threshold', type=float, default=20.0, help='MAF threshold percentage')
    parser.add_argument('--max-annotate', type=int, default=50000, help='Maximum variants to annotate')
    
    args = parser.parse_args()
    main(args)
