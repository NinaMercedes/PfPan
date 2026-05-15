#!/usr/bin/env python3
"""
Mosdepth CNV Validation Using Flanking Region Normalization
Per-sample normalization using flanking regions instead of population-based approach

Corrections vs original script:
1. Removed MDR1_199kb pattern (left flank sits in low-coverage subtelomeric region,
   artificially inflating coverage ratios)
2. Added minimum flank coverage threshold (MIN_FLANK_COVERAGE) to skip unreliable patterns
3. Added flank consistency check (MAX_FLANK_RATIO) to flag coverage gradients
4. Raised low-confidence aggregation requirement from 1 to 2 patterns

Output terminology uses concordance/discordance rather than TP/TN/FP/FN, as
MalariaGen linear mapping calls are used as a reference comparator rather than
a gold standard truth set. Discordant calls reflect method differences rather
than definitive errors in either pipeline.
"""

import pandas as pd
import numpy as np
import glob
import os
from pathlib import Path
import argparse

# ── Quality thresholds ────────────────────────────────────────────────────────
MIN_FLANK_COVERAGE    = 5   # minimum mean flank coverage to trust a ratio
MAX_FLANK_RATIO       = 1.5 # maximum left/right flank imbalance before flagging
MIN_LOW_CONF_PATTERNS = 2   # number of low-confidence patterns required to call amplification
# ─────────────────────────────────────────────────────────────────────────────


def define_breakpoints_with_flanks():
    """Define breakpoint patterns with flanking regions.

    MDR1_199kb removed: its left flank (760906-780906) falls in a subtelomeric
    low-coverage region, producing spuriously high coverage ratios (~2.17) that
    are indistinguishable from genuine amplification but appear identically in
    both linear and pangenome alignments, indicating a mappability artefact.
    GCH1 and CRT excluded: GCH1 is structurally complex and not collapsed to a
    single copy in assemblies; CRT amplification is rare and not the focus here.
    """

    patterns = [
        # MDR1 patterns (199kb removed)
        {
            'gene': 'MDR1', 'pattern': 'MDR1_17kb', 'chrom': 'Pf3D7#0#Pf3D7_05_v3',
            'target_start': 949514, 'target_end': 967266,
            'left_start': 929514, 'left_end': 949514,
            'right_start': 967266, 'right_end': 987266
        },
        {
            'gene': 'MDR1', 'pattern': 'MDR1_42kb', 'chrom': 'Pf3D7#0#Pf3D7_05_v3',
            'target_start': 938329, 'target_end': 980040,
            'left_start': 918329, 'left_end': 938329,
            'right_start': 980040, 'right_end': 1000040
        },
        {
            'gene': 'MDR1', 'pattern': 'MDR1_15kb', 'chrom': 'Pf3D7#0#Pf3D7_05_v3',
            'target_start': 947790, 'target_end': 962454,
            'left_start': 927790, 'left_end': 947790,
            'right_start': 962454, 'right_end': 982454
        },
        {
            'gene': 'MDR1', 'pattern': 'MDR1_19kb', 'chrom': 'Pf3D7#0#Pf3D7_05_v3',
            'target_start': 953961, 'target_end': 973034,
            'left_start': 933961, 'left_end': 953961,
            'right_start': 973034, 'right_end': 993034
        },

        # Plasmepsin patterns
        {
            'gene': 'plasmepsin', 'pattern': 'PM_17kb', 'chrom': 'Pf3D7#0#Pf3D7_14_v3',
            'target_start': 283034, 'target_end': 300522,
            'left_start': 263034, 'left_end': 283034,
            'right_start': 300522, 'right_end': 320522
        },
        {
            'gene': 'plasmepsin', 'pattern': 'PM_80kb', 'chrom': 'Pf3D7#0#Pf3D7_14_v3',
            'target_start': 283034, 'target_end': 363020,
            'left_start': 263034, 'left_end': 283034,
            'right_start': 363020, 'right_end': 383020
        },
        {
            'gene': 'plasmepsin', 'pattern': 'PM_9kb', 'chrom': 'Pf3D7#0#Pf3D7_14_v3',
            'target_start': 289611, 'target_end': 298792,
            'left_start': 269611, 'left_end': 289611,
            'right_start': 298792, 'right_end': 318792
        }
    ]

    return pd.DataFrame(patterns)


def extract_region_coverage(coverage_df, chrom, start, end):
    """Extract weighted average coverage for a genomic region."""

    overlapping = coverage_df[
        (coverage_df['chrom'] == chrom) &
        (coverage_df['end']   > start)  &
        (coverage_df['start'] < end)
    ].copy()

    if len(overlapping) == 0:
        return 0

    overlapping['overlap_start']     = np.maximum(overlapping['start'], start)
    overlapping['overlap_end']       = np.minimum(overlapping['end'],   end)
    overlapping['overlap_length']    = overlapping['overlap_end'] - overlapping['overlap_start']
    overlapping['weighted_coverage'] = overlapping['coverage'] * overlapping['overlap_length']

    total_length   = overlapping['overlap_length'].sum()
    total_weighted = overlapping['weighted_coverage'].sum()

    return total_weighted / total_length if total_length > 0 else 0


def process_mosdepth_file_with_flanks(file_path, patterns_df):
    """Process one mosdepth file extracting target and flanking region coverage."""

    name = Path(file_path).name
    for suffix in ['.regions.bed.gz', '_regions_bed.gz', '_sort_bam_regions_bed.gz',
                   '_sort_bam.regions.bed.gz', '.bed.gz']:
        name = name.replace(suffix, '')
    for infix in ['_sort.bam', '_sort_bam', '_sort']:
        name = name.replace(infix, '')
    sample_id = name

    try:
        coverage_df = pd.read_csv(
            file_path, sep='\t', header=None,
            names=['chrom', 'start', 'end', 'coverage'],
            compression='gzip'
        )

        results = []

        for _, pattern in patterns_df.iterrows():
            chrom = pattern['chrom']

            target_cov = extract_region_coverage(
                coverage_df, chrom, pattern['target_start'], pattern['target_end'])
            left_cov   = extract_region_coverage(
                coverage_df, chrom, pattern['left_start'],   pattern['left_end'])
            right_cov  = extract_region_coverage(
                coverage_df, chrom, pattern['right_start'],  pattern['right_end'])

            flank_baseline = (left_cov + right_cov) / 2 if (left_cov + right_cov) > 0 else 0

            # ── Quality checks ────────────────────────────────────────────────
            if flank_baseline < MIN_FLANK_COVERAGE:
                results.append({
                    'sample_id':            sample_id,
                    'gene':                 pattern['gene'],
                    'pattern':              pattern['pattern'],
                    'target_coverage':      target_cov,
                    'left_flank_coverage':  left_cov,
                    'right_flank_coverage': right_cov,
                    'flank_baseline':       flank_baseline,
                    'coverage_ratio':       None,
                    'cnv_call':             'unreliable_flank',
                    'confidence':           'none',
                    'qc_flag':              'low_flank_coverage',
                    'target_size':          pattern['target_end'] - pattern['target_start']
                })
                continue

            if left_cov > 0 and right_cov > 0:
                flank_imbalance = max(left_cov, right_cov) / min(left_cov, right_cov)
            else:
                flank_imbalance = 999

            if flank_imbalance > MAX_FLANK_RATIO:
                results.append({
                    'sample_id':            sample_id,
                    'gene':                 pattern['gene'],
                    'pattern':              pattern['pattern'],
                    'target_coverage':      target_cov,
                    'left_flank_coverage':  left_cov,
                    'right_flank_coverage': right_cov,
                    'flank_baseline':       flank_baseline,
                    'coverage_ratio':       None,
                    'cnv_call':             'uneven_flanks',
                    'confidence':           'none',
                    'qc_flag':              f'flank_imbalance_{flank_imbalance:.2f}',
                    'target_size':          pattern['target_end'] - pattern['target_start']
                })
                continue
            # ─────────────────────────────────────────────────────────────────

            coverage_ratio = target_cov / flank_baseline if flank_baseline > 0 else 0

            # ── CNV calling thresholds ────────────────────────────────────────
            if target_cov == 0 or flank_baseline == 0:
                cnv_call, confidence = 'no_coverage', 'low'
            elif coverage_ratio >= 2.0:
                cnv_call, confidence = 'amplification', 'high'
            elif coverage_ratio >= 1.5:
                cnv_call, confidence = 'amplification', 'medium'
            elif coverage_ratio >= 1.3:
                cnv_call, confidence = 'amplification', 'low'
            elif coverage_ratio <= 0.5:
                cnv_call, confidence = 'deletion', 'high'
            elif coverage_ratio <= 0.7:
                cnv_call, confidence = 'deletion', 'medium'
            else:
                cnv_call, confidence = 'normal', 'high'
            # ─────────────────────────────────────────────────────────────────

            results.append({
                'sample_id':            sample_id,
                'gene':                 pattern['gene'],
                'pattern':              pattern['pattern'],
                'target_coverage':      target_cov,
                'left_flank_coverage':  left_cov,
                'right_flank_coverage': right_cov,
                'flank_baseline':       flank_baseline,
                'coverage_ratio':       coverage_ratio,
                'cnv_call':             cnv_call,
                'confidence':           confidence,
                'qc_flag':              'pass',
                'target_size':          pattern['target_end'] - pattern['target_start']
            })

        return pd.DataFrame(results)

    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        return pd.DataFrame()


def call_cnvs_from_ratios(coverage_data):
    """Call CNVs using flanking region ratios (pattern level)."""

    print("CALLING CNVs FROM FLANKING RATIOS")
    print("=" * 50)

    passing = coverage_data[coverage_data['qc_flag'] == 'pass'].copy()
    failed  = coverage_data[coverage_data['qc_flag'] != 'pass']

    print(f"Patterns passing QC:  {len(passing):,}")
    print(f"Patterns failing QC:  {len(failed):,}")
    if len(failed) > 0:
        print("QC failure reasons:")
        for reason, count in failed['qc_flag'].value_counts().items():
            print(f"  {reason}: {count:,}")

    print()
    total = len(passing)
    if total == 0:
        print("No passing patterns — check QC thresholds.")
        return coverage_data

    amps   = len(passing[passing['cnv_call'] == 'amplification'])
    dels   = len(passing[passing['cnv_call'] == 'deletion'])
    normal = len(passing[passing['cnv_call'] == 'normal'])

    print(f"Total passing calls: {total:,}")
    print(f"Amplifications:      {amps:,}  ({100*amps/total:.1f}%)")
    print(f"Deletions:           {dels:,}  ({100*dels/total:.1f}%)")
    print(f"Normal:              {normal:,} ({100*normal/total:.1f}%)")
    print()
    print("Per-gene amplification rates:")
    for gene in passing['gene'].unique():
        gd      = passing[passing['gene'] == gene]
        g_amps  = len(gd[gd['cnv_call'] == 'amplification'])
        g_total = len(gd)
        if g_total > 0:
            print(f"  {gene}: {g_amps}/{g_total} ({100*g_amps/g_total:.1f}%)")

    return coverage_data


def aggregate_to_gene_calls(cnv_df):
    """Aggregate pattern-level calls to gene-level calls."""

    print(f"\nAGGREGATING TO GENE-LEVEL CALLS")
    print("=" * 40)
    print(f"Low-confidence threshold: >= {MIN_LOW_CONF_PATTERNS} patterns required")
    print()

    gene_calls = []

    for sample_id in cnv_df['sample_id'].unique():
        sample_data = cnv_df[cnv_df['sample_id'] == sample_id]

        for gene in sample_data['gene'].unique():
            gene_data  = sample_data[sample_data['gene'] == gene]
            reliable   = gene_data[gene_data['qc_flag'] == 'pass']
            n_excluded = len(gene_data) - len(reliable)

            high_amps = len(reliable[(reliable['cnv_call'] == 'amplification') &
                                     (reliable['confidence'] == 'high')])
            med_amps  = len(reliable[(reliable['cnv_call'] == 'amplification') &
                                     (reliable['confidence'] == 'medium')])
            low_amps  = len(reliable[(reliable['cnv_call'] == 'amplification') &
                                     (reliable['confidence'] == 'low')])

            if high_amps >= 1:
                final_call, final_confidence = 1, 'high'
            elif med_amps >= 1:
                final_call, final_confidence = 1, 'medium'
            elif low_amps >= MIN_LOW_CONF_PATTERNS:
                final_call, final_confidence = 1, 'low'
            else:
                final_call, final_confidence = 0, 'high'

            gene_calls.append({
                'sample_id':          sample_id,
                'gene':               gene,
                'final_call':         final_call,
                'confidence':         final_confidence,
                'high_conf_patterns': high_amps,
                'med_conf_patterns':  med_amps,
                'low_conf_patterns':  low_amps,
                'total_patterns':     len(gene_data),
                'patterns_excluded':  n_excluded,
                'reliable_patterns':  len(reliable)
            })

    gene_calls_df = pd.DataFrame(gene_calls)

    print("Gene-level amplification rates:")
    for gene in gene_calls_df['gene'].unique():
        gd      = gene_calls_df[gene_calls_df['gene'] == gene]
        g_amps  = len(gd[gd['final_call'] == 1])
        g_total = len(gd)
        if g_total > 0:
            print(f"  {gene}: {g_amps}/{g_total} ({100*g_amps/g_total:.1f}%)")

    return gene_calls_df


def compare_against_linear(gene_calls_df, metadata_path, pf8_path=None):
    """Compare pangenome CNV calls against MalariaGen linear mapping calls.

    Uses concordance/discordance terminology rather than TP/TN/FP/FN, since
    MalariaGen linear calls serve as a reference comparator rather than a
    gold standard truth set. Discordant calls reflect methodological differences
    rather than definitive errors in either pipeline.

    Concordance categories:
      concordant_amplification     — both pipelines call amplification
      concordant_non_amplification — both pipelines call non-amplification
      discordant_pan_only          — pangenome calls amplification, linear does not
      discordant_linear_only       — linear calls amplification, pangenome does not

    When pf8_path is provided, concordance is additionally broken down by
    MalariaGen evidence type:
      Breakpoint-confirmed calls (highest confidence linear calls)
      Coverage-only calls (no breakpoint evidence in linear pipeline)
    """

    print(f"\nCOMPARING PANGENOME vs LINEAR MAPPING CNV CALLS")
    print("=" * 55)

    col_map = {
        'MDR1': {
            'final':      'MDR1_final_amplification_call',
            'coverage':   'MDR1_curated_coverage_only',
            'breakpoint': 'MDR1_breakpoint',
        },
        'plasmepsin': {
            'final':      'PM2_PM3_final_amplification_call',
            'coverage':   'PM2_PM3_curated_coverage_only',
            'breakpoint': 'PM2_PM3_breakpoint',
        },
    }

    try:
        metadata_df = pd.read_csv(metadata_path)
        print(f"Loaded metadata: {len(metadata_df)} samples")
    except Exception as e:
        print(f"Error loading metadata: {e}")
        return pd.DataFrame()

    pf8_df = None
    if pf8_path:
        try:
            pf8_df = pd.read_csv(pf8_path, sep='\t')
            print(f"Loaded Pf8 reference: {len(pf8_df)} samples")
        except Exception as e:
            print(f"Warning: could not load Pf8 file ({e}) — skipping split comparison")

    comparison_results = []

    for gene in gene_calls_df['gene'].unique():
        cols = col_map.get(gene)
        if not cols or cols['final'] not in metadata_df.columns:
            print(f"  No column mapping for {gene}")
            continue

        print(f"\nComparing {gene}")
        gene_data      = gene_calls_df[gene_calls_df['gene'] == gene]
        compared_count = 0

        for _, row in gene_data.iterrows():
            sample_id = row['sample_id']
            meta_row  = metadata_df[metadata_df['wgs_id'] == sample_id]
            if meta_row.empty:
                continue

            linear_call = meta_row[cols['final']].iloc[0]
            if pd.isna(linear_call) or linear_call == -1:
                continue

            pan_call      = row['final_call']
            linear_binary = 1 if linear_call == 1 else 0

            if pan_call == 1 and linear_binary == 1:
                concordance = 'concordant_amplification'
            elif pan_call == 0 and linear_binary == 0:
                concordance = 'concordant_non_amplification'
            elif pan_call == 1 and linear_binary == 0:
                concordance = 'discordant_pan_only'
            else:
                concordance = 'discordant_linear_only'

            record = {
                'gene':           gene,
                'sample_id':      sample_id,
                'sample_name':    meta_row['Sample'].iloc[0],
                'region':         meta_row.get('Region', ['Unknown']).iloc[0],
                'pan_call':       pan_call,
                'linear_call':    linear_binary,
                'confidence':     row['confidence'],
                'concordance':    concordance,
                'has_breakpoint': None,
                'coverage_only':  None,
            }

            if pf8_df is not None:
                pf8_row = pf8_df[pf8_df['Sample'] == meta_row['Sample'].iloc[0]]
                if not pf8_row.empty:
                    bp_val  = pf8_row[cols['breakpoint']].iloc[0]
                    cov_val = pf8_row[cols['coverage']].iloc[0]
                    record['has_breakpoint'] = (bp_val != '-') and not pd.isna(bp_val)
                    record['coverage_only']  = (cov_val == 1)

            comparison_results.append(record)
            compared_count += 1

        print(f"  Compared {compared_count} samples")

    comparison_df = pd.DataFrame(comparison_results)

    if len(comparison_df) == 0:
        return comparison_df

    _print_concordance(comparison_df, label="All calls (pangenome vs linear)")

    if pf8_df is not None and 'has_breakpoint' in comparison_df.columns:
        bp_subset  = comparison_df[comparison_df['has_breakpoint'] == True]
        cov_subset = comparison_df[
            (comparison_df['has_breakpoint'] == False) &
            (comparison_df['coverage_only']  == True)
        ]
        if len(bp_subset) > 0:
            _print_concordance(bp_subset,
                               label="vs Breakpoint-confirmed linear calls only")
        if len(cov_subset) > 0:
            _print_concordance(cov_subset,
                               label="vs Coverage-only linear calls (no breakpoint)")

    return comparison_df


def _print_concordance(df, label=""):
    """Print per-gene concordance metrics for a comparison subset."""

    print(f"\n  --- {label} (N={len(df)}) ---")
    for gene in df['gene'].unique():
        gd        = df[df['gene'] == gene]
        n         = len(gd)
        conc_amp  = len(gd[gd['concordance'] == 'concordant_amplification'])
        conc_non  = len(gd[gd['concordance'] == 'concordant_non_amplification'])
        disc_pan  = len(gd[gd['concordance'] == 'discordant_pan_only'])
        disc_lin  = len(gd[gd['concordance'] == 'discordant_linear_only'])

        lin_amps  = conc_amp + disc_lin
        pan_amps  = conc_amp + disc_pan
        overall   = round((conc_amp + conc_non) / n,     3) if n       > 0 else None
        det_rate  = round(conc_amp / lin_amps,            3) if lin_amps > 0 else None
        pos_conc  = round(conc_amp / pan_amps,            3) if pan_amps > 0 else None

        print(f"  {gene:12s}: N={n:4d}  "
              f"Conc.amp={conc_amp:3d}  Conc.non-amp={conc_non:4d}  "
              f"Disc.pan-only={disc_pan:3d}  Disc.linear-only={disc_lin:3d}  "
              f"Overall={overall}  Detection={det_rate}  Pos.concordance={pos_conc}")


def build_concordance_summary(comparison_df):
    """Build a concordance summary table across genes and confidence levels.

    Columns:
      Gene                     — target gene
      Confidence               — call confidence (All / High / Medium / Low)
      N                        — total samples compared
      Concordant_amp           — both pipelines call amplification
      Concordant_non_amp       — both pipelines call non-amplification
      Discordant_pan_only      — pangenome calls amplification, linear does not
      Discordant_linear_only   — linear calls amplification, pangenome does not
      Overall_concordance      — (concordant_amp + concordant_non_amp) / N
      Detection_rate           — concordant_amp / (concordant_amp + discordant_linear_only)
      Pos_concordance          — concordant_amp / (concordant_amp + discordant_pan_only)
      Neg_concordance          — concordant_non_amp / (concordant_non_amp + discordant_linear_only)
      F1                       — harmonic mean of detection rate and pos. concordance
    """

    rows = []

    for gene in comparison_df['gene'].unique():
        gd = comparison_df[comparison_df['gene'] == gene]

        for label, subset in [
            ('All',    gd),
            ('High',   gd[gd['confidence'] == 'high']),
            ('Medium', gd[gd['confidence'] == 'medium']),
            ('Low',    gd[gd['confidence'] == 'low']),
        ]:
            if len(subset) == 0:
                continue

            n         = len(subset)
            conc_amp  = len(subset[subset['concordance'] == 'concordant_amplification'])
            conc_non  = len(subset[subset['concordance'] == 'concordant_non_amplification'])
            disc_pan  = len(subset[subset['concordance'] == 'discordant_pan_only'])
            disc_lin  = len(subset[subset['concordance'] == 'discordant_linear_only'])

            lin_amps  = conc_amp + disc_lin
            pan_amps  = conc_amp + disc_pan
            non_amps  = conc_non + disc_lin

            overall_conc = round((conc_amp + conc_non) / n,        3) if n        > 0 else None
            det_rate     = round(conc_amp / lin_amps,               3) if lin_amps > 0 else None
            pos_conc     = round(conc_amp / pan_amps,               3) if pan_amps > 0 else None
            neg_conc     = round(conc_non / non_amps,               3) if non_amps > 0 else None
            f1           = round(2*conc_amp / (2*conc_amp + disc_pan + disc_lin), 3) \
                           if (2*conc_amp + disc_pan + disc_lin) > 0 else None

            rows.append({
                'Gene':                   gene,
                'Confidence':             label,
                'N':                      n,
                'Concordant_amp':         conc_amp,
                'Concordant_non_amp':     conc_non,
                'Discordant_pan_only':    disc_pan,
                'Discordant_linear_only': disc_lin,
                'Overall_concordance':    overall_conc,
                'Detection_rate':         det_rate,
                'Pos_concordance':        pos_conc,
                'Neg_concordance':        neg_conc,
                'F1':                     f1,
            })

    summary_df = pd.DataFrame(rows)

    print()
    print('=' * 110)
    print('CONCORDANCE SUMMARY (pangenome vs MalariaGen linear mapping)')
    print('=' * 110)
    print(summary_df[summary_df['Confidence'] == 'All'].to_string(index=False))
    print()
    print('By confidence level:')
    print(summary_df[summary_df['Confidence'] != 'All'].to_string(index=False))

    return summary_df


def main():
    parser = argparse.ArgumentParser(
        description='Flanking Region CNV Concordance — pangenome vs linear mapping'
    )
    parser.add_argument('--mosdepth_dir',
                        default='/mnt/storage13/nbillows/pangenome/shortreads/pan_files')
    parser.add_argument('--metadata',
                        default='/mnt/storage13/nbillows/pangenome/final_no_mix/analyse/cnv/meta_cnv_calls.csv')
    parser.add_argument('--max_files', type=int, default=None)
    parser.add_argument('--min_flank_coverage', type=float, default=MIN_FLANK_COVERAGE,
                        help=f'Minimum mean flank coverage to trust ratio (default: {MIN_FLANK_COVERAGE})')
    parser.add_argument('--max_flank_ratio', type=float, default=MAX_FLANK_RATIO,
                        help=f'Maximum left/right flank imbalance (default: {MAX_FLANK_RATIO})')
    parser.add_argument('--pf8',
                        default=None,
                        help='Path to Pf8 CNV calls TSV (enables breakpoint vs coverage split)')
    parser.add_argument('--min_low_conf_patterns', type=int, default=MIN_LOW_CONF_PATTERNS,
                        help=f'Low-confidence patterns required to call amplification (default: {MIN_LOW_CONF_PATTERNS})')

    args = parser.parse_args()

    import sys
    this = sys.modules[__name__]
    setattr(this, 'MIN_FLANK_COVERAGE',    args.min_flank_coverage)
    setattr(this, 'MAX_FLANK_RATIO',       args.max_flank_ratio)
    setattr(this, 'MIN_LOW_CONF_PATTERNS', args.min_low_conf_patterns)

    print("PANGENOME CNV CONCORDANCE — FLANKING REGION APPROACH")
    print("=" * 70)
    print(f"QC thresholds:")
    print(f"  Min flank coverage:       {MIN_FLANK_COVERAGE}x")
    print(f"  Max flank imbalance:      {MAX_FLANK_RATIO}x")
    print(f"  Min low-conf patterns:    {MIN_LOW_CONF_PATTERNS}")
    print(f"  MDR1_199kb pattern:       EXCLUDED (subtelomeric flank artefact)")
    print()

    patterns_df = define_breakpoints_with_flanks()
    print(f"Defined {len(patterns_df)} patterns:")
    for gene in patterns_df['gene'].unique():
        count = len(patterns_df[patterns_df['gene'] == gene])
        print(f"  {gene}: {count} patterns")

    # Handle multiple mosdepth output naming conventions
    mosdepth_files = []
    for pattern in ['*.regions.bed.gz', '*_regions_bed.gz',
                    '*_sort_bam_regions_bed.gz', '*_sort.bam.regions.bed.gz']:
        mosdepth_files.extend(glob.glob(os.path.join(args.mosdepth_dir, pattern)))
    mosdepth_files = sorted(set(mosdepth_files))

    if not mosdepth_files:
        print(f"No mosdepth files found in {args.mosdepth_dir}")
        return

    if args.max_files:
        mosdepth_files = mosdepth_files[:args.max_files]

    print(f"\nProcessing {len(mosdepth_files)} mosdepth files...")

    all_coverage_data = []
    for i, file_path in enumerate(mosdepth_files):
        if (i + 1) % 100 == 0:
            print(f"  Processed {i+1}/{len(mosdepth_files)} files...")
        file_data = process_mosdepth_file_with_flanks(file_path, patterns_df)
        if not file_data.empty:
            all_coverage_data.append(file_data)

    if not all_coverage_data:
        print("No coverage data extracted")
        return

    coverage_df   = pd.concat(all_coverage_data, ignore_index=True)
    print(f"Combined coverage data: {len(coverage_df):,} region-sample pairs")

    cnv_df        = call_cnvs_from_ratios(coverage_df)
    gene_calls_df = aggregate_to_gene_calls(cnv_df)
    comparison_df = compare_against_linear(
        gene_calls_df, args.metadata, pf8_path=args.pf8)

    if len(comparison_df) > 0:
        summary_df = build_concordance_summary(comparison_df)
    else:
        summary_df = pd.DataFrame()

    coverage_df.to_csv('flanking_coverage_data.csv',        index=False)
    cnv_df.to_csv('flanking_cnv_calls.csv',                 index=False)
    gene_calls_df.to_csv('flanking_gene_calls.csv',         index=False)
    comparison_df.to_csv('flanking_comparison_results.csv', index=False)
    if len(summary_df) > 0:
        summary_df.to_csv('concordance_summary.csv',        index=False)

    print("\nRESULTS SAVED:")
    print("  flanking_coverage_data.csv       — raw coverage per pattern per sample")
    print("  flanking_cnv_calls.csv           — pattern-level CNV calls")
    print("  flanking_gene_calls.csv          — gene-level aggregated calls")
    print("  flanking_comparison_results.csv  — per-sample concordance with linear calls")
    print("  concordance_summary.csv          — concordance summary table")


if __name__ == "__main__":
    main()