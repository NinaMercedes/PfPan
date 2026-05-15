#!/usr/bin/env python3
"""
Annotate BLAST results with gene information.
Run after BLAST steps are complete.

Usage:
    python3 annotate_svs.py --outdir sv_results --gff annotation.gff3 [--flank 2000]
"""
import os
import re
import sys
import argparse
import datetime
from collections import defaultdict

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--outdir",  required=True)
    p.add_argument("--gff",     required=True)
    p.add_argument("--flank",   type=int, default=2000,
                   help="bp window around BLAST hit for gene overlap (default 2000)")
    return p.parse_args()

# =============================================================================
# 1. Parse GFF into gene dict
# =============================================================================
def parse_gff(gff_path):
    genes_by_chrom = defaultdict(list)
    gene_info = {}

    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'gene':
                continue
            chrom  = fields[0]
            start  = int(fields[3]) - 1   # convert to 0-based
            end    = int(fields[4])
            attrs  = fields[8]

            gid   = re.search(r'ID=([^;]+)',         attrs)
            gname = re.search(r'Name=([^;]+)',        attrs)
            desc  = re.search(r'description=([^;]+)', attrs)

            gid   = gid.group(1).replace('gene:', '').strip()   if gid   else 'unknown'
            gname = gname.group(1).replace('gene:', '').strip()  if gname else gid
            desc  = desc.group(1)                                if desc  else ''

            gene_info[gid] = {'name': gname, 'desc': desc,
                              'chrom': chrom, 'start': start, 'end': end}
            genes_by_chrom[chrom].append((start, end, gid))

    print(f"  Loaded {len(gene_info)} genes from GFF")
    return genes_by_chrom, gene_info

# =============================================================================
# 2. Write gene BED
# =============================================================================
def write_gene_bed(genes_by_chrom, gene_info, out_path):
    with open(out_path, 'w') as f:
        for chrom, entries in genes_by_chrom.items():
            for start, end, gid in entries:
                name = gene_info[gid]['name']
                desc = gene_info[gid]['desc']
                f.write(f"{chrom}\t{start}\t{end}\t{gid}\t{name}\t{desc}\n")
    print(f"  Gene BED written: {out_path}")

# =============================================================================
# 3. Find overlapping genes for a BLAST hit
# =============================================================================
def find_genes(chrom, h_start, h_end, genes_by_chrom, gene_info, flank):
    """
    Find genes overlapping [h_start, h_end] on chrom.
    Uses flank window to catch promoter/downstream hits.
    Genes found only in flank zone are tagged with ~.
    Returns (gene_string, is_genic)
    """
    hits = []
    for gs, ge, gid in genes_by_chrom.get(chrom, []):
        # Check overlap with flank
        if gs <= h_end + flank and ge >= h_start - flank:
            true_ov   = max(0, min(h_end, ge) - max(h_start, gs))
            gene_len  = ge - gs
            pct       = true_ov / gene_len * 100 if gene_len > 0 else 0
            in_gene   = gs <= h_end and ge >= h_start
            tag       = '' if in_gene else '~'
            name      = gene_info[gid]['name']
            hits.append((pct, in_gene, f"{tag}{gid}({name},{pct:.0f}%)"))

    if not hits:
        return 'intergenic', False

    hits.sort(key=lambda x: (x[1], x[0]), reverse=True)
    gene_str = ';'.join(h[2] for h in hits[:5])
    is_genic = any(h[1] for h in hits)   # True if any hit is directly in a gene
    return gene_str, is_genic

# =============================================================================
# 4. Annotate BLAST results
# =============================================================================
def annotate_blast(blast_file, sv_type, out_file, genes_by_chrom, gene_info, flank):
    header = ("sv_id\tsv_type\tn_hits\tmax_bitscore\ttotal_bitscore\t"
              "score_ratio\tcopy_number_in_insert\thit_chrom\thit_start\thit_end\t"
              "pident\tqcovs\toverlapping_genes\tis_genic\tis_same_chrom\n")

    if not os.path.exists(blast_file) or os.path.getsize(blast_file) == 0:
        print(f"  No BLAST results for {sv_type}")
        with open(out_file, 'w') as f:
            f.write(header)
        return []

    # Group hits by query
    query_hits = defaultdict(list)
    with open(blast_file) as f:
        for line in f:
            ff = line.strip().split('\t')
            if len(ff) < 13:
                continue
            q, s, pi, ln, ql, sl, qc, qs, qe, ss, se, ev, bs = ff
            query_hits[q].append(dict(
                sseqid  = s,
                pident  = float(pi),
                qcovs   = float(qc),
                sstart  = int(ss),
                send    = int(se),
                bitscore= float(bs)
            ))

    results = []
    for qid, hits in query_hits.items():
        hits.sort(key=lambda x: x['bitscore'], reverse=True)
        max_s = hits[0]['bitscore']
        tot_s = sum(h['bitscore'] for h in hits)
        ratio = tot_s / max_s if max_s > 0 else 1.0
        cn    = round(ratio)
        best  = hits[0]

        # Handle reverse-strand hits (sstart > send)
        h_start = min(best['sstart'], best['send'])
        h_end   = max(best['sstart'], best['send'])
        h_chrom = best['sseqid']

        genes, is_genic = find_genes(h_chrom, h_start, h_end,
                                     genes_by_chrom, gene_info, flank)

        # Extract SV chromosome from ID: Pf3D7_12_v3_968881_INS10418
        # chromosome is everything before the last numeric_SVtype block
        # Pattern: chrom ends at last underscore before position digits
        parts = qid.split('_')
        # find where the position number starts
        sv_chrom = h_chrom   # fallback
        for i, p in enumerate(parts):
            if p.isdigit():
                sv_chrom = '_'.join(parts[:i])
                break

        same_chrom = (h_chrom == sv_chrom)

        results.append(dict(
            sv_id      = qid,
            sv_type    = sv_type,
            n_hits     = len(hits),
            max_bs     = max_s,
            tot_bs     = tot_s,
            ratio      = ratio,
            cn         = cn,
            hit_chrom  = h_chrom,
            hit_start  = h_start,
            hit_end    = h_end,
            pident     = best['pident'],
            qcovs      = best['qcovs'],
            genes      = genes,
            is_genic   = is_genic,
            same_chrom = same_chrom,
        ))

    with open(out_file, 'w') as f:
        f.write(header)
        for r in sorted(results, key=lambda x: x['tot_bs'], reverse=True):
            f.write(f"{r['sv_id']}\t{r['sv_type']}\t{r['n_hits']}\t"
                    f"{r['max_bs']:.1f}\t{r['tot_bs']:.1f}\t{r['ratio']:.2f}\t{r['cn']}\t"
                    f"{r['hit_chrom']}\t{r['hit_start']}\t{r['hit_end']}\t"
                    f"{r['pident']:.1f}\t{r['qcovs']:.1f}\t"
                    f"{r['genes']}\t{r['is_genic']}\t{r['same_chrom']}\n")

    g = sum(1 for r in results if r['is_genic'])
    t = sum(1 for r in results if r['is_genic'] and r['same_chrom'])
    print(f"  {sv_type}: {len(results)} SVs | {g} genic | {t} same-chrom (tandem dup candidates)")
    return results

# =============================================================================
# 5. Whole/partial gene deletions from bedtools overlap
# =============================================================================
def detect_gene_deletions(overlap_file, gene_info):
    whole, partial = [], []

    if not os.path.exists(overlap_file) or os.path.getsize(overlap_file) == 0:
        return whole, partial

    with open(overlap_file) as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 10:
                continue
            sv_chrom, sv_start, sv_end, sv_id = fields[:4]
            g_chrom, g_start, g_end, g_id, g_name, g_desc = fields[4:10]

            g_id   = g_id.replace('gene:', '').strip()
            g_name = g_name.replace('gene:', '').strip()

            sv_start, sv_end = int(sv_start), int(sv_end)
            g_start,  g_end  = int(g_start),  int(g_end)
            gene_len = g_end - g_start
            overlap  = max(0, min(sv_end, g_end) - max(sv_start, g_start))
            ovlp_pct = overlap / gene_len * 100 if gene_len > 0 else 0

            rec = dict(
                sv_id=sv_id, sv_chrom=sv_chrom,
                sv_start=sv_start, sv_end=sv_end,
                sv_len=sv_end - sv_start,
                gene_id=g_id, gene_name=g_name,
                gene_desc=g_desc, gene_len=gene_len,
                overlap_bp=overlap, overlap_pct=ovlp_pct
            )
            if ovlp_pct >= 90:
                whole.append(rec)
            elif ovlp_pct >= 20:
                partial.append(rec)

    return whole, partial

# =============================================================================
# 6. Summary report
# =============================================================================
RESISTANCE_GENES = {
    'PF3D7_1224000': 'GCH1 (antifolate)',
    'PF3D7_0523000': 'MDR1 (multidrug)',
    'PF3D7_0810800': 'DHPS (antifolate)',
    'PF3D7_0417200': 'DHFR (antifolate)',
    'PF3D7_1343700': 'Kelch13 (artemisinin)',
    'PF3D7_1408100': 'Plasmepsin III (piperaquine)',
    'PF3D7_1408000': 'Plasmepsin II (piperaquine)',
    'PF3D7_0709000': 'CRT (chloroquine)',
}

def flag_res(s):
    for gid, drug in RESISTANCE_GENES.items():
        if gid in s:
            return drug
    return ''

def write_summary(outdir, flank, ins_results, del_results, whole, partial):
    lines = [
        "=" * 65,
        "PANGENOME SV PIPELINE - SUMMARY REPORT",
        f"Output:    {outdir}",
        f"Generated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}",
        f"Gene flank window: {flank}bp (~ = promoter/downstream)",
        "=" * 65,
    ]

    for label, results in [("INSERTIONS", ins_results), ("DELETIONS", del_results)]:
        genic  = [r for r in results if r['is_genic']]
        tandem = [r for r in results if r['is_genic'] and r['same_chrom']]
        res    = [r for r in results if flag_res(r['genes'])]
        lines += [
            f"\n{label}",
            f"  Total with BLAST hits:        {len(results)}",
            f"  Hitting annotated genes:      {len(genic)}",
            f"  Putative tandem duplications: {len(tandem)}",
            f"  Resistance gene hits:         {len(res)}",
        ]
        if res:
            lines.append(f"\n  *** RESISTANCE GENE {label} ***")
            for r in res:
                drug = flag_res(r['genes'])
                lines.append(f"    {r['sv_id']:<45} [{drug}]  cn={r['cn']}")
                lines.append(f"      genes: {r['genes'][:80]}")

    lines += [
        f"\nWHOLE GENE DELETIONS (>90% gene body covered)",
        f"  Total:                        {len(whole)}",
        f"  Resistance gene hits:         {sum(1 for r in whole if flag_res(r['gene_id']))}",
    ]
    if whole:
        lines.append(f"\n  {'SV_ID':<40} {'GENE':<20} {'OVERLAP%':>9}  FLAG")
        for r in sorted(whole, key=lambda x: x['overlap_pct'], reverse=True)[:25]:
            tag  = f"  *** {flag_res(r['gene_id'])}" if flag_res(r['gene_id']) else ''
            lines.append(f"  {r['sv_id']:<40} {r['gene_name']:<20} {r['overlap_pct']:>9.1f}{tag}")

    lines += [
        f"\nPARTIAL GENE DELETIONS (20-90% gene body covered)",
        f"  Total:                        {len(partial)}",
    ]

    lines += [
        f"\nOUTPUT FILES",
        f"  {outdir}/results/insertions_annotated.tsv",
        f"  {outdir}/results/deletions_annotated.tsv",
        f"  {outdir}/results/whole_gene_deletions.tsv",
        f"  {outdir}/results/partial_gene_deletions.tsv",
        f"  {outdir}/results/summary_report.txt",
        "=" * 65,
    ]

    report = '\n'.join(lines)
    print(report)
    with open(f"{outdir}/results/summary_report.txt", 'w') as f:
        f.write(report)

# =============================================================================
# Main
# =============================================================================
def main():
    args = parse_args()
    outdir = args.outdir
    flank  = args.flank

    os.makedirs(f"{outdir}/results", exist_ok=True)
    os.makedirs(f"{outdir}/tmp",     exist_ok=True)

    print(f"\n[Step 1] Parsing GFF: {args.gff}")
    genes_by_chrom, gene_info = parse_gff(args.gff)
    write_gene_bed(genes_by_chrom, gene_info, f"{outdir}/tmp/genes.bed")

    print(f"\n[Step 2] Annotating insertions...")
    ins_results = annotate_blast(
        f"{outdir}/blast/insertions_blast_raw.txt", "INS",
        f"{outdir}/results/insertions_annotated.tsv",
        genes_by_chrom, gene_info, flank
    )

    print(f"\n[Step 3] Annotating deletions...")
    del_results = annotate_blast(
        f"{outdir}/blast/deletions_blast_raw.txt", "DEL",
        f"{outdir}/results/deletions_annotated.tsv",
        genes_by_chrom, gene_info, flank
    )

    print(f"\n[Step 4] Detecting whole/partial gene deletions...")
    header = ("sv_id\tsv_chrom\tsv_start\tsv_end\tsv_len\t"
              "gene_id\tgene_name\tgene_desc\tgene_len\toverlap_bp\toverlap_pct\n")

    whole, partial = detect_gene_deletions(
        f"{outdir}/tmp/deletions_gene_overlap.bed", gene_info
    )
    print(f"  Whole   (>90%): {len(whole)}")
    print(f"  Partial (>20%): {len(partial)}")

    for fname, data in [
        (f"{outdir}/results/whole_gene_deletions.tsv",   whole),
        (f"{outdir}/results/partial_gene_deletions.tsv", partial),
    ]:
        with open(fname, 'w') as f:
            f.write(header)
            for r in sorted(data, key=lambda x: x['overlap_pct'], reverse=True):
                f.write(f"{r['sv_id']}\t{r['sv_chrom']}\t{r['sv_start']}\t{r['sv_end']}\t"
                        f"{r['sv_len']}\t{r['gene_id']}\t{r['gene_name']}\t{r['gene_desc']}\t"
                        f"{r['gene_len']}\t{r['overlap_bp']}\t{r['overlap_pct']:.1f}\n")

    print(f"\n[Step 5] Writing summary report...")
    write_summary(outdir, flank, ins_results, del_results, whole, partial)

if __name__ == "__main__":
    main()
