#!/usr/bin/env bash
# =============================================================================
# Pangenome SV Pipeline - vg/ODGI format
# Designed for pangenome VCFs with AT allele traversal tags (no SVTYPE/SVLEN).
# SVs are identified by REF/ALT length difference.
# Filters to core genome, extracts insertion/deletion sequences via MAFFT,
# BLASTs against reference, intersects with gene annotations,
# detects whole gene deletions, generates summary report.
#
# Usage:
#   bash pangenome_sv_pipeline.sh \
#     -v pangenome.vcf.gz \
#     -b core_genome.bed \
#     -r reference.fasta \
#     -g annotation.gff3 \
#     -o output_dir \
#     [-t threads] [-m min_svlen] [-p min_pident] [-q min_qcov] [-s from_step]
#
# Steps:
#   1  Filter VCF to core genome
#   2  Classify SVs by REF/ALT length
#   3  Build BLAST database
#   4  Extract insertion sequences (MAFFT)
#   5  Extract deletion sequences (MAFFT)
#   6  BLAST insertions
#   7  BLAST deletions
#   8  Parse GFF to gene BED
#   9  Intersect SVs with genes (bedtools)
#   10 Detect whole/partial gene deletions
#   11 Annotate BLAST results
#   12 Summary report
#
# To restart from after BLAST (skip steps 1-7):
#   bash pangenome_sv_pipeline.sh [all args] -s 8
#
# Dependencies: bcftools, mafft, blastn, makeblastdb, bedtools, python3
# =============================================================================

set -euo pipefail

# --- Defaults ----------------------------------------------------------------
THREADS=8
MIN_SVLEN=2000
MIN_PIDENT=85
MIN_QCOV=20
FROM_STEP=1
FLANK=2000   # bp window around BLAST hit for gene overlap (captures promoter hits)

# --- Argument parsing --------------------------------------------------------
usage() {
    sed -n '/^# Usage/,/^# Steps/p' "$0"
    exit 1
}

while getopts "v:b:r:g:o:t:m:p:q:s:" opt; do
    case $opt in
        v) VCF="$OPTARG"       ;;
        b) BED="$OPTARG"       ;;
        r) REF="$OPTARG"       ;;
        g) GFF="$OPTARG"       ;;
        o) OUTDIR="$OPTARG"    ;;
        t) THREADS="$OPTARG"   ;;
        m) MIN_SVLEN="$OPTARG" ;;
        p) MIN_PIDENT="$OPTARG";;
        q) MIN_QCOV="$OPTARG"  ;;
        s) FROM_STEP="$OPTARG" ;;
        *) usage ;;
    esac
done

# When restarting from a later step, VCF/BED/REF are not required
if [ "$FROM_STEP" -le 1 ]; then
    for var in VCF BED REF GFF OUTDIR; do
        if [ -z "${!var:-}" ]; then
            echo "ERROR: Missing required argument for $var"
            usage
        fi
    done
else
    for var in GFF OUTDIR; do
        if [ -z "${!var:-}" ]; then
            echo "ERROR: Missing required argument for $var"
            usage
        fi
    done
fi

export VCF=${VCF:-""} BED=${BED:-""} REF=${REF:-""} GFF OUTDIR
export THREADS MIN_SVLEN MIN_PIDENT MIN_QCOV FROM_STEP FLANK

# --- Helpers -----------------------------------------------------------------
should_run() {
    # Returns 0 (true) if the step number >= FROM_STEP
    [ "$1" -ge "$FROM_STEP" ]
}

# --- Setup -------------------------------------------------------------------
mkdir -p "$OUTDIR"/{logs,vcf,sequences,blast,results,tmp}
LOG="$OUTDIR/logs/pipeline.log"
exec > >(tee -a "$LOG") 2>&1

echo "============================================="
echo "Pangenome SV Pipeline (vg/ODGI format)"
echo "Started:   $(date)"
echo "GFF:       $GFF"
echo "OUTDIR:    $OUTDIR"
echo "MIN_SVLEN: $MIN_SVLEN"
echo "FROM_STEP: $FROM_STEP"
echo "FLANK:     ${FLANK}bp gene window"
if [ "$FROM_STEP" -le 3 ]; then
    echo "VCF:       $VCF"
    echo "BED:       $BED"
    echo "REF:       $REF"
fi
echo "============================================="

# =============================================================================
# STEP 1: Filter VCF to core genome
# =============================================================================
if should_run 1; then
    echo ""
    echo "[STEP 1] Filtering VCF to core genome..."
    bcftools view -R "$BED" "$VCF" -O z -o "$OUTDIR/vcf/core_genome.vcf.gz"
    bcftools index "$OUTDIR/vcf/core_genome.vcf.gz"
    TOTAL=$(bcftools view -H "$OUTDIR/vcf/core_genome.vcf.gz" | wc -l)
    echo "  Total records in core genome: $TOTAL"
else
    echo "[STEP 1] Skipping (FROM_STEP=$FROM_STEP)"
fi

# =============================================================================
# STEP 2: Classify SVs by REF/ALT length difference
# =============================================================================
if should_run 2; then
    echo ""
    echo "[STEP 2] Classifying SVs by REF/ALT length difference..."

    python3 << 'PYEOF'
import gzip, os, sys

outdir    = os.environ['OUTDIR']
min_svlen = int(os.environ['MIN_SVLEN'])
vcf_in    = f"{outdir}/vcf/core_genome.vcf.gz"
ins_tmp   = f"{outdir}/tmp/insertions.vcf"
del_tmp   = f"{outdir}/tmp/deletions.vcf"

n_ins = n_del = n_skip_short = n_skip_symbolic = 0
sizes_ins, sizes_del = [], []

with gzip.open(vcf_in, 'rt') as fin, \
     open(ins_tmp, 'w') as fins, \
     open(del_tmp, 'w') as fdel:

    for line in fin:
        if line.startswith('#'):
            fins.write(line)
            fdel.write(line)
            continue
        fields = line.strip().split('\t')
        if len(fields) < 5:
            continue
        ref, alt = fields[3], fields[4]
        if alt.startswith('<') or ',' in alt:
            n_skip_symbolic += 1
            continue
        diff = len(alt) - len(ref)
        if abs(diff) < min_svlen:
            n_skip_short += 1
            continue
        if diff > 0:
            fins.write(line)
            n_ins += 1
            sizes_ins.append(diff)
        else:
            fdel.write(line)
            n_del += 1
            sizes_del.append(abs(diff))

print(f"  Insertions (ALT-REF >= {min_svlen}bp): {n_ins}")
print(f"  Deletions  (REF-ALT >= {min_svlen}bp): {n_del}")
print(f"  Skipped (too short):                   {n_skip_short}")
print(f"  Skipped (symbolic/multi-allelic):       {n_skip_symbolic}")
print(f"\n  SV size distribution:")
for label, sizes in [("INS", sizes_ins), ("DEL", sizes_del)]:
    if not sizes:
        continue
    sizes.sort()
    n = len(sizes)
    print(f"    {label}: n={n} | min={sizes[0]} | "
          f"median={sizes[n//2]} | max={sizes[-1]} | "
          f">10kb={sum(1 for s in sizes if s > 10000)}")
PYEOF

    for svtype in insertions deletions; do
        if [ -s "$OUTDIR/tmp/${svtype}.vcf" ]; then
            bcftools view "$OUTDIR/tmp/${svtype}.vcf" \
                -O z -o "$OUTDIR/vcf/core_${svtype}.vcf.gz"
            bcftools index "$OUTDIR/vcf/core_${svtype}.vcf.gz"
        else
            bcftools view -h "$OUTDIR/vcf/core_genome.vcf.gz" \
                -O z -o "$OUTDIR/vcf/core_${svtype}.vcf.gz"
            bcftools index "$OUTDIR/vcf/core_${svtype}.vcf.gz"
        fi
    done
else
    echo "[STEP 2] Skipping (FROM_STEP=$FROM_STEP)"
fi

# =============================================================================
# STEP 3: Build BLAST database
# =============================================================================
if should_run 3; then
    echo ""
    echo "[STEP 3] Building BLAST database..."
    if [ ! -f "${REF}.nhr" ] && [ ! -f "${REF}.nin" ]; then
        makeblastdb -in "$REF" -dbtype nucl -out "$REF" -title "Pf_reference"
        echo "  BLAST database built"
    else
        echo "  BLAST database already exists, skipping"
    fi
else
    echo "[STEP 3] Skipping (FROM_STEP=$FROM_STEP)"
fi

# =============================================================================
# STEP 4: Extract insertion sequences using MAFFT
# =============================================================================
if should_run 4; then
    echo ""
    echo "[STEP 4] Extracting insertion sequences via MAFFT..."

    python3 << 'PYEOF'
import subprocess, os, gzip, tempfile, sys

outdir    = os.environ['OUTDIR']
vcf_path  = f"{outdir}/vcf/core_insertions.vcf.gz"
min_svlen = int(os.environ['MIN_SVLEN'])

def open_vcf(p):
    return gzip.open(p, 'rt') if p.endswith('.gz') else open(p, 'r')

def run_mafft(ref_seq, alt_seq, tmpdir):
    in_fa = f"{tmpdir}/input.fa"
    with open(in_fa, 'w') as f:
        f.write(f">ref\n{ref_seq}\n>alt\n{alt_seq}\n")
    r = subprocess.run(
        ["mafft", "--auto", "--quiet", "--thread", "1", in_fa],
        capture_output=True, text=True)
    if r.returncode != 0 or not r.stdout.strip():
        return None, None
    seqs, cur = {}, None
    for line in r.stdout.strip().split('\n'):
        if line.startswith('>'):
            cur = line[1:].split()[0]; seqs[cur] = []
        elif cur:
            seqs[cur].append(line.strip())
    aligned = {k: ''.join(v).upper() for k, v in seqs.items()}
    return aligned.get('ref'), aligned.get('alt')

def gap_blocks(seq_gaps, other, min_len):
    blocks, i = [], 0
    while i < len(seq_gaps):
        if seq_gaps[i] == '-':
            j = i
            while j < len(seq_gaps) and seq_gaps[j] == '-':
                j += 1
            ext = other[i:j].replace('-', '')
            if len(ext) >= min_len:
                blocks.append((i, j, ext, len(ext)))
            i = j
        else:
            i += 1
    return blocks

records = []
with open_vcf(vcf_path) as f:
    for line in f:
        if line.startswith('#'): continue
        fields = line.strip().split('\t')
        if len(fields) < 5: continue
        chrom, pos, ref, alt = fields[0], int(fields[1]), fields[3], fields[4]
        if alt.startswith('<') or ',' in alt: continue
        svlen = len(alt) - len(ref)
        if svlen < min_svlen: continue
        records.append({'id': f"{chrom}_{pos}_INS{svlen}",
                        'chrom': chrom, 'pos': pos,
                        'ref': ref, 'alt': alt, 'svlen': svlen})

print(f"  Processing {len(records)} INS records...")
if not records:
    print("  No records to process.")
    sys.exit(0)

extracted, failed = [], []
with tempfile.TemporaryDirectory() as tmpdir:
    for i, rec in enumerate(records):
        aln_ref, aln_alt = run_mafft(rec['ref'], rec['alt'], tmpdir)
        if aln_ref is None:
            failed.append(f"{rec['id']} (mafft_failed)"); continue
        blocks = gap_blocks(aln_ref, aln_alt, min_svlen)
        if not blocks:
            all_b = gap_blocks(aln_ref, aln_alt, 1)
            lg = max(all_b, key=lambda x: x[3])[3] if all_b else 0
            failed.append(f"{rec['id']} (largest_gap={lg}bp)"); continue
        blocks.sort(key=lambda x: x[3], reverse=True)
        aln_s, aln_e, sv_seq, sv_len = blocks[0]
        extracted.append({'id': rec['id'], 'chrom': rec['chrom'],
                          'pos': rec['pos'], 'svlen': rec['svlen'],
                          'aln_pos': aln_s, 'ext_len': sv_len, 'seq': sv_seq})
        if (i+1) % 10 == 0:
            print(f"  Processed {i+1}/{len(records)}...", flush=True)

with open(f"{outdir}/sequences/insertion_sequences.fa", 'w') as f:
    for r in extracted:
        f.write(f">{r['id']}\n{r['seq']}\n")
with open(f"{outdir}/sequences/insertion_summary.tsv", 'w') as f:
    f.write("sv_id\tchrom\tpos\tsvlen\taln_pos\textracted_len\n")
    for r in extracted:
        f.write(f"{r['id']}\t{r['chrom']}\t{r['pos']}\t{r['svlen']}\t{r['aln_pos']}\t{r['ext_len']}\n")
if failed:
    with open(f"{outdir}/logs/insertion_extraction_failed.txt", 'w') as f:
        f.write('\n'.join(failed) + '\n')
print(f"  Successfully extracted: {len(extracted)}")
print(f"  Failed:                 {len(failed)}")
PYEOF

else
    echo "[STEP 4] Skipping (FROM_STEP=$FROM_STEP)"
fi

# =============================================================================
# STEP 5: Extract deletion sequences using MAFFT
# =============================================================================
if should_run 5; then
    echo ""
    echo "[STEP 5] Extracting deletion sequences via MAFFT..."

    python3 << 'PYEOF'
import subprocess, os, gzip, tempfile, sys

outdir    = os.environ['OUTDIR']
vcf_path  = f"{outdir}/vcf/core_deletions.vcf.gz"
min_svlen = int(os.environ['MIN_SVLEN'])

def open_vcf(p):
    return gzip.open(p, 'rt') if p.endswith('.gz') else open(p, 'r')

def run_mafft(ref_seq, alt_seq, tmpdir):
    in_fa = f"{tmpdir}/input.fa"
    with open(in_fa, 'w') as f:
        f.write(f">ref\n{ref_seq}\n>alt\n{alt_seq}\n")
    r = subprocess.run(
        ["mafft", "--auto", "--quiet", "--thread", "1", in_fa],
        capture_output=True, text=True)
    if r.returncode != 0 or not r.stdout.strip():
        return None, None
    seqs, cur = {}, None
    for line in r.stdout.strip().split('\n'):
        if line.startswith('>'):
            cur = line[1:].split()[0]; seqs[cur] = []
        elif cur:
            seqs[cur].append(line.strip())
    aligned = {k: ''.join(v).upper() for k, v in seqs.items()}
    return aligned.get('ref'), aligned.get('alt')

def gap_blocks(seq_gaps, other, min_len):
    blocks, i = [], 0
    while i < len(seq_gaps):
        if seq_gaps[i] == '-':
            j = i
            while j < len(seq_gaps) and seq_gaps[j] == '-':
                j += 1
            ext = other[i:j].replace('-', '')
            if len(ext) >= min_len:
                blocks.append((i, j, ext, len(ext)))
            i = j
        else:
            i += 1
    return blocks

records = []
with open_vcf(vcf_path) as f:
    for line in f:
        if line.startswith('#'): continue
        fields = line.strip().split('\t')
        if len(fields) < 5: continue
        chrom, pos, ref, alt = fields[0], int(fields[1]), fields[3], fields[4]
        if alt.startswith('<') or ',' in alt: continue
        svlen = len(ref) - len(alt)
        if svlen < min_svlen: continue
        records.append({'id': f"{chrom}_{pos}_DEL{svlen}",
                        'chrom': chrom, 'pos': pos,
                        'ref': ref, 'alt': alt, 'svlen': svlen})

print(f"  Processing {len(records)} DEL records...")
if not records:
    print("  No records to process.")
    sys.exit(0)

extracted, failed = [], []
with tempfile.TemporaryDirectory() as tmpdir:
    for i, rec in enumerate(records):
        aln_ref, aln_alt = run_mafft(rec['ref'], rec['alt'], tmpdir)
        if aln_ref is None:
            failed.append(f"{rec['id']} (mafft_failed)"); continue
        blocks = gap_blocks(aln_alt, aln_ref, min_svlen)
        if not blocks:
            all_b = gap_blocks(aln_alt, aln_ref, 1)
            lg = max(all_b, key=lambda x: x[3])[3] if all_b else 0
            failed.append(f"{rec['id']} (largest_gap={lg}bp)"); continue
        blocks.sort(key=lambda x: x[3], reverse=True)
        aln_s, aln_e, sv_seq, sv_len = blocks[0]
        extracted.append({'id': rec['id'], 'chrom': rec['chrom'],
                          'pos': rec['pos'], 'svlen': rec['svlen'],
                          'aln_pos': aln_s, 'ext_len': sv_len, 'seq': sv_seq})
        if (i+1) % 10 == 0:
            print(f"  Processed {i+1}/{len(records)}...", flush=True)

with open(f"{outdir}/sequences/deletion_sequences.fa", 'w') as f:
    for r in extracted:
        f.write(f">{r['id']}\n{r['seq']}\n")
with open(f"{outdir}/sequences/deletion_summary.tsv", 'w') as f:
    f.write("sv_id\tchrom\tpos\tsvlen\taln_pos\textracted_len\n")
    for r in extracted:
        f.write(f"{r['id']}\t{r['chrom']}\t{r['pos']}\t{r['svlen']}\t{r['aln_pos']}\t{r['ext_len']}\n")
if failed:
    with open(f"{outdir}/logs/deletion_extraction_failed.txt", 'w') as f:
        f.write('\n'.join(failed) + '\n')
print(f"  Successfully extracted: {len(extracted)}")
print(f"  Failed:                 {len(failed)}")
PYEOF

else
    echo "[STEP 5] Skipping (FROM_STEP=$FROM_STEP)"
fi

# =============================================================================
# STEP 6: BLAST insertions
# =============================================================================
if should_run 6; then
    echo ""
    echo "[STEP 6] BLASTing insertion sequences against reference..."
    if [ -s "$OUTDIR/sequences/insertion_sequences.fa" ]; then
        blastn \
            -query         "$OUTDIR/sequences/insertion_sequences.fa" \
            -db            "$REF" \
            -outfmt        "6 qseqid sseqid pident length qlen slen qcovs qstart qend sstart send evalue bitscore" \
            -perc_identity $MIN_PIDENT \
            -qcov_hsp_perc $MIN_QCOV \
            -num_threads   $THREADS \
            -max_target_seqs 5 \
            -max_hsps      10 \
            -out "$OUTDIR/blast/insertions_blast_raw.txt"
        echo "  Hits: $(wc -l < "$OUTDIR/blast/insertions_blast_raw.txt")"
    else
        echo "  No insertion sequences to BLAST"
        touch "$OUTDIR/blast/insertions_blast_raw.txt"
    fi
else
    echo "[STEP 6] Skipping (FROM_STEP=$FROM_STEP)"
fi

# =============================================================================
# STEP 7: BLAST deletions
# =============================================================================
if should_run 7; then
    echo ""
    echo "[STEP 7] BLASTing deletion sequences against reference..."
    if [ -s "$OUTDIR/sequences/deletion_sequences.fa" ]; then
        blastn \
            -query         "$OUTDIR/sequences/deletion_sequences.fa" \
            -db            "$REF" \
            -outfmt        "6 qseqid sseqid pident length qlen slen qcovs qstart qend sstart send evalue bitscore" \
            -perc_identity $MIN_PIDENT \
            -qcov_hsp_perc $MIN_QCOV \
            -num_threads   $THREADS \
            -max_target_seqs 5 \
            -max_hsps      10 \
            -out "$OUTDIR/blast/deletions_blast_raw.txt"
        echo "  Hits: $(wc -l < "$OUTDIR/blast/deletions_blast_raw.txt")"
    else
        echo "  No deletion sequences to BLAST"
        touch "$OUTDIR/blast/deletions_blast_raw.txt"
    fi
else
    echo "[STEP 7] Skipping (FROM_STEP=$FROM_STEP)"
fi

# =============================================================================
# STEP 8: Parse GFF into gene BED
# =============================================================================
if should_run 8; then
    echo ""
    echo "[STEP 8] Parsing gene annotations from GFF..."

    python3 << 'PYEOF'
import os, re

gff    = os.environ['GFF']
outdir = os.environ['OUTDIR']

genes = []
with open(gff) as f:
    for line in f:
        if line.startswith('#'): continue
        fields = line.strip().split('\t')
        if len(fields) < 9 or fields[2] != 'gene': continue
        chrom, start, end, attrs = fields[0], int(fields[3]), int(fields[4]), fields[8]
        gid   = re.search(r'ID=([^;]+)',         attrs)
        gname = re.search(r'Name=([^;]+)',        attrs)
        desc  = re.search(r'description=([^;]+)', attrs)
        # Strip "gene:" prefix if present e.g. "gene:PF3D7_1224000" -> "PF3D7_1224000"
        gid   = gid.group(1).replace('gene:', '').strip()   if gid   else 'unknown'
        gname = gname.group(1).replace('gene:', '').strip()  if gname else gid
        desc  = desc.group(1)                                if desc  else ''
        genes.append((chrom, start-1, end, gid, gname, desc))

with open(f"{outdir}/tmp/genes.bed", 'w') as f:
    for g in genes:
        f.write('\t'.join(str(x) for x in g) + '\n')

print(f"  Extracted {len(genes)} genes from GFF")
print(f"  Sample entries:")
for g in genes[:3]:
    print(f"    {g[0]}:{g[1]}-{g[2]}  {g[3]}  {g[4]}")
PYEOF

else
    echo "[STEP 8] Skipping (FROM_STEP=$FROM_STEP)"
fi

# =============================================================================
# STEP 9: Make SV BED files and intersect with genes
# =============================================================================
if should_run 9; then
    echo ""
    echo "[STEP 9] Intersecting SVs with gene coordinates..."

    python3 << 'PYEOF'
import gzip, os

outdir    = os.environ['OUTDIR']
min_svlen = int(os.environ['MIN_SVLEN'])

for svtype, vcf_path, out_bed in [
    ("INS", f"{outdir}/vcf/core_insertions.vcf.gz", f"{outdir}/tmp/insertions.bed"),
    ("DEL", f"{outdir}/vcf/core_deletions.vcf.gz",  f"{outdir}/tmp/deletions.bed"),
]:
    n = 0
    with gzip.open(vcf_path, 'rt') as fin, open(out_bed, 'w') as fout:
        for line in fin:
            if line.startswith('#'): continue
            fields = line.strip().split('\t')
            if len(fields) < 5: continue
            chrom, pos, ref, alt = fields[0], int(fields[1]), fields[3], fields[4]
            if alt.startswith('<') or ',' in alt: continue
            if svtype == "INS":
                svlen = len(alt) - len(ref)
                end   = pos + svlen
            else:
                svlen = len(ref) - len(alt)
                end   = pos + svlen
            if svlen < min_svlen: continue
            sv_id = f"{chrom}_{pos}_{svtype}{svlen}"
            fout.write(f"{chrom}\t{pos-1}\t{end}\t{sv_id}\n")
            n += 1
    print(f"  {svtype} BED: {n} records written")
PYEOF

    bedtools intersect \
        -a "$OUTDIR/tmp/insertions.bed" \
        -b "$OUTDIR/tmp/genes.bed" \
        -wa -wb > "$OUTDIR/tmp/insertions_gene_overlap.bed" || true

    bedtools intersect \
        -a "$OUTDIR/tmp/deletions.bed" \
        -b "$OUTDIR/tmp/genes.bed" \
        -wa -wb > "$OUTDIR/tmp/deletions_gene_overlap.bed" || true

    echo "  Insertion-gene overlaps: $(wc -l < "$OUTDIR/tmp/insertions_gene_overlap.bed")"
    echo "  Deletion-gene overlaps:  $(wc -l < "$OUTDIR/tmp/deletions_gene_overlap.bed")"
else
    echo "[STEP 9] Skipping (FROM_STEP=$FROM_STEP)"
fi

# =============================================================================
# STEP 10: Whole and partial gene deletions
# =============================================================================
if should_run 10; then
    echo ""
    echo "[STEP 10] Detecting whole gene deletions..."

    python3 << 'PYEOF'
import os

outdir = os.environ['OUTDIR']
whole, partial = [], []

overlap_file = f"{outdir}/tmp/deletions_gene_overlap.bed"
if os.path.exists(overlap_file) and os.path.getsize(overlap_file) > 0:
    with open(overlap_file) as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 10: continue
            sv_chrom, sv_start, sv_end, sv_id = fields[:4]
            g_chrom, g_start, g_end, g_id, g_name, g_desc = fields[4:10]
            # Strip gene: prefix
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
                gene_id=g_id, gene_name=g_name, gene_desc=g_desc,
                gene_len=gene_len, overlap_bp=overlap, overlap_pct=ovlp_pct
            )
            if ovlp_pct >= 90:   whole.append(rec)
            elif ovlp_pct >= 20: partial.append(rec)
else:
    print("  No deletion-gene overlaps found")

header = ("sv_id\tsv_chrom\tsv_start\tsv_end\tsv_len\t"
          "gene_id\tgene_name\tgene_desc\tgene_len\toverlap_bp\toverlap_pct\n")

for fname, data, label in [
    (f"{outdir}/results/whole_gene_deletions.tsv",   whole,   "Whole   (>90%)"),
    (f"{outdir}/results/partial_gene_deletions.tsv", partial, "Partial (>20%)")
]:
    with open(fname, 'w') as f:
        f.write(header)
        for r in sorted(data, key=lambda x: x['overlap_pct'], reverse=True):
            f.write(f"{r['sv_id']}\t{r['sv_chrom']}\t{r['sv_start']}\t{r['sv_end']}\t"
                    f"{r['sv_len']}\t{r['gene_id']}\t{r['gene_name']}\t{r['gene_desc']}\t"
                    f"{r['gene_len']}\t{r['overlap_bp']}\t{r['overlap_pct']:.1f}\n")
    print(f"  {label}: {len(data)}")
PYEOF

else
    echo "[STEP 10] Skipping (FROM_STEP=$FROM_STEP)"
fi

# =============================================================================
# STEP 11: Annotate BLAST results with gene info and copy number
# =============================================================================
if should_run 11; then
    echo ""
    echo "[STEP 11] Annotating BLAST results..."

    python3 << 'PYEOF'
import os
from collections import defaultdict

outdir = os.environ['OUTDIR']
flank  = int(os.environ['FLANK'])   # window around BLAST hit for gene lookup

genes_by_chrom = defaultdict(list)
gene_info = {}
with open(f"{outdir}/tmp/genes.bed") as f:
    for line in f:
        fields = line.strip().split('\t')
        if len(fields) < 6: continue
        chrom, start, end, gene_id, gene_name, desc = fields[:6]
        # Strip gene: prefix just in case
        gene_id   = gene_id.replace('gene:', '').strip()
        gene_name = gene_name.replace('gene:', '').strip()
        gene_info[gene_id] = {'name': gene_name, 'desc': desc}
        genes_by_chrom[chrom].append((int(start), int(end), gene_id))

def find_genes(chrom, start, end):
    """Find genes overlapping the interval, with flank window for promoter hits."""
    hits = []
    for gs, ge, gid in genes_by_chrom.get(chrom, []):
        if gs <= end + flank and ge >= start - flank:
            # Calculate true overlap (without flank) for percentage
            true_overlap = max(0, min(end, ge) - max(start, gs))
            pct = true_overlap / (ge - gs) * 100 if (ge - gs) > 0 else 0
            # Flag if hit is in flank only (promoter/downstream)
            in_gene = gs <= end and ge >= start
            tag = '' if in_gene else '~'  # ~ = flanking/promoter
            hits.append((pct, f"{tag}{gid}({gene_info[gid]['name']},{pct:.0f}%)"))
    hits.sort(key=lambda x: x[0], reverse=True)
    gene_str = ';'.join(h[1] for h in hits[:5]) if hits else 'intergenic'
    is_genic = any(h[0] > 0 or '~' not in h[1] for h in hits) if hits else False
    return gene_str, len(hits) > 0

def annotate(blast_file, sv_type, out_file):
    header = ("sv_id\tsv_type\tn_hits\tmax_bitscore\ttotal_bitscore\t"
              "score_ratio\tcopy_number_in_insert\thit_chrom\thit_start\thit_end\t"
              "pident\tqcovs\toverlapping_genes\tis_genic\tis_same_chrom\n")

    if not os.path.exists(blast_file) or os.path.getsize(blast_file) == 0:
        print(f"  No BLAST results for {sv_type}")
        with open(out_file, 'w') as f: f.write(header)
        return

    query_hits = defaultdict(list)
    with open(blast_file) as f:
        for line in f:
            ff = line.strip().split('\t')
            if len(ff) < 13: continue
            q, s, pi, ln, ql, sl, qc, qs, qe, ss, se, ev, bs = ff
            query_hits[q].append(dict(
                sseqid=s, pident=float(pi), qcovs=float(qc),
                sstart=int(ss), send=int(se), bitscore=float(bs)
            ))

    results = []
    for qid, hits in query_hits.items():
        hits.sort(key=lambda x: x['bitscore'], reverse=True)
        max_s = hits[0]['bitscore']
        tot_s = sum(h['bitscore'] for h in hits)
        ratio = tot_s / max_s if max_s > 0 else 1
        cn    = round(ratio)
        best  = hits[0]
        # Handle reverse-strand hits (sstart > send)
        h_start = min(best['sstart'], best['send'])
        h_end   = max(best['sstart'], best['send'])
        genes, is_genic = find_genes(best['sseqid'], h_start, h_end)
        # same_chrom: SV chromosome matches hit chromosome
        sv_chrom = qid.split('_INS')[0].split('_DEL')[0]
        same_chrom = best['sseqid'] == sv_chrom
        results.append(dict(
            sv_id=qid, sv_type=sv_type, n_hits=len(hits),
            max_bs=max_s, tot_bs=tot_s, ratio=ratio, cn=cn,
            hit_chrom=best['sseqid'], hit_start=h_start, hit_end=h_end,
            pident=best['pident'], qcovs=best['qcovs'],
            genes=genes, is_genic=is_genic, same_chrom=same_chrom
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
    print(f"  {sv_type}: {len(results)} total | {g} genic | {t} putative tandem dups")

os.makedirs(f"{outdir}/results", exist_ok=True)
annotate(f"{outdir}/blast/insertions_blast_raw.txt", "INS",
         f"{outdir}/results/insertions_annotated.tsv")
annotate(f"{outdir}/blast/deletions_blast_raw.txt",  "DEL",
         f"{outdir}/results/deletions_annotated.tsv")
PYEOF

else
    echo "[STEP 11] Skipping (FROM_STEP=$FROM_STEP)"
fi

# =============================================================================
# STEP 12: Summary report
# =============================================================================
if should_run 12; then
    echo ""
    echo "[STEP 12] Generating summary report..."

    python3 << 'PYEOF'
import os, datetime

outdir = os.environ['OUTDIR']
flank  = int(os.environ['FLANK'])

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
        if gid in s: return drug
    return ''

lines = [
    "=" * 65,
    "PANGENOME SV PIPELINE - SUMMARY REPORT",
    f"Output:    {outdir}",
    f"Generated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}",
    f"Gene flank window: {flank}bp (~ = promoter/downstream hit)",
    "=" * 65,
]

for label, fname in [("INSERTIONS", "insertions_annotated.tsv"),
                      ("DELETIONS",  "deletions_annotated.tsv")]:
    fp = f"{outdir}/results/{fname}"
    if not os.path.exists(fp): continue
    with open(fp) as f:
        rows = [l.strip().split('\t') for l in f if not l.startswith('sv_id')]
    genic  = [r for r in rows if len(r) > 13 and r[13] == 'True']
    tandem = [r for r in rows if len(r) > 14 and r[13] == 'True' and r[14] == 'True']
    res    = [r for r in rows if len(r) > 12 and flag_res(r[12])]
    lines += [
        f"\n{label}",
        f"  Total with BLAST hits:        {len(rows)}",
        f"  Hitting annotated genes:      {len(genic)}",
        f"  Putative tandem duplications: {len(tandem)}",
        f"  Resistance gene hits:         {len(res)}",
    ]
    if res:
        lines.append(f"\n  *** RESISTANCE GENE {label} ***")
        for r in res:
            lines.append(f"    {r[0]:<45} [{flag_res(r[12])}]  cn={r[6]}  genes={r[12][:60]}")

fp = f"{outdir}/results/whole_gene_deletions.tsv"
if os.path.exists(fp):
    with open(fp) as f:
        rows = [l.strip().split('\t') for l in f if not l.startswith('sv_id')]
    res = [r for r in rows if len(r) > 5 and flag_res(r[5])]
    lines += [
        f"\nWHOLE GENE DELETIONS (>90% gene body covered)",
        f"  Total:                        {len(rows)}",
        f"  Resistance gene hits:         {len(res)}",
    ]
    if rows:
        lines.append(f"\n  {'SV_ID':<42} {'GENE_NAME':<22} {'OVERLAP%':>9}  FLAG")
        for r in rows[:25]:
            tag  = f"  *** {flag_res(r[5])}" if len(r) > 5 and flag_res(r[5]) else ''
            name = r[6] if len(r) > 6 else ''
            lines.append(f"  {r[0]:<42} {name:<22} {r[10]:>9}{tag}")

fp = f"{outdir}/results/partial_gene_deletions.tsv"
if os.path.exists(fp):
    with open(fp) as f:
        rows = [l.strip().split('\t') for l in f if not l.startswith('sv_id')]
    lines += [
        f"\nPARTIAL GENE DELETIONS (20-90% gene body covered)",
        f"  Total:                        {len(rows)}",
    ]

lines += [
    f"\nOUTPUT FILES",
    f"  {outdir}/results/insertions_annotated.tsv",
    f"  {outdir}/results/deletions_annotated.tsv",
    f"  {outdir}/results/whole_gene_deletions.tsv",
    f"  {outdir}/results/partial_gene_deletions.tsv",
    f"  {outdir}/logs/pipeline.log",
    "=" * 65,
]

report = '\n'.join(lines)
print(report)
with open(f"{outdir}/results/summary_report.txt", 'w') as f:
    f.write(report)
PYEOF

else
    echo "[STEP 12] Skipping (FROM_STEP=$FROM_STEP)"
fi

# Cleanup
rm -f "$OUTDIR/tmp/"*.fa 2>/dev/null || true

echo ""
echo "============================================="
echo "Pipeline complete: $(date)"
echo "Results: $OUTDIR/results/"
echo "============================================="