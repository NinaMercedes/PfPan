#!/usr/bin/env python3
"""
Extract precise insertion/deletion sequences from pangenome VCF.
Aligns REF vs ALT using MAFFT, extracts gap blocks to find
inserted (in ALT) or deleted (in REF) sequences.

Usage:
    python3 extract_sv_sequences.py \
        --vcf sv_results/vcf/core_insertions.vcf.gz \
        --outdir sv_results \
        --svtype INS \
        --min-svlen 1000

    python3 extract_sv_sequences.py \
        --vcf sv_results/vcf/core_deletions.vcf.gz \
        --outdir sv_results \
        --svtype DEL \
        --min-svlen 1000
"""

import subprocess
import os
import re
import gzip
import sys
import argparse
import tempfile
from pathlib import Path

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--vcf",      required=True,  help="VCF/VCF.gz file")
    p.add_argument("--outdir",   required=True,  help="Output directory")
    p.add_argument("--svtype",   required=True,  choices=["INS","DEL"])
    p.add_argument("--min-svlen",type=int, default=1000)
    p.add_argument("--threads",  type=int, default=4)
    return p.parse_args()

def open_vcf(path):
    if str(path).endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')

def run_mafft(ref_seq, alt_seq, tmpdir):
    """
    Align ref and alt with MAFFT.
    Returns (aligned_ref, aligned_alt) strings.
    """
    in_fa = os.path.join(tmpdir, "input.fa")
    out_fa = os.path.join(tmpdir, "output.fa")

    with open(in_fa, 'w') as f:
        f.write(f">ref\n{ref_seq}\n>alt\n{alt_seq}\n")

    result = subprocess.run(
        ["mafft", "--auto", "--quiet", "--thread", "1", in_fa],
        capture_output=True, text=True
    )

    if result.returncode != 0 or not result.stdout.strip():
        return None, None

    # Parse fasta output
    seqs = {}
    current = None
    for line in result.stdout.strip().split('\n'):
        if line.startswith('>'):
            current = line[1:].strip()
            seqs[current] = []
        elif current:
            seqs[current].append(line.strip())

    aligned = {k: ''.join(v).upper() for k, v in seqs.items()}

    if 'ref' not in aligned or 'alt' not in aligned:
        return None, None

    return aligned['ref'], aligned['alt']

def extract_gap_blocks(seq1, seq2, min_len):
    """
    Find gap blocks in seq1 (gaps = '-').
    These represent sequence present in seq2 but absent in seq1.
    For INS: seq1=ref (has gaps where ALT has insertion)
    For DEL: seq1=alt (has gaps where REF has deletion)
    Returns list of (aln_start, aln_end, sequence_from_seq2, length)
    """
    blocks = []
    i = 0
    while i < len(seq1):
        if seq1[i] == '-':
            # Start of gap block
            j = i
            while j < len(seq1) and seq1[j] == '-':
                j += 1
            gap_len = j - i
            extracted_seq = seq2[i:j].replace('-', '')
            if len(extracted_seq) >= min_len:
                blocks.append((i, j, extracted_seq, len(extracted_seq)))
            i = j
        else:
            i += 1
    return blocks

def process_vcf(vcf_path, outdir, svtype, min_svlen, threads):

    sequences_dir = os.path.join(outdir, "sequences")
    logs_dir      = os.path.join(outdir, "logs")
    os.makedirs(sequences_dir, exist_ok=True)
    os.makedirs(logs_dir, exist_ok=True)

    # Output files
    label = svtype.lower()
    out_fa      = os.path.join(sequences_dir, f"{label}_sequences.fa")
    out_summary = os.path.join(sequences_dir, f"{label}_summary.tsv")
    out_failed  = os.path.join(logs_dir,      f"{label}_extraction_failed.txt")

    # Read VCF records
    records = []
    with open_vcf(vcf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue

            chrom  = fields[0]
            pos    = int(fields[1])
            ref    = fields[3]
            alt    = fields[4]
            info   = fields[7]

            # Parse SVLEN and SVTYPE from INFO
            svlen   = None
            sv_type = None
            for item in info.split(';'):
                if item.startswith('SVLEN='):
                    try:
                        svlen = int(item.split('=')[1])
                    except ValueError:
                        pass
                if item.startswith('SVTYPE='):
                    sv_type = item.split('=')[1]

            if svlen is None:
                # Fall back to length difference
                svlen = len(alt) - len(ref)

            if abs(svlen) < min_svlen:
                continue

            # Skip if REF or ALT is symbolic (<INS>, <DEL>, etc.)
            if alt.startswith('<') or ref.startswith('<'):
                continue

            # Skip if sequences are too short to be meaningful
            if len(ref) < 10 or len(alt) < 10:
                continue

            sv_id = f"{chrom}_{pos}_{svtype}{abs(svlen)}"
            records.append({
                'id':    sv_id,
                'chrom': chrom,
                'pos':   pos,
                'ref':   ref,
                'alt':   alt,
                'svlen': svlen,
            })

    print(f"  Found {len(records)} {svtype} records >= {min_svlen}bp")

    if not records:
        print("  No records to process. Check SVTYPE/SVLEN INFO fields in VCF.")
        # Debug: show first few INFO fields
        print("\n  DEBUG - First 3 INFO fields from VCF:")
        with open_vcf(vcf_path) as f:
            n = 0
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) >= 8:
                    print(f"    {fields[7][:200]}")
                    n += 1
                    if n >= 3:
                        break
        return

    extracted = []
    failed    = []

    with tempfile.TemporaryDirectory() as tmpdir:
        for i, rec in enumerate(records):
            sv_id = rec['id']

            # For INS: gaps in REF alignment = sequence in ALT not in REF
            # For DEL: gaps in ALT alignment = sequence in REF not in ALT
            aln_ref, aln_alt = run_mafft(rec['ref'], rec['alt'], tmpdir)

            if aln_ref is None:
                failed.append(sv_id)
                continue

            if svtype == "INS":
                # Look for gaps in REF (= insertions in ALT)
                blocks = extract_gap_blocks(aln_ref, aln_alt, min_svlen)
            else:
                # Look for gaps in ALT (= deletions relative to REF)
                blocks = extract_gap_blocks(aln_alt, aln_ref, min_svlen)

            if not blocks:
                # No gap block >= min_svlen found
                # Could be alignment didn't resolve it cleanly
                # Try lower threshold to debug
                if svtype == "INS":
                    all_blocks = extract_gap_blocks(aln_ref, aln_alt, 1)
                else:
                    all_blocks = extract_gap_blocks(aln_alt, aln_ref, 1)

                if all_blocks:
                    largest = max(all_blocks, key=lambda x: x[3])
                    failed.append(f"{sv_id} (largest_gap={largest[3]}bp)")
                else:
                    failed.append(f"{sv_id} (no_gaps_found)")
                continue

            # Take the largest gap block as the SV sequence
            blocks.sort(key=lambda x: x[3], reverse=True)
            aln_start, aln_end, sv_seq, sv_len = blocks[0]

            extracted.append({
                'id':       sv_id,
                'chrom':    rec['chrom'],
                'pos':      rec['pos'],
                'svlen':    rec['svlen'],
                'aln_pos':  aln_start,
                'ext_len':  sv_len,
                'seq':      sv_seq,
            })

            if (i + 1) % 20 == 0:
                print(f"  Processed {i+1}/{len(records)}...", flush=True)

    # Write outputs
    with open(out_fa, 'w') as f:
        for r in extracted:
            f.write(f">{r['id']}\n{r['seq']}\n")

    with open(out_summary, 'w') as f:
        f.write("sv_id\tchrom\tpos\tsvlen\taln_pos\textracted_len\n")
        for r in extracted:
            f.write(f"{r['id']}\t{r['chrom']}\t{r['pos']}\t"
                    f"{r['svlen']}\t{r['aln_pos']}\t{r['ext_len']}\n")

    if failed:
        with open(out_failed, 'w') as f:
            f.write('\n'.join(failed) + '\n')

    print(f"\n  Successfully extracted: {len(extracted)}")
    print(f"  Failed:                 {len(failed)}")
    print(f"  Sequences: {out_fa}")
    print(f"  Summary:   {out_summary}")
    if failed:
        print(f"  Failed:    {out_failed}")

if __name__ == "__main__":
    args = parse_args()
    process_vcf(
        vcf_path  = args.vcf,
        outdir    = args.outdir,
        svtype    = args.svtype,
        min_svlen = args.min_svlen,
        threads   = args.threads,
    )
