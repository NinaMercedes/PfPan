#!/usr/bin/env python3
"""
PfPan_SVIM_ASM.py
======================
Assembly-based SV and SNP calling pipeline using SVIM-asm and paftools.

Generates the truth set used for benchmarking in Pipeline 2:
  - SVIM-asm SV calls (used as the SV truth baseline by Truvari)
  - paftools SNP/indel calls (used as the short variant baseline by vcfeval)

For each assembly:
  1. Aligns the assembly to Pf3D7 using minimap2 (asm5 preset) → sorted BAM
  2. Calls SVs from the BAM using svim-asm (haploid mode)
  3. Aligns the assembly to Pf3D7 using minimap2 (PAF output) → sorted PAF
  4. Calls SNPs/indels from the PAF using paftools, normalises with bcftools

Usage:
    python PfPan_SVIM_ASM.py

    # Process a single strain only
    python PfPan_SVIM_ASM.py --strain PfDd2

    # Skip copying/unzipping FASTA files (already present)
    python PfPan_SVIM_ASM.py --skip-copy

Inputs:
    Assembly FASTA files (gzipped) — paths defined in ASSEMBLY_SRCS config
    Pf3D7 reference FASTA (gzipped) — path defined in REF_SRC config

Outputs written to ./svim_asm_results/<strain>/:
    <strain>.vs.Pf3D7.sorted.bam        — assembly-to-reference BAM
    <strain>.vs.Pf3D7.sorted.bam.bai    — BAM index
    svim/                               — SVIM-asm output directory
        variants.vcf                    — raw SV calls
        ...
    <strain>.paftools.snps.vcf.gz       — paftools SNP/indel calls (normalised)
    <strain>.paftools.snps.vcf.gz.tbi   — tabix index
"""

import argparse
import logging
import subprocess
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler("PfPan_SVIM_ASM.log"),
    ],
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Config — edit paths to match your environment
# ---------------------------------------------------------------------------

THREADS  = 16
SORT_MEM = "4G"
OUT_DIR  = "svim_asm_results"

REF_SRC = "/mnt/storage13/nbillows/pangenome/genomes/Pf3D7.April2018.fasta.gz"

ASSEMBLY_SRCS = [
    "/mnt/storage13/nbillows/pangenome/genomes/Pf7G8.April2018.fasta.gz",
    "/mnt/storage13/nbillows/pangenome/genomes/PfCD01.April2018.fasta.gz",
    "/mnt/storage13/nbillows/pangenome/genomes/PfDd2.April2018.fasta.gz",
    "/mnt/storage13/nbillows/pangenome/genomes/PfGA01.April2018.fasta.gz",
    "/mnt/storage13/nbillows/pangenome/genomes/PfGB4.April2018.fasta.gz",
    "/mnt/storage13/nbillows/pangenome/genomes/PfGN01.April2018.fasta.gz",
    "/mnt/storage13/nbillows/pangenome/genomes/PfHB3.April2018.fasta.gz",
    "/mnt/storage13/nbillows/pangenome/genomes/PfIT.April2018.fasta.gz",
    "/mnt/storage13/nbillows/pangenome/genomes/PfKE01.April2018.fasta.gz",
    "/mnt/storage13/nbillows/pangenome/genomes/PfKH01.April2018.fasta.gz",
    "/mnt/storage13/nbillows/pangenome/genomes/PfKH02.April2018.fasta.gz",
    "/mnt/storage13/nbillows/pangenome/genomes/PfML01.April2018.fasta.gz",
    "/mnt/storage13/nbillows/pangenome/genomes/PfSD01.April2018.fasta.gz",
    "/mnt/storage13/nbillows/pangenome/genomes/PfSN01.April2018.fasta.gz",
    "/mnt/storage13/nbillows/pangenome/genomes/PfTG01.April2018.fasta.gz",
]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run(cmd: str) -> None:
    log.info(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)


def ensure_dir(path: str) -> None:
    Path(path).mkdir(parents=True, exist_ok=True)


def copy_and_unzip(src: str) -> str:
    """Copy a gzipped FASTA to the working directory and decompress it."""
    gz   = Path(src).name
    fa   = gz.replace(".gz", "")
    if not Path(fa).exists():
        run(f"cp {src} . && gunzip -k {gz}")
    else:
        log.info(f"  Already present: {fa}")
    return fa


# ---------------------------------------------------------------------------
# Per-sample steps
# ---------------------------------------------------------------------------

def align_bam(sample: str, asm_fa: str, ref_fa: str) -> str:
    """Align assembly to reference with minimap2 (asm5), produce sorted BAM."""
    log.info(f"  --- BAM alignment (minimap2 asm5) ---")
    sample_dir = f"{OUT_DIR}/{sample}"
    ensure_dir(sample_dir)

    sam     = f"{sample_dir}/{sample}.vs.Pf3D7.sam"
    bam     = f"{sample_dir}/{sample}.vs.Pf3D7.sorted.bam"

    run(
        f"minimap2 -a -x asm5 --cs -r2k -t {THREADS} "
        f"{ref_fa} {asm_fa} > {sam}"
    )
    run(
        f"samtools sort -m {SORT_MEM} -@ {THREADS} "
        f"-o {bam} {sam}"
    )
    run(f"samtools index {bam}")
    Path(sam).unlink()
    return bam


def svim_asm_call(sample: str, bam: str, ref_fa: str) -> None:
    """Call SVs from assembly BAM using svim-asm in haploid mode."""
    log.info(f"  --- SV calling (svim-asm haploid) ---")
    svim_dir = f"{OUT_DIR}/{sample}/svim"
    run(f"svim-asm haploid {svim_dir} {bam} {ref_fa}")


def paftools_snp_call(sample: str, asm_fa: str, ref_fa: str) -> None:
    """
    Call SNPs/indels using paftools:
      - Align assembly to reference with minimap2 (PAF + cs string)
      - Sort PAF by reference position
      - Call variants with paftools.js
      - Normalise with bcftools and compress
    """
    log.info(f"  --- SNP/indel calling (paftools) ---")
    sample_dir = f"{OUT_DIR}/{sample}"
    paf        = f"{sample_dir}/{sample}.vs.Pf3D7.paf"
    paf_sorted = f"{sample_dir}/{sample}.vs.Pf3D7.sorted.paf"
    snp_vcf    = f"{sample_dir}/{sample}.paftools.snps.vcf.gz"

    run(
        f"minimap2 -c --cs -x asm5 -t {THREADS} "
        f"{ref_fa} {asm_fa} > {paf}"
    )
    run(f"sort -k6,6 -k8,8n {paf} > {paf_sorted}")
    Path(paf).unlink()

    run(
        f"paftools.js call -f {ref_fa} {paf_sorted} "
        f"| bcftools norm -m -any "
        f"| bcftools view -O z -o {snp_vcf}"
    )
    run(f"tabix -f -p vcf {snp_vcf}")
    Path(paf_sorted).unlink()


def process_sample(sample: str, asm_fa: str, ref_fa: str) -> None:
    log.info(f"=== Processing: {sample} ===")
    bam = align_bam(sample, asm_fa, ref_fa)
    svim_asm_call(sample, bam, ref_fa)
    paftools_snp_call(sample, asm_fa, ref_fa)
    log.info(f"=== Finished: {sample} ===")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Assembly-based SV and SNP truth set generation (SVIM-asm + paftools)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--strain", default=None,
        help="Process only this strain (e.g. PfDd2). Default: all strains."
    )
    parser.add_argument(
        "--skip-copy", action="store_true",
        help="Skip copying/unzipping FASTA files (use already-present local copies)"
    )
    args = parser.parse_args()

    ensure_dir(OUT_DIR)

    # Set up reference
    if args.skip_copy:
        ref_fa = Path(REF_SRC).name.replace(".gz", "")
        if not Path(ref_fa).exists():
            raise FileNotFoundError(
                f"Reference FASTA not found locally: {ref_fa}\n"
                "Run without --skip-copy to copy and decompress it."
            )
    else:
        log.info("=== Copying and decompressing FASTA files ===")
        ref_fa = copy_and_unzip(REF_SRC)

    # Filter to single strain if requested
    sources = ASSEMBLY_SRCS
    if args.strain:
        sources = [s for s in ASSEMBLY_SRCS if args.strain in s]
        if not sources:
            raise ValueError(
                f"Strain '{args.strain}' not found in ASSEMBLY_SRCS. "
                f"Available: {[Path(s).name.split('.')[0] for s in ASSEMBLY_SRCS]}"
            )

    # Process each assembly
    for src in sources:
        gz     = Path(src).name
        sample = gz.split(".")[0]   # e.g. PfDd2 from PfDd2.April2018.fasta.gz

        if args.skip_copy:
            asm_fa = gz.replace(".gz", "")
            if not Path(asm_fa).exists():
                log.warning(f"Assembly FASTA not found locally: {asm_fa}, skipping")
                continue
        else:
            asm_fa = copy_and_unzip(src)

        process_sample(sample, asm_fa, ref_fa)

    log.info("=== All samples complete ===")


if __name__ == "__main__":
    main()
