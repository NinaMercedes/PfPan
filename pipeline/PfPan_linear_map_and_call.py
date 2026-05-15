#!/usr/bin/env python3
"""
PfPan_linear_map_and_call.py
=============================
Linear reference mapping and variant calling pipeline.

Runs BWA-based mapping, BQSR, and GATK variant calling per sample via
fastq2matrix (fastq2vcf.py), followed by a multi-sample DELLY SV calling
workflow.

Requirements:
    conda activate <fastq2matrix_env>   — for fastq2vcf.py (Step 1)
    conda activate delly                — for DELLY steps (Steps 2-5)

    Note: Steps 1 and 2-5 use different conda environments. Run them
    separately or update FASTQ2MATRIX_ENV and DELLY_ENV in the config
    to activate environments automatically.

Steps:
  1. fastq2vcf.py all  — BWA mapping + BQSR + GATK HaplotypeCaller per sample
  2. delly call        — per-sample SV calling against the reference
  3. delly merge       — merge per-sample BCFs into a unified sites file
  4. delly call -v     — re-genotype all samples at merged sites
  5. bcftools merge    — merge all re-genotyped BCFs into a single multi-sample VCF

Usage:
    # Full pipeline
    python PfPan_linear_map_and_call.py --samples-file fastqs.txt

    # Skip individual steps
    python PfPan_linear_map_and_call.py --samples-file fastqs.txt --skip-fastq2vcf
    python PfPan_linear_map_and_call.py --samples-file fastqs.txt --skip-delly

Inputs:
    fastqs.txt          — one sample name per line (no extensions)
    <sample>_1.trimmed.fastq.gz  — trimmed reads (in ./reads/
    <sample>_2.trimmed.fastq.gz

Outputs written to ./linear/:
    <sample>.bqsr.bam              — mapped, recalibrated BAM (from fastq2vcf)
    <sample>_gatk.g.vcf.gz        — GATK gVCF (from fastq2vcf)
    <sample>_delly.bcf             — per-sample DELLY calls
    <sample>_delly_sites.bcf       — per-sample re-genotyped at merged sites
    sites.bcf                      — merged DELLY sites across all samples
    pan_delly.bcf                  — merged multi-sample DELLY BCF
    pan_delly.vcf.gz               — merged multi-sample DELLY VCF
    pan_delly.vcf.gz.tbi           — tabix index
"""

import argparse
import logging
import subprocess
import sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler("PfPan_linear_map_and_call.log"),
    ],
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Config — edit paths to match your environment
# ---------------------------------------------------------------------------

REF = "/mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/Pfalciparum.genome.fasta"

FASTQ_DIR     = str(Path(__file__).resolve().parent.parent / "reads")
FASTQ2VCF_PY  = str(Path(__file__).resolve().parent / "fastq2vcf.py")
OUT_DIR       = str(Path(__file__).resolve().parent.parent / "linear")

THREADS          = 20
DELLY_PARALLEL   = 10   # samples processed in parallel for DELLY steps
MAPPER           = "bwa"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run(cmd: str) -> None:
    log.info(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)


def run_sample(cmd: str) -> None:
    """Used for parallel per-sample execution — errors are logged not raised."""
    log.info(f"Running: {cmd}")
    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        log.error(f"Command failed (exit {result.returncode}): {cmd}")


def ensure_dir(path: str) -> None:
    Path(path).mkdir(parents=True, exist_ok=True)


def load_samples(samples_file: str) -> list:
    lines = Path(samples_file).read_text().splitlines()
    samples = [l.strip() for l in lines if l.strip()]
    log.info(f"Loaded {len(samples)} samples from {samples_file}")
    return samples


# ---------------------------------------------------------------------------
# Step 1: fastq2vcf — BWA mapping + BQSR + GATK per sample
# ---------------------------------------------------------------------------

def step_fastq2vcf(samples: list) -> None:
    """
    Run fastq2vcf.py for each sample sequentially (-P 1).
    Requires the fastq2matrix conda environment.
    Each sample produces a BQSR BAM and GATK gVCF in OUT_DIR.
    """
    log.info("=== Step 1: BWA mapping + BQSR + GATK (fastq2vcf.py) ===")
    log.info("NOTE: Requires fastq2matrix conda environment")

    for sample in samples:
        fq1 = f"{FASTQ_DIR}/{sample}_1.trimmed.fastq.gz"
        fq2 = f"{FASTQ_DIR}/{sample}_2.trimmed.fastq.gz"
        cmd = (
            f'python "{FASTQ2VCF_PY}" all '
            f"--read1 {fq1} "
            f"--read2 {fq2} "
            f"--prefix {sample} "
            f'--ref "{REF}" '
            f"--threads {THREADS} "
            f"--mapper {MAPPER}"
        )
        log.info(f"  -> {sample}")
        run_sample(cmd)


# ---------------------------------------------------------------------------
# Step 2: DELLY per-sample call
# ---------------------------------------------------------------------------

def step_delly_call(samples: list) -> None:
    """
    Run delly call per sample against the reference BAM in parallel.
    Requires the delly conda environment.
    """
    log.info("=== Step 2: DELLY per-sample SV calling ===")
    log.info("NOTE: Requires delly conda environment")

    cmds = [
        f"delly call "
        f'-g "{REF}" '
        f"./{sample}.bqsr.bam "
        f"-o {sample}_delly.bcf"
        for sample in samples
    ]

    with ThreadPoolExecutor(max_workers=DELLY_PARALLEL) as pool:
        list(pool.map(run_sample, cmds))


# ---------------------------------------------------------------------------
# Step 3: DELLY merge sites
# ---------------------------------------------------------------------------

def step_delly_merge() -> None:
    """Merge per-sample DELLY BCFs into a single unified sites file."""
    log.info("=== Step 3: DELLY merge sites ===")
    run("delly merge -o sites.bcf *_delly.bcf")


# ---------------------------------------------------------------------------
# Step 4: DELLY re-genotype at merged sites
# ---------------------------------------------------------------------------

def step_delly_regenotype(samples: list) -> None:
    """
    Re-genotype each sample at the merged site list (-v sites.bcf).
    This ensures all samples are called at the same loci for joint analysis.
    """
    log.info("=== Step 4: DELLY re-genotyping at merged sites ===")

    cmds = [
        f"delly call "
        f'-g "{REF}" '
        f"-v sites.bcf "
        f"./{sample}.bqsr.bam "
        f"-o {sample}_delly_sites.bcf"
        for sample in samples
    ]

    with ThreadPoolExecutor(max_workers=DELLY_PARALLEL) as pool:
        list(pool.map(run_sample, cmds))


# ---------------------------------------------------------------------------
# Step 5: Merge all re-genotyped BCFs into multi-sample VCF
# ---------------------------------------------------------------------------

def step_delly_merge_final() -> None:
    """
    Merge all per-sample re-genotyped BCFs into a single multi-sample BCF/VCF.
    Output: pan_delly.vcf.gz
    """
    log.info("=== Step 5: Merge all DELLY BCFs into multi-sample VCF ===")
    run("bcftools merge -m id -Ob -o pan_delly.bcf *_delly_sites.bcf")
    run("bcftools view -Oz -o pan_delly.vcf.gz pan_delly.bcf")
    run("tabix -p vcf pan_delly.vcf.gz")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Linear reference mapping and multi-sample DELLY SV pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--samples-file", required=True,
        help="Path to file listing sample names, one per line (e.g. fastqs.txt)"
    )
    parser.add_argument(
        "--skip-fastq2vcf", action="store_true",
        help="Skip Step 1 (fastq2vcf BWA mapping + GATK)"
    )
    parser.add_argument(
        "--skip-delly", action="store_true",
        help="Skip Steps 2-5 (all DELLY SV calling steps)"
    )
    args = parser.parse_args()

    samples = load_samples(args.samples_file)
    ensure_dir(OUT_DIR)

    import os
    os.chdir(OUT_DIR)
    log.info(f"Working directory: {Path(OUT_DIR).resolve()}")

    if not args.skip_fastq2vcf:
        step_fastq2vcf(samples)
    else:
        log.info("Skipping Step 1 (fastq2vcf)")

    if not args.skip_delly:
        step_delly_call(samples)
        step_delly_merge()
        step_delly_regenotype(samples)
        step_delly_merge_final()
    else:
        log.info("Skipping Steps 2-5 (DELLY)")

    log.info("=== PfPan_linear_map_and_call complete ===")


if __name__ == "__main__":
    main()
