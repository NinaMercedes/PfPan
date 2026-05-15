#!/usr/bin/env python3
"""
PfPan_PfPan_map_and_call.py
=========================
Per-sample pangenome read mapping and variant calling pipeline.

For each sample:
  1. Maps trimmed paired-end Illumina reads to the pangenome graph using vg giraffe (GAF output)
  2. Surjects graph alignments onto Pf3D7 linear reference paths to produce a BAM
  3. Sorts and indexes the BAM with samtools; runs flagstat
  4. Calls short variants (SNPs/indels) via GATK HaplotypeCaller in GVCF mode
  5. Calls SVs with DELLY from the sorted BAM
  6. Builds a pack coverage index from graph alignments (vg pack)
  7. Calls variants from the graph (vg call) — diploid
  8. Calls variants from the graph (vg call) — haploid (--ploidy 1)

Usage:
    python PfPan_PfPan_map_and_call.py --sample <sample_name>

    # Skip individual steps
    python PfPan_PfPan_map_and_call.py --sample <sample_name> --skip-delly
    python PfPan_PfPan_map_and_call.py --sample <sample_name> --skip-haploid
    python PfPan_PfPan_map_and_call.py --sample <sample_name> --skip-delly --skip-haploid

Expected inputs (in current directory):
    <sample>_1.trimmed.fastq.gz
    <sample>_2.trimmed.fastq.gz

Required reference files (edit config block below):
    PfPan.gbz        — pangenome graph
    PfPan.snarls     — snarl index for vg call
    Pf3D7.paths.txt        — list of Pf3D7 reference paths for surjection
    Pf3D7_graph.fa         — Pf3D7 linear reference FASTA for GATK + DELLY

Outputs written to ./output/:
    <sample>.gaf.gz                  — graph alignments (vg giraffe)
    <sample>.bam                     — unsorted surjected BAM
    <sample>_sort.bam                — sorted, indexed BAM
    <sample>_sort.bam.bai            — BAM index
    <sample>_stat.txt                — samtools flagstat summary
    <sample>_gatk.g.vcf.gz          — GATK HaplotypeCaller GVCF
    <sample>_delly_sites.bcf         — DELLY SV calls (BCF)
    <sample>_delly_sites.bcf.csi     — DELLY BCF index
    <sample>.pack                    — vg pack coverage index
    <sample>.SV.call.vcf.gz          — vg call diploid variant calls
    <sample>.SV.call_haploid_.vcf.gz — vg call haploid variant calls (--ploidy 1)
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
        logging.FileHandler("PfPan_map_and_call.log"),
    ],
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Config — edit paths to match your environment
# ---------------------------------------------------------------------------

GBZ        = "./PfPan_Pf3D7_pan/PfPan_test.gbz"
SNARLS     = "./PfPan_Pf3D7_pan/PfPan_test.snarls"
PATHS_FILE = "./Pf3D7.paths.txt"  # generate with: vg paths -x PfPan_Pf3D7_pan/PfPan_test.gbz -L | grep Pf3D7 > Pf3D7.paths.txt
REF_FA     = "/mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/Pfalciparum.genome.fasta"
OUT_DIR    = "./output"
THREADS    = 8
VG_THREADS = 16

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run(cmd: str) -> None:
    log.info(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)


def ensure_dir(path: str) -> None:
    Path(path).mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Pipeline steps
# ---------------------------------------------------------------------------

def graph_align(sample: str) -> None:
    """Map trimmed reads to the pangenome graph using vg giraffe."""
    log.info("--- Step 1: Graph alignment (vg giraffe) ---")
    fq1 = f"./reads/{sample}_1.trimmed.fastq.gz"
    fq2 = f"./reads/{sample}_2.trimmed.fastq.gz"
    out = f"{OUT_DIR}/{sample}.gaf.gz"
    run(
        f"vg giraffe -p -t {VG_THREADS} -Z {GBZ} "
        f"-f {fq1} -f {fq2} "
        f"-o gaf | bgzip > {out}"
    )


def surject_to_bam(sample: str) -> None:
    """
    Surject graph alignments onto Pf3D7 linear reference paths.
    Pipes through samtools reheader to strip pangenome-style contig prefixes
    (e.g. Pf3D7#0#Pf3D7_01_v3 -> Pf3D7_01_v3) so the BAM is compatible
    with the linear reference used by GATK and DELLY.
    """
    log.info("--- Step 2: Surjection to BAM (vg surject) ---")
    gaf = f"{OUT_DIR}/{sample}.gaf.gz"
    bam = f"{OUT_DIR}/{sample}.bam"
    run(
        f"vg surject -x {GBZ} -G {gaf} --interleaved "
        f"-F {PATHS_FILE} -b "
        f"-N {sample} "
        f"-R 'ID:1 LB:lib1 SM:{sample} PL:illumina PU:unit1' "
        f"| samtools reheader -c 'sed s/Pf3D7#0#//g' - "
        f"> {bam}"
    )


def sort_and_index_bam(sample: str) -> None:
    """Sort and index the surjected BAM, then run flagstat."""
    log.info("--- Step 3: Sort and index BAM (samtools) ---")
    bam      = f"{OUT_DIR}/{sample}.bam"
    sort_bam = f"{OUT_DIR}/{sample}_sort.bam"
    stat     = f"{OUT_DIR}/{sample}_stat.txt"
    run(f"samtools sort {bam} -O BAM -o {sort_bam} --threads {THREADS}")
    run(f"samtools index {sort_bam} --threads {THREADS}")
    run(f"samtools flagstat {sort_bam} > {stat}")


def gatk_haplotype_caller(sample: str) -> None:
    """Call short variants in GVCF mode using GATK HaplotypeCaller."""
    log.info("--- Step 4: Short variant calling (GATK HaplotypeCaller) ---")
    sort_bam = f"{OUT_DIR}/{sample}_sort.bam"
    gvcf     = f"{OUT_DIR}/{sample}_gatk.g.vcf.gz"
    run(
        f"gatk HaplotypeCaller "
        f"-I {sort_bam} "
        f"-R {REF_FA} "
        f"-O {gvcf} "
        f"-ERC GVCF"
    )


def delly_sv_call(sample: str) -> None:
    """
    Call structural variants from the sorted BAM using DELLY.
    Produces a BCF file which is filtered in Pipeline 2 (Step 1)
    to PASS SVs >= 50bp.
    """
    log.info("--- Step 5: SV calling (DELLY) ---")
    sort_bam = f"{OUT_DIR}/{sample}_sort.bam"
    bcf      = f"{OUT_DIR}/{sample}_delly_sites.bcf"
    run(
        f"delly call "
        f"-g {REF_FA} "
        f"-o {bcf} "
        f"{sort_bam}"
    )


def vg_pack(sample: str) -> None:
    """Build pack coverage index from graph alignments."""
    log.info("--- Step 6: Build pack index (vg pack) ---")
    gaf  = f"{OUT_DIR}/{sample}.gaf.gz"
    pack = f"{OUT_DIR}/{sample}.pack"
    run(f"vg pack -x {GBZ} -Q5 -a {gaf} -o {pack}")


def vg_call_diploid(sample: str) -> None:
    """Call variants from the graph in diploid mode (default ploidy=2)."""
    log.info("--- Step 7: Graph variant calling — diploid (vg call) ---")
    pack = f"{OUT_DIR}/{sample}.pack"
    vcf  = f"{OUT_DIR}/{sample}.SV.call.vcf.gz"
    run(
        f"vg call {GBZ} -r {SNARLS} -k {pack} "
        f"-t {VG_THREADS} "
        f"-s {sample} -S Pf3D7 -az | bgzip > {vcf}"
    )


def vg_call_haploid(sample: str) -> None:
    """
    Call variants from the graph in haploid mode (--ploidy 1).
    Useful for P. falciparum which is haploid during the blood stage.
    """
    log.info("--- Step 8: Graph variant calling — haploid (vg call --ploidy 1) ---")
    pack = f"{OUT_DIR}/{sample}.pack"
    vcf  = f"{OUT_DIR}/{sample}.SV.call_haploid_.vcf.gz"
    run(
        f"vg call {GBZ} -r {SNARLS} -k {pack} "
        f"-t {VG_THREADS} "
        f"--ploidy 1 "
        f"-s {sample} -S Pf3D7 -az | bgzip > {vcf}"
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run_pipeline(sample: str, skip_delly: bool, skip_haploid: bool) -> None:
    ensure_dir(OUT_DIR)
    log.info(f"=== Starting pipeline for sample: {sample} ===")

    graph_align(sample)
    surject_to_bam(sample)
    sort_and_index_bam(sample)
    gatk_haplotype_caller(sample)

    if not skip_delly:
        delly_sv_call(sample)
    else:
        log.info("Skipping Step 5 (DELLY SV calling)")

    vg_pack(sample)
    vg_call_diploid(sample)

    if not skip_haploid:
        vg_call_haploid(sample)
    else:
        log.info("Skipping Step 8 (haploid vg call)")

    log.info(f"=== Pipeline complete for sample: {sample} ===")


def main():
    parser = argparse.ArgumentParser(
        description="Pangenome read mapping and variant calling pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--sample",       required=True,      help="Sample name")
    parser.add_argument("--skip-delly",   action="store_true", help="Skip DELLY SV calling (Step 5)")
    parser.add_argument("--skip-haploid", action="store_true", help="Skip haploid vg call (Step 8)")
    args = parser.parse_args()
    run_pipeline(args.sample, args.skip_delly, args.skip_haploid)


if __name__ == "__main__":
    main()