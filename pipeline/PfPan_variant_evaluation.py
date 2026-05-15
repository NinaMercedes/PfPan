#!/usr/bin/env python3
"""
PfPan_PfPan_variant_evaluation.py
================================
Post-calling variant filtering, benchmarking, and evaluation pipeline.

Runs the following steps in order:

  Step 1 — filter_delly_gatk:
      Filter DELLY BCFs to PASS SVs >=50bp.
      Process GATK gVCFs through the full hard-filter pipeline,
      producing normalised SNP+indel VCFs <50bp with genotype recoding.

  Step 2 — filter_pangenome_vcfs:
      Filter diploid vg call VCFs into short variants (<50bp) and SVs (>=50bp).

  Step 3 — filter_pangenome_vcfs_haploid:
      Same as Step 2 but for haploid vg call VCFs.

  Step 4 — recode_haploid_genotypes:
      Recode haploid GTs (0 → 0/0, 1 → 1/1) in vcfeval/truvari output directories
      so downstream tools expecting diploid FORMAT fields work correctly.

  Step 5 — prebench:
      Prepare comparison VCFs for benchmarking:
        - Rename and filter SVIM-asm truth VCFs to SVs >=50bp
        - Split pangenome SV VCFs into SNP/indel and SV subsets (hom-alt only)
        - Copy paftools SNP/indel baselines

  Step 6 — decompose:
      Decompose MNPs and complex indels in all VCFs using rtg vcfdecompose.

  Step 7 — vcfeval:
      Run RTG vcfeval to benchmark short variants (SNPs/indels) against
      paftools baselines, across all BED regions and callsets.
      Optionally applies DP>5 filter to specified callsets.

  Step 8 — truvari:
      Run Truvari bench to benchmark SVs (>=50bp) against SVIM-asm truth,
      across all BED regions and callsets. Filters to ALT-only calls
      (keeping mixed genotypes) before benchmarking.

  Step 9 — plot:
      Calls plot_evaluation.R to generate F-measure boxplot figures and
      summary CSVs for small variants and SVs across genomic regions.

Usage:
    # Full pipeline
    python PfPan_PfPan_variant_evaluation.py

    # Individual steps
    python PfPan_PfPan_variant_evaluation.py --steps filter_delly_gatk filter_pangenome_vcfs
    python PfPan_PfPan_variant_evaluation.py --steps prebench decompose vcfeval truvari plot
    python PfPan_PfPan_variant_evaluation.py --steps vcfeval truvari plot

    # Skip steps
    python PfPan_PfPan_variant_evaluation.py --skip recode_haploid_genotypes plot
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
        logging.FileHandler("PfPan_variant_evaluation.log"),
    ],
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Config — edit paths to match your environment
# ---------------------------------------------------------------------------

REF          = "/mnt/storage13/nbillows/Pf_09_24/Pfalciparum_09_24_v2/reference/Pf3D7_v3/Pfalciparum.genome.fasta"
REF_SDF      = "/mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/Pfalciparum.genome.fasta.sdf"
REF_FASTA    = "/mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/Pfalciparum.genome.fasta"

ILLUMINA_BASE = "/mnt/storage13/nbillows/pangenome/analysis/pan_GT/og_reads/merged_illumina/nina_test/output"
SVIM_BASE     = "/mnt/storage13/nbillows/pangenome/analysis/pan_GT/svim_asm_results"

VCFEVAL_INDIR  = "/mnt/storage13/nbillows/pangenome/analysis/pan_GT/illumina_comparison/vcfeval_decomposed"
TRUVARI_INDIR  = "/mnt/storage13/nbillows/pangenome/analysis/pan_GT/illumina_comparison/truvari"
BED_DIR        = "/mnt/storage13/nbillows/pangenome/analysis/bed_files"

PREBENCH_OUTDIR   = "comparison_vcfs_illumina"
DECOMPOSE_INDIR   = "vcfeval"
DECOMPOSE_OUTDIR  = "vcfeval_decomposed"
VCFEVAL_OUTDIR    = "rtg_vcfeval_multi_compare_dp5"
TRUVARI_OUTDIR    = "truvari_sv_eval_mixed_calls"
PLOT_SCRIPT       = str(Path(__file__).resolve().parent / "plot_evaluation.R")

STRAINS = [
    "Pf7G8", "PfCD01", "PfDd2", "PfGA01", "PfGB4",
    "PfGN01", "PfHB3", "PfIT", "PfKE01", "PfKH01",
    "PfKH02", "PfML01", "PfSN01", "PfTG01",
]

# vcfeval callset suffixes
VCFEVAL_CALLSETS = {
    "pan_snps":          "pan.snps_indels.decomposed.vcf.gz",
    "pan_shortvars":     "pan.shortvars.lt50bp.PASS.decomposed.vcf.gz",
    "pan_shortvars_hap": "pan.shortvars.lt50bp.PASS.hap.decomposed.vcf.gz",
    "gatk_shortvars":    "GATK.shortvars.lt50bp.PASS.GT.decomposed.vcf.gz",
}
# These callsets get an additional DP>5 filter before vcfeval
DP_FILTER_CALLSETS = {"pan_shortvars", "pan_shortvars_hap", "gatk_shortvars"}

# Truvari callset suffixes (>=50bp SVs)
TRUVARI_CALLSETS = {
    "delly":       "delly.SV.ge50bp.PASS.vcf.gz",
    "pan_SV":      "pan.SV.ge50bp.PASS.vcf.gz",
    "pan_SV_hap":  "pan.SV.ge50bp.PASS.hap.vcf.gz",
    "pan_svs":     "pan.svs.vcf.gz",
}

# Directories containing haploid VCFs to recode
RECODE_DIRS = [
    "/mnt/storage13/nbillows/pangenome/analysis/pan_GT/illumina_comparison/vcfeval_decomposed",
    "/mnt/storage13/nbillows/pangenome/analysis/pan_GT/illumina_comparison/truvari",
]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run(cmd: str) -> None:
    log.info(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)


def ensure_dir(path: str) -> None:
    Path(path).mkdir(parents=True, exist_ok=True)


def remove(*paths: str) -> None:
    for p in paths:
        for f in [p, p + ".tbi"]:
            try:
                Path(f).unlink()
            except FileNotFoundError:
                pass


def tabix(vcf: str) -> None:
    run(f"tabix -f -p vcf {vcf}")


# ---------------------------------------------------------------------------
# Step 1: Filter DELLY + GATK variants
# ---------------------------------------------------------------------------

def filter_delly(bcf: str) -> None:
    sample = Path(bcf).name.replace("_delly_sites.bcf", "")
    log.info(f"  DELLY -> {sample}")
    out = f"{sample}.delly.SV.ge50bp.PASS.vcf.gz"
    run(
        f"bcftools view {bcf} -f PASS "
        f"-i 'INFO/SVTYPE!=\"\" && abs(INFO/SVLEN)>=50' "
        f"-Oz -o {out}"
    )
    tabix(out)


def filter_gatk(gvcf: str) -> None:
    sample = Path(gvcf).name.replace(".g.vcf.gz", "")
    log.info(f"  GATK -> {sample}")

    raw    = f"{sample}.gatk.raw.vcf.gz"
    snps   = f"{sample}.gatk.snps.vcf.gz"
    indels = f"{sample}.gatk.indels.vcf.gz"
    snps_f = f"{sample}.gatk.snps.filtered.vcf.gz"
    ind_f  = f"{sample}.gatk.indels.filtered.vcf.gz"
    comb   = f"{sample}.gatk.filtered.vcf.gz"
    passed = f"{sample}.gatk.PASS.vcf.gz"
    normed = f"{sample}.gatk.norm.vcf.gz"
    short  = f"{sample}.GATK.shortvars.lt50bp.PASS.vcf.gz"
    final  = f"{sample}.GATK.shortvars.lt50bp.PASS.GT.vcf.gz"

    # Genotype gVCF
    run(f"gatk GenotypeGVCFs -R {REF} -V {gvcf} -O {raw}")
    tabix(raw)

    # Separate SNPs and indels
    run(f"gatk SelectVariants -R {REF} -V {raw} --select-type-to-include SNP -O {snps}")
    tabix(snps)
    run(f"gatk SelectVariants -R {REF} -V {raw} --select-type-to-include INDEL -O {indels}")
    tabix(indels)

    # Hard-filter SNPs (GATK recommended thresholds)
    run(
        f"gatk VariantFiltration -R {REF} -V {snps} "
        f"--filter-expression 'QD < 2.0'            --filter-name 'QD2' "
        f"--filter-expression 'QUAL < 30.0'          --filter-name 'QUAL30' "
        f"--filter-expression 'SOR > 4.0'            --filter-name 'SOR4' "
        f"--filter-expression 'FS > 60.0'            --filter-name 'FS60' "
        f"--filter-expression 'MQ < 40.0'            --filter-name 'MQ40' "
        f"--filter-expression 'MQRankSum < -15.0'    --filter-name 'MQRankSum-15' "
        f"--filter-expression 'ReadPosRankSum < -5.0' --filter-name 'ReadPosRankSum-5' "
        f"-O {snps_f}"
    )
    tabix(snps_f)

    # Hard-filter indels (GATK recommended thresholds)
    run(
        f"gatk VariantFiltration -R {REF} -V {indels} "
        f"--filter-expression 'QD < 2.0'                --filter-name 'QD2' "
        f"--filter-expression 'FS > 200.0'               --filter-name 'FS200' "
        f"--filter-expression 'ReadPosRankSum < -20.0'  --filter-name 'ReadPosRankSum-20' "
        f"-O {ind_f}"
    )
    tabix(ind_f)

    # Merge filtered SNPs + indels, keep PASS only, normalise
    run(f"bcftools concat -a {snps_f} {ind_f} -Oz -o {comb}")
    tabix(comb)
    run(f"gatk SelectVariants -R {REF} -V {comb} --exclude-filtered -O {passed}")
    tabix(passed)
    run(f"bcftools norm -f {REF} -m -both {passed} -Oz -o {normed}")
    tabix(normed)

    # Keep SNPs + indels <50bp
    run(
        f"bcftools view "
        f"-i 'TYPE=\"snp\" || (TYPE=\"indel\" && abs(strlen(ALT)-strlen(REF))<50)' "
        f"{normed} -Oz -o {short}"
    )
    tabix(short)

    # Recode genotypes by allele fraction (threshold 0.7)
    run(
        f"bcftools view {short} | setGT.py --fraction 0.7 "
        f"| bcftools view -O z -c 1 -o {final}"
    )

    # Cleanup intermediates
    remove(raw, snps, indels, snps_f, ind_f, comb, passed, normed)
    log.info(f"  Output: {final}")


def step_filter_delly_gatk() -> None:
    log.info("=== Step 1: Filter DELLY BCFs and GATK gVCFs ===")
    bcfs  = sorted(Path(".").glob("*_delly_sites.bcf"))
    gvcfs = sorted(Path(".").glob("*.g.vcf.gz"))
    if not bcfs:
        log.warning("No DELLY BCF files found (*_delly_sites.bcf)")
    for bcf in bcfs:
        filter_delly(str(bcf))
    if not gvcfs:
        log.warning("No GATK gVCF files found (*.g.vcf.gz)")
    for gvcf in gvcfs:
        filter_gatk(str(gvcf))


# ---------------------------------------------------------------------------
# Step 2: Filter pangenome diploid VCFs
# ---------------------------------------------------------------------------

def filter_pan_vcf(vcf: str, suffix_in: str, suffix_out: str) -> None:
    sample = Path(vcf).name.replace(suffix_in, "")
    log.info(f"  Pangenome (diploid) -> {sample}")

    norm_unsorted = f"{sample}.pan.norm.unsorted{suffix_out}.vcf.gz"
    norm          = f"{sample}.pan.norm{suffix_out}.vcf.gz"
    snps          = f"{sample}.pan.snps.PASS{suffix_out}.vcf.gz"
    ind           = f"{sample}.pan.indels.lt50bp.PASS{suffix_out}.vcf.gz"
    short         = f"{sample}.pan.shortvars.lt50bp.PASS{suffix_out}.vcf.gz"
    svs           = f"{sample}.pan.SV.ge50bp.PASS{suffix_out}.vcf.gz"

    run(f"tabix -f {vcf}")
    run(f"bcftools norm -f {REF} -m -both {vcf} -Oz -o {norm_unsorted}")
    run(f"bcftools sort {norm_unsorted} -Oz -o {norm}")
    tabix(norm)
    Path(norm_unsorted).unlink(missing_ok=True)

    run(f"bcftools view -f PASS -v snps {norm} -Oz -o {snps}")
    tabix(snps)

    run(
        f"bcftools view -f PASS -v indels "
        f"-i 'abs(strlen(ALT)-strlen(REF))<50' "
        f"{norm} -Oz -o {ind}"
    )
    tabix(ind)

    run(f"bcftools concat -a {snps} {ind} -Oz -o {short}")
    tabix(short)

    run(
        f"bcftools view -f PASS "
        f"-i 'abs(strlen(ALT)-strlen(REF))>=50' "
        f"{norm} -Oz -o {svs}"
    )
    tabix(svs)

    remove(norm, snps, ind)


def step_filter_pangenome_vcfs() -> None:
    log.info("=== Step 2: Filter pangenome diploid VCFs ===")
    vcfs = sorted(Path(".").glob("*.SV.call.vcf.gz"))
    if not vcfs:
        log.warning("No diploid pangenome VCFs found (*.SV.call.vcf.gz)")
    for vcf in vcfs:
        filter_pan_vcf(str(vcf), ".SV.call.vcf.gz", "")


# ---------------------------------------------------------------------------
# Step 3: Filter pangenome haploid VCFs
# ---------------------------------------------------------------------------

def step_filter_pangenome_vcfs_haploid() -> None:
    log.info("=== Step 3: Filter pangenome haploid VCFs ===")
    vcfs = sorted(Path(".").glob("*.SV.call_haploid_.vcf.gz"))
    if not vcfs:
        log.warning("No haploid pangenome VCFs found (*.SV.call_haploid_.vcf.gz)")
    for vcf in vcfs:
        filter_pan_vcf(str(vcf), ".SV.call_haploid_.vcf.gz", ".hap")


# ---------------------------------------------------------------------------
# Step 4: Recode haploid genotypes (0 → 0/0, 1 → 1/1)
# ---------------------------------------------------------------------------

def step_recode_haploid_genotypes() -> None:
    """
    Recode haploid GTs to diploid style in vcfeval and truvari directories.
    Required because vcfeval and truvari expect diploid FORMAT/GT fields.
    """
    log.info("=== Step 4: Recode haploid genotypes ===")
    for d in RECODE_DIRS:
        vcfs = sorted(Path(d).glob("*.hap*.vcf.gz"))
        if not vcfs:
            log.warning(f"No haploid VCFs found in {d}")
            continue
        for vcf in vcfs:
            vcf = str(vcf)
            tmp = vcf.replace(".vcf.gz", ".tmp.vcf.gz")
            log.info(f"  Recoding: {vcf}")
            run(
                f"bcftools +setGT {vcf} -- -t q -n 'c:0/0' -i 'GT=\"0\"' "
                f"| bcftools +setGT -- -t q -n 'c:1/1' -i 'GT=\"1\"' "
                f"| bgzip -c > {tmp}"
            )
            Path(tmp).replace(vcf)
            tabix(vcf)


# ---------------------------------------------------------------------------
# Step 5: Prebench — prepare comparison VCFs
# ---------------------------------------------------------------------------

def step_prebench() -> None:
    """
    Prepare VCFs for benchmarking:
      - Rename SVIM-asm VCFs and fix sample names
      - Filter SVIM-asm to SVs >=50bp
      - Split pangenome SV VCFs into SNP/indel and SV subsets (hom-alt only)
      - Copy paftools SNP/indel baselines
    """
    log.info("=== Step 5: Prebench — prepare comparison VCFs ===")
    ensure_dir(PREBENCH_OUTDIR)

    # Rename SVIM-asm VCFs
    log.info("  Renaming SVIM-asm VCFs...")
    for strain in STRAINS:
        in_vcf  = f"{SVIM_BASE}/{strain}.April2018/svim/variants.vcf.gz"
        out_vcf = f"{SVIM_BASE}/{strain}.April2018/svim/{strain}.svim.vcf.gz"
        if not Path(in_vcf).exists():
            log.warning(f"  WARNING: {in_vcf} not found, skipping")
            continue
        if not Path(out_vcf).exists():
            rename_txt = "rename.txt"
            Path(rename_txt).write_text(f"Sample\t{strain}\n")
            run(f"bcftools reheader -s {rename_txt} {in_vcf} -o {out_vcf}")
            tabix(out_vcf)
            Path(rename_txt).unlink(missing_ok=True)

    # Filter SVIM-asm to SVs >=50bp
    log.info("  Filtering SVIM-asm to SVs >=50bp...")
    for strain in STRAINS:
        in_vcf  = f"{SVIM_BASE}/{strain}.April2018/svim/{strain}.svim.vcf.gz"
        out_vcf = f"{PREBENCH_OUTDIR}/{strain}.svim.svs.vcf.gz"
        if not Path(in_vcf).exists():
            log.warning(f"  WARNING: {in_vcf} not found, skipping")
            continue
        run(
            f"bcftools view "
            f"-i 'abs(strlen(ALT)-strlen(REF))>=50' "
            f"-Oz -o {out_vcf} {in_vcf}"
        )
        tabix(out_vcf)

    # Process pangenome SV VCFs — split into SNP/indel and SV subsets (hom-alt only)
    log.info("  Processing pangenome SV VCFs...")
    for strain in STRAINS:
        in_vcf = f"{ILLUMINA_BASE}/{strain}.SV.call.vcf.gz"
        if not Path(in_vcf).exists():
            log.warning(f"  WARNING: {in_vcf} not found, skipping")
            continue
        tmp = f"{PREBENCH_OUTDIR}/{strain}.illumina.norm.11.vcf.gz"
        run(
            f"bcftools norm -m -any {in_vcf} "
            f"| bcftools view -i 'GT=\"1/1\"' -Oz -o {tmp}"
        )
        tabix(tmp)

        # SNPs + indels <50bp
        snp_out = f"{PREBENCH_OUTDIR}/{strain}.illumina.snps_indels.vcf.gz"
        run(
            f"bcftools view -v snps,indels "
            f"-i 'abs(strlen(ALT)-strlen(REF))<50' "
            f"-Oz -o {snp_out} {tmp}"
        )
        tabix(snp_out)

        # SVs >=50bp
        sv_out = f"{PREBENCH_OUTDIR}/{strain}.illumina.svs.vcf.gz"
        run(
            f"bcftools view "
            f"-i 'abs(strlen(ALT)-strlen(REF))>=50' "
            f"-Oz -o {sv_out} {tmp}"
        )
        tabix(sv_out)

        remove(tmp)

    # Copy paftools SNP/indel baselines
    log.info("  Copying paftools baselines...")
    for strain in STRAINS:
        src = f"{SVIM_BASE}/{strain}.April2018/{strain}.April2018.paftools.snps.vcf.gz"
        dst = f"{PREBENCH_OUTDIR}/{strain}.paftools.snps_indels.vcf.gz"
        if not Path(src).exists():
            log.warning(f"  WARNING: {src} not found, skipping")
            continue
        run(f"cp {src} {dst}")
        run(f"cp {src}.tbi {dst}.tbi")


# ---------------------------------------------------------------------------
# Step 6: Decompose MNPs and complex indels
# ---------------------------------------------------------------------------

def step_decompose() -> None:
    """Decompose MNPs and complex indels using rtg vcfdecompose."""
    log.info("=== Step 6: Decompose VCFs (rtg vcfdecompose) ===")
    ensure_dir(DECOMPOSE_OUTDIR)
    vcfs = sorted(Path(DECOMPOSE_INDIR).glob("*.vcf.gz"))
    if not vcfs:
        log.warning(f"No VCFs found in {DECOMPOSE_INDIR}")
        return
    for vcf in vcfs:
        basename = vcf.name.replace(".vcf.gz", "")
        out = f"{DECOMPOSE_OUTDIR}/{basename}.decomposed.vcf.gz"
        log.info(f"  -> {basename}")
        run(
            f"rtg vcfdecompose "
            f"--break-mnps --break-indels "
            f"-i {vcf} -o {out}"
        )
        tabix(out)


# ---------------------------------------------------------------------------
# Step 7: RTG vcfeval — short variant benchmarking
# ---------------------------------------------------------------------------

def step_vcfeval() -> None:
    """
    Benchmark short variants against paftools baselines using RTG vcfeval.
    Applies DP>5 filter to specified callsets, then filters all to 1/1 GTs.
    Runs across all BED regions and callsets.
    """
    log.info("=== Step 7: RTG vcfeval benchmarking ===")
    ensure_dir(VCFEVAL_OUTDIR)

    # Reheader paftools baselines to match strain sample names
    log.info("  Reheadering paftools baselines...")
    for strain in STRAINS:
        in_vcf  = f"{VCFEVAL_INDIR}/{strain}.paftools.snps_indels.decomposed.vcf.gz"
        out_vcf = f"{VCFEVAL_INDIR}/{strain}.paftools.snps_indels.decomposed.reheader.vcf.gz"
        if not Path(in_vcf).exists():
            log.warning(f"  WARNING: Missing {in_vcf}")
            continue
        if not Path(out_vcf).exists():
            rename_txt = "rename.txt"
            Path(rename_txt).write_text(f"sample\t{strain}\n")
            run(f"bcftools reheader -s {rename_txt} -o {out_vcf} {in_vcf}")
            tabix(out_vcf)
            Path(rename_txt).unlink(missing_ok=True)

    # Run vcfeval across all BED regions, callsets, and strains
    for bed in sorted(Path(BED_DIR).glob("*.bed")):
        bedname = bed.stem
        for caller, call_suffix in VCFEVAL_CALLSETS.items():
            caller_outdir = f"{VCFEVAL_OUTDIR}/{bedname}/{caller}"
            ensure_dir(caller_outdir)
            for strain in STRAINS:
                baseline       = f"{VCFEVAL_INDIR}/{strain}.paftools.snps_indels.decomposed.reheader.vcf.gz"
                calls          = f"{VCFEVAL_INDIR}/{strain}.{call_suffix}"
                sample_out     = f"{caller_outdir}/{strain}"

                if not Path(baseline).exists() or not Path(calls).exists():
                    log.warning(f"  WARNING: Missing files for {strain} ({caller})")
                    continue

                # DP>5 filter for specified callsets
                if caller in DP_FILTER_CALLSETS:
                    dp_filtered = f"{VCFEVAL_INDIR}/{strain}.{call_suffix.replace('.vcf.gz', '.dp5.vcf.gz')}"
                    if not Path(dp_filtered).exists():
                        log.info(f"  Filtering {strain} ({caller}) to DP>5...")
                        run(f"bcftools view -i 'FORMAT/DP>5' {calls} -Oz -o {dp_filtered}")
                        tabix(dp_filtered)
                    input_for_gt = dp_filtered
                    filtered_calls = f"{VCFEVAL_INDIR}/{strain}.{call_suffix.replace('.vcf.gz', '.dp5.1_1_only.vcf.gz')}"
                else:
                    input_for_gt   = calls
                    filtered_calls = f"{VCFEVAL_INDIR}/{strain}.{call_suffix.replace('.vcf.gz', '.1_1_only.vcf.gz')}"

                # Filter to homozygous ALT genotypes
                if not Path(filtered_calls).exists():
                    log.info(f"  Filtering {strain} ({caller}) to 1/1 genotypes...")
                    run(f"bcftools view -i 'GT=\"1/1\" || GT=\"1\"' {input_for_gt} -Oz -o {filtered_calls}")
                    tabix(filtered_calls)

                log.info(f"  -> {strain} | {caller} | {bedname}")
                run(
                    f"rtg vcfeval "
                    f"-m annotate "
                    f"--all-records "
                    f"--ref-overlap "
                    f"--no-roc "
                    f"--bed-regions {bed} "
                    f"-b {baseline} "
                    f"-c {filtered_calls} "
                    f"-t {REF_SDF} "
                    f"--sample {strain} "
                    f"-o {sample_out}"
                )


# ---------------------------------------------------------------------------
# Step 8: Truvari — SV benchmarking
# ---------------------------------------------------------------------------

def step_truvari() -> None:
    """
    Benchmark SVs (>=50bp) against SVIM-asm truth using Truvari bench.
    Filters to ALT-only calls (keeping mixed genotypes) before benchmarking.
    Runs across all BED regions and callsets.
    """
    log.info("=== Step 8: Truvari SV benchmarking ===")
    ensure_dir(TRUVARI_OUTDIR)

    for bed in sorted(Path(BED_DIR).glob("*.bed")):
        bedname = bed.stem
        for caller, call_suffix in TRUVARI_CALLSETS.items():
            caller_outdir = f"{TRUVARI_OUTDIR}/{bedname}/{caller}"
            ensure_dir(caller_outdir)
            for strain in STRAINS:
                baseline    = f"{TRUVARI_INDIR}/{strain}.svim.svs.vcf.gz"
                calls       = f"{TRUVARI_INDIR}/{strain}.{call_suffix}"
                alt_baseline = f"{TRUVARI_INDIR}/{strain}.svim.svs.alt_only.mixed.vcf.gz"
                alt_calls    = f"{TRUVARI_INDIR}/{strain}.{call_suffix.replace('.vcf.gz', '.alt_only.mixed.vcf.gz')}"
                sample_out  = f"{caller_outdir}/{strain}"

                if not Path(baseline).exists() or not Path(calls).exists():
                    log.warning(f"  WARNING: Missing VCFs for {strain} ({caller})")
                    continue

                # Filter to ALT-only, keeping mixed genotypes
                if not Path(alt_baseline).exists():
                    run(f"bcftools view -i 'GT!=\"ref\"' {baseline} -Oz -o {alt_baseline}")
                    tabix(alt_baseline)
                if not Path(alt_calls).exists():
                    run(f"bcftools view -i 'GT!=\"ref\"' {calls} -Oz -o {alt_calls}")
                    tabix(alt_calls)

                log.info(f"  -> {strain} | {caller} | {bedname}")
                run(
                    f"truvari bench "
                    f"-b {alt_baseline} "
                    f"-c {alt_calls} "
                    f"-o {sample_out} "
                    f"-f {REF_FASTA} "
                    f"-r 1000 -C 1000 "
                    f"-O 0.0 -p 0.0 -P 0.3 "
                    f"-s 50 -S 15 "
                    f"--sizemax 10000 "
                    f"--includebed {bed}"
                )


# ---------------------------------------------------------------------------
# Step 9: Plot evaluation results
# ---------------------------------------------------------------------------

def step_plot() -> None:
    """
    Call plot_evaluation.R to generate F-measure boxplots and summary CSVs
    for small variants (vcfeval) and SVs (truvari) across genomic regions.
    """
    log.info("=== Step 9: Plot evaluation results (plot_evaluation.R) ===")
    if not Path(PLOT_SCRIPT).exists():
        log.error(f"Plot script not found: {PLOT_SCRIPT}")
        return
    run(f"Rscript {PLOT_SCRIPT}")


# ---------------------------------------------------------------------------
# Step registry and CLI
# ---------------------------------------------------------------------------

ALL_STEPS = [
    "filter_delly_gatk",
    "filter_pangenome_vcfs",
    "filter_pangenome_vcfs_haploid",
    "recode_haploid_genotypes",
    "prebench",
    "decompose",
    "vcfeval",
    "truvari",
    "plot",
]

STEP_FUNCTIONS = {
    "filter_delly_gatk":           step_filter_delly_gatk,
    "filter_pangenome_vcfs":       step_filter_pangenome_vcfs,
    "filter_pangenome_vcfs_haploid": step_filter_pangenome_vcfs_haploid,
    "recode_haploid_genotypes":    step_recode_haploid_genotypes,
    "prebench":                    step_prebench,
    "decompose":                   step_decompose,
    "vcfeval":                     step_vcfeval,
    "truvari":                     step_truvari,
    "plot":                        step_plot,
}


def main():
    parser = argparse.ArgumentParser(
        description="Variant filtering, benchmarking, and evaluation pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--steps", nargs="+", choices=ALL_STEPS, default=None,
        help="Run only these steps (default: all steps in order)",
    )
    parser.add_argument(
        "--skip", nargs="+", choices=ALL_STEPS, default=[],
        help="Skip these steps",
    )
    args = parser.parse_args()

    steps_to_run = args.steps if args.steps else ALL_STEPS
    steps_to_run = [s for s in steps_to_run if s not in args.skip]

    log.info(f"Steps to run: {steps_to_run}")

    for step in steps_to_run:
        STEP_FUNCTIONS[step]()

    log.info("=== Pipeline 2 complete ===")


if __name__ == "__main__":
    main()
