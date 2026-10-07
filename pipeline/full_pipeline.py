#!/usr/bin/env python3
"""
Pf pangenome LOO validation pipeline -- Python port of full_pipeline.sh

DESIGN NOTE: external bioinformatics tools (gatk, bcftools, vg, delly,
manta, dysgu, truvari, rtg, bedtools, awk, samtools...) are still
invoked as shell command strings via subprocess. Reimplementing bcftools
filter expressions or awk one-liners in native Python would add real risk
of subtly changing behaviour for no real benefit. This script ports the
*orchestration* -- control flow, stage toggles, conda-env handling,
per-strain looping, idempotency checks -- to Python; the underlying tool
commands are functionally unchanged from full_pipeline.sh.

CHANGE FROM full_pipeline.sh: SV comparisons (run_truvari_all) now compare
against a single core-genome bed (CORE_BED below) instead of looping over
every category bed in BEDDIR. Short-variant vcfeval (run_vcfeval_all)
still loops over all of BEDDIR, unchanged.
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

############################################################
# CONFIG
# No hardcoded values live here anymore -- every default lives once, in
# build_arg_parser() below, and apply_args() populates these names before
# main() runs. Running with no CLI arguments reproduces the original
# full_pipeline.sh values exactly; every one of them can be overridden.
############################################################
STRAINS = None
SEQ_FILE = None
READS_DIR = None
LINEAR_GVCF_DIR = None
LOO_BASE = None
COMPARISON_DIR = None

CACTUS_VENV = None
REF = None
REF_GATK = None
REF_SDF = None
VT_BIN = None
GFF = None

BEDDIR = None
VSA_BED = None
TRF_BED = None
CONF_BED_DIR = None
CORE_BED = None

OUTROOT = None

VSA_PAD = None
VT_MIN_FZ_RL = None
BOUNDARY_BUFFER = None
DENSITY_WINDOW = None
DENSITY_PERCENTILE = None
MIN_BP = None

RUN_MAPPING = None
RUN_SMALL_VAR = None
RUN_SV = None
RUN_SV_EXTRA = None
RUN_SV_DYSGU = None

# Set by build_shared_beds()
GENOME_FILE = None

############################################################
# SHELL / ENV HELPERS
############################################################

def sh(cmd: str, cwd=None, check=True, capture=False):
    """Run a shell command string via bash -c (set -euo pipefail applied),
    mirroring the original script's direct shell calls (pipes, redirection,
    awk etc. kept as-is)."""
    full_cmd = "set -euo pipefail\n" + cmd
    kwargs = dict(cwd=str(cwd) if cwd else None, check=check)
    if capture:
        kwargs.update(capture_output=True, text=True)
    result = subprocess.run(["bash", "-c", full_cmd], **kwargs)
    return result.stdout if capture else None


def conda_sh(env: str, cmd: str, cwd=None, check=True, capture=False):
    """Equivalent of the bash conda_activate() helper: source conda.sh,
    activate ENV (with -u relaxed around activation itself, matching the
    original workaround for conda's activation script referencing unset
    vars), then run CMD in that shell."""
    full = (
        'source "$(conda info --base)/etc/profile.d/conda.sh"\n'
        'set +u\n'
        f'conda activate {env}\n'
        'set -u\n'
        f'{cmd}'
    )
    return sh(full, cwd=cwd, check=check, capture=capture)


def is_file(p) -> bool:
    return Path(p).is_file()


def is_nonempty(p) -> bool:
    p = Path(p)
    return p.is_file() and p.stat().st_size > 0


def sum_bed_bp(bed_path) -> int:
    """awk '{s+=$3-$2}END{print s+0}' equivalent."""
    total = 0
    with open(bed_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            total += int(parts[2]) - int(parts[1])
    return total


def rmtree_if_exists(p: Path):
    if p.exists():
        shutil.rmtree(p)


def unlink_if_exists(p: Path):
    """p.unlink(missing_ok=True) equivalent, compatible with Python < 3.8."""
    try:
        Path(p).unlink()
    except FileNotFoundError:
        pass


############################################################
# STEP 0: shared, genome-wide beds (build once)
############################################################
def build_shared_beds():
    global GENOME_FILE
    if not is_file(VSA_BED):
        print("=== building VSA bed ===")
        tmp = VSA_BED.parent / "vsa_tmp"
        tmp.mkdir(parents=True, exist_ok=True)
        sh("""awk -F'\\t' '$3=="gene"{print}' %(gff)s | \
            grep -P 'Name=(VAR|VAR-like|VAR2CSA|RIF|RIFA|MC-2TM|SURF)' | \
            awk -F'\\t' 'BEGIN{OFS="\\t"}{print $1,$4-1,$5}' > %(out)s""" % {
            "gff": GFF, "out": tmp / "vsa_genes.bed"})
        sh("""awk -F'\\t' '$3=="repeat_region"{print $1"\\t"($4-1)"\\t"$5}' %(gff)s > %(out)s""" % {
            "gff": GFF, "out": tmp / "repeat_regions.bed"})
        sh("""awk -F'\\t' '$3=="centromere"{print $1"\\t"($4-1)"\\t"$5}' %(gff)s > %(out)s""" % {
            "gff": GFF, "out": tmp / "centromeres.bed"})
        sh("""cat %(tmp)s/*.bed | sort -k1,1 -k2,2n > %(out)s""" % {
            "tmp": tmp, "out": tmp / "combined_raw.bed"})
        conda_sh("truvari", """awk -v p=%(pad)s 'BEGIN{OFS="\\t"}{s=$2-p; if(s<0)s=0; print $1,s,$3+p}' %(inp)s | \
            sort -k1,1 -k2,2n | bedtools merge -i - > %(out)s""" % {
            "pad": VSA_PAD, "inp": tmp / "combined_raw.bed", "out": VSA_BED})

    if not is_file(TRF_BED):
        sys.exit("ERROR: TRF bed not found")

    GENOME_FILE = OUTROOT / "genome.txt"
    if not is_file(GENOME_FILE):
        conda_sh("fastq2matrix", """samtools faidx %(ref)s""" % {"ref": REF})
        sh("""cut -f1,2 %(fai)s > %(out)s""" % {"fai": str(REF) + ".fai", "out": GENOME_FILE})


############################################################
# STEP 1: LOO graph build + read mapping + vg call + GATK-on-surject
############################################################
def run_loo(strain: str):
    work = LOO_BASE / f"loo_{strain}"
    work.mkdir(parents=True, exist_ok=True)

    print(f"=== [{strain}] LOO seq file ===")
    sh("""grep -v "^%(strain)s\\b" %(seq)s > pf3k_seq_no_%(strain)s.txt""" % {
        "strain": strain, "seq": SEQ_FILE}, cwd=work)

    graph_name = f"PfPan_no_{strain}"
    out_dir = work / f"{graph_name}_pan"
    gbz = out_dir / f"{graph_name}.gbz"
    snarls = out_dir / f"{graph_name}.snarls"

    print(f"=== [{strain}] cactus-pangenome (build if missing) ===")
    if not is_file(gbz):
        job_store = work / f"{graph_name}_js"
        # --restart resumes an existing Toil job store after a crash; it
        # errors out if the job store doesn't exist yet (e.g. after a
        # from-scratch rerun where the strain's whole loo_ dir was wiped),
        # so only pass it when there's actually something to resume.
        restart_flag = "--restart" if job_store.exists() else ""
        sh("""set +u
            source %(venv)s
            set -u
            cactus-pangenome %(js)s %(seqfile)s \
            --outDir %(outdir)s --outName %(name)s --reference Pf3D7 \
            --filter 2 --haplo --giraffe clip filter \
            --gbz clip filter full --gfa clip full \
            --vcf --vcfReference Pf3D7 \
            --logFile %(log)s \
            --workDir %(work)s --consCores 8 --mgMemory 128Gi \
            %(restart)s
            deactivate""" % {
            "venv": CACTUS_VENV, "js": job_store,
            "seqfile": work / f"pf3k_seq_no_{strain}.txt",
            "outdir": out_dir, "name": graph_name,
            "log": work / f"{graph_name}.log", "work": work,
            "restart": restart_flag,
        })
    if not is_file(gbz):
        print(f"ERROR: GBZ build failed for {strain}")
        return

    print(f"=== [{strain}] vg giraffe + surject ===")
    fq1 = READS_DIR / f"{strain}_1.trimmed.fastq.gz"
    fq2 = READS_DIR / f"{strain}_2.trimmed.fastq.gz"
    gaf = work / f"{strain}.loo.gaf.gz"
    sort_bam = work / f"{strain}.loo.sort.bam"

    if not is_nonempty(gaf):
        conda_sh("cactus", """vg giraffe -p -t 16 -Z %(gbz)s -f %(fq1)s -f %(fq2)s -o gaf | bgzip > %(gaf)s""" % {
            "gbz": gbz, "fq1": fq1, "fq2": fq2, "gaf": gaf})
    conda_sh("cactus", """vg paths -x %(gbz)s -L | grep "^Pf3D7" > %(out)s""" % {
        "gbz": gbz, "out": work / "Pf3D7.paths.txt"})

    if not is_nonempty(sort_bam):
        conda_sh("cactus", """vg surject -x %(gbz)s -G %(gaf)s --interleaved \
            -F %(paths)s -b \
            -N %(strain)s -R "ID:1 LB:lib1 SM:%(strain)s PL:illumina PU:unit1" \
            | samtools reheader -c 'sed s/Pf3D7#0#//g' - > %(bam)s""" % {
            "gbz": gbz, "gaf": gaf, "paths": work / "Pf3D7.paths.txt",
            "strain": strain, "bam": work / f"{strain}.loo.bam"})
        conda_sh("cactus", """samtools sort %(bam)s -O BAM -o %(sorted)s --threads 8""" % {
            "bam": work / f"{strain}.loo.bam", "sorted": sort_bam})
        conda_sh("cactus", """samtools index %(sorted)s -@ 8""" % {"sorted": sort_bam})

    print(f"=== [{strain}] vg pack + vg call (haploid) ===")
    pack = work / f"{strain}.loo.pack"
    hap_vcf = work / f"{strain}.loo.SV.call_haploid_.vcf.gz"
    if not is_nonempty(pack):
        conda_sh("cactus", """vg pack -x %(gbz)s -Q5 -a %(gaf)s -o %(pack)s""" % {
            "gbz": gbz, "gaf": gaf, "pack": pack})
    if not is_nonempty(hap_vcf):
        conda_sh("cactus", """vg call %(gbz)s -r %(snarls)s -k %(pack)s -t 16 --ploidy 1 \
            -s %(strain)s -S Pf3D7 -az | bgzip > %(vcf)s""" % {
            "gbz": gbz, "snarls": snarls, "pack": pack, "strain": strain, "vcf": hap_vcf})

    print(f"=== [{strain}] fix contig naming on raw vg call output ===")
    fix_contigs(hap_vcf, f"loo_{strain}")

    print(f"=== [{strain}] GATK on LOO surjected BAM ===")
    sample = f"{strain}.loo.gatk_surject"
    gvcf = work / f"{sample}.g.vcf.gz"
    if not is_nonempty(gvcf):
        conda_sh("fastq2matrix", """gatk HaplotypeCaller -I %(bam)s -R %(ref)s -O %(gvcf)s -ERC GVCF""" % {
            "bam": sort_bam, "ref": REF_GATK, "gvcf": gvcf})


############################################################
# fix_contigs: strip Pf3D7#0# prefix from vg call output (header + body)
############################################################
def fix_contigs(vcf: Path, tag: str):
    work = vcf.parent
    vcf_stem = str(vcf)
    assert vcf_stem.endswith(".vcf.gz")
    vcf_stem = vcf_stem[: -len(".vcf.gz")]
    renamed = Path(vcf_stem + ".renamed.vcf.gz")

    conda_sh("cactus", """bcftools index -f -t %(vcf)s""" % {"vcf": vcf}, cwd=work)
    conda_sh("cactus", """bcftools view -h %(vcf)s | grep "^##contig" | \
        sed -E 's/##contig=<ID=([^,]+),.*/\\1/' | \
        awk '{new=$0; gsub("Pf3D7#0#","",new); print $0"\\t"new}' > rename_map_%(tag)s.txt""" % {
        "vcf": vcf, "tag": tag}, cwd=work)
    conda_sh("cactus", """bcftools annotate --rename-chrs rename_map_%(tag)s.txt %(vcf)s -Oz -o body_renamed_%(tag)s.vcf.gz""" % {
        "tag": tag, "vcf": vcf}, cwd=work)
    conda_sh("cactus", """bcftools view -h body_renamed_%(tag)s.vcf.gz | sed 's/Pf3D7#0#//g' > header_fixed_%(tag)s.txt""" % {
        "tag": tag}, cwd=work)
    conda_sh("cactus", """bcftools reheader -h header_fixed_%(tag)s.txt body_renamed_%(tag)s.vcf.gz -o %(renamed)s""" % {
        "tag": tag, "renamed": renamed}, cwd=work)
    conda_sh("cactus", """bcftools index -f -t %(renamed)s""" % {"renamed": renamed}, cwd=work)

    renamed.rename(vcf)
    Path(str(renamed) + ".tbi").rename(Path(str(vcf) + ".tbi"))
    for f in [work / f"body_renamed_{tag}.vcf.gz", work / f"header_fixed_{tag}.txt", work / f"rename_map_{tag}.txt"]:
        unlink_if_exists(f)


############################################################
# GATK filter chain, reused for linear and for LOO-surject GVCFs
############################################################
def filter_gatk(gvcf: Path, sample: str, workdir: Path, final_out: Path):
    workdir.mkdir(parents=True, exist_ok=True)

    conda_sh("fastq2matrix", """gatk GenotypeGVCFs -R %(ref)s -V %(gvcf)s -O %(sample)s.raw.vcf.gz""" % {
        "ref": REF_GATK, "gvcf": gvcf, "sample": sample}, cwd=workdir)
    conda_sh("fastq2matrix", """tabix -f -p vcf %(sample)s.raw.vcf.gz""" % {"sample": sample}, cwd=workdir)

    conda_sh("fastq2matrix", """gatk SelectVariants -R %(ref)s -V %(sample)s.raw.vcf.gz --select-type-to-include SNP -O %(sample)s.snps.vcf.gz""" % {
        "ref": REF_GATK, "sample": sample}, cwd=workdir)
    conda_sh("fastq2matrix", """tabix -f -p vcf %(sample)s.snps.vcf.gz""" % {"sample": sample}, cwd=workdir)
    conda_sh("fastq2matrix", """gatk SelectVariants -R %(ref)s -V %(sample)s.raw.vcf.gz --select-type-to-include INDEL -O %(sample)s.indels.vcf.gz""" % {
        "ref": REF_GATK, "sample": sample}, cwd=workdir)
    conda_sh("fastq2matrix", """tabix -f -p vcf %(sample)s.indels.vcf.gz""" % {"sample": sample}, cwd=workdir)

    conda_sh("fastq2matrix", """gatk VariantFiltration -R %(ref)s -V %(sample)s.snps.vcf.gz \
        --filter-expression "QD < 2.0"             --filter-name "QD2" \
        --filter-expression "QUAL < 30.0"           --filter-name "QUAL30" \
        --filter-expression "SOR > 4.0"             --filter-name "SOR4" \
        --filter-expression "FS > 60.0"             --filter-name "FS60" \
        --filter-expression "MQ < 40.0"             --filter-name "MQ40" \
        --filter-expression "MQRankSum < -15.0"     --filter-name "MQRankSum-15" \
        --filter-expression "ReadPosRankSum < -5.0" --filter-name "ReadPosRankSum-5" \
        -O %(sample)s.snps.filtered.vcf.gz""" % {"ref": REF_GATK, "sample": sample}, cwd=workdir)
    conda_sh("fastq2matrix", """tabix -f -p vcf %(sample)s.snps.filtered.vcf.gz""" % {"sample": sample}, cwd=workdir)

    conda_sh("fastq2matrix", """gatk VariantFiltration -R %(ref)s -V %(sample)s.indels.vcf.gz \
        --filter-expression "QD < 2.0"               --filter-name "QD2" \
        --filter-expression "FS > 200.0"             --filter-name "FS200" \
        --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
        -O %(sample)s.indels.filtered.vcf.gz""" % {"ref": REF_GATK, "sample": sample}, cwd=workdir)
    conda_sh("fastq2matrix", """tabix -f -p vcf %(sample)s.indels.filtered.vcf.gz""" % {"sample": sample}, cwd=workdir)

    conda_sh("fastq2matrix", """bcftools concat -a %(sample)s.snps.filtered.vcf.gz %(sample)s.indels.filtered.vcf.gz -Oz -o %(sample)s.filtered.vcf.gz""" % {
        "sample": sample}, cwd=workdir)
    conda_sh("fastq2matrix", """tabix -f -p vcf %(sample)s.filtered.vcf.gz""" % {"sample": sample}, cwd=workdir)
    conda_sh("fastq2matrix", """gatk SelectVariants -R %(ref)s -V %(sample)s.filtered.vcf.gz --exclude-filtered -O %(sample)s.PASS.vcf.gz""" % {
        "ref": REF_GATK, "sample": sample}, cwd=workdir)
    conda_sh("fastq2matrix", """tabix -f -p vcf %(sample)s.PASS.vcf.gz""" % {"sample": sample}, cwd=workdir)
    conda_sh("fastq2matrix", """bcftools norm -f %(ref)s -m -both %(sample)s.PASS.vcf.gz -Oz -o %(sample)s.norm.vcf.gz""" % {
        "ref": REF_GATK, "sample": sample}, cwd=workdir)
    conda_sh("fastq2matrix", """tabix -f -p vcf %(sample)s.norm.vcf.gz""" % {"sample": sample}, cwd=workdir)

    conda_sh("fastq2matrix", """bcftools view -i 'TYPE="snp" || (TYPE="indel" && abs(strlen(ALT)-strlen(REF))<50)' \
        %(sample)s.norm.vcf.gz -Oz -o %(sample)s.shortvars.lt50bp.PASS.vcf.gz""" % {"sample": sample}, cwd=workdir)
    conda_sh("fastq2matrix", """tabix -f -p vcf %(sample)s.shortvars.lt50bp.PASS.vcf.gz""" % {"sample": sample}, cwd=workdir)

    conda_sh("fastq2matrix", """bcftools view %(sample)s.shortvars.lt50bp.PASS.vcf.gz | setGT.py --fraction 0.7 | \
        bcftools view -O z -c 1 -o %(sample)s.shortvars.lt50bp.PASS.GT.vcf.gz""" % {"sample": sample}, cwd=workdir)
    conda_sh("fastq2matrix", """tabix -f -p vcf %(sample)s.shortvars.lt50bp.PASS.GT.vcf.gz""" % {"sample": sample}, cwd=workdir)

    final_out.parent.mkdir(parents=True, exist_ok=True)
    unlink_if_exists(final_out)
    unlink_if_exists(Path(str(final_out) + ".tbi"))
    conda_sh("truvari", """rtg vcfdecompose --break-mnps --break-indels -i %(sample)s.shortvars.lt50bp.PASS.GT.vcf.gz -o %(out)s""" % {
        "sample": sample, "out": final_out}, cwd=workdir)
    conda_sh("truvari", """tabix -f -p vcf %(out)s""" % {"out": final_out}, cwd=workdir)


############################################################
# STEP 2: linear GATK (raw gvcfs -> filtered), reused across strains
############################################################
def run_gatk_linear(strain: str):
    work = OUTROOT / strain / "gatk_linear"
    work.mkdir(parents=True, exist_ok=True)
    gvcf = LINEAR_GVCF_DIR / f"{strain}.g.vcf.gz"
    if not is_file(gvcf):
        print(f"WARNING: no linear gvcf for {strain}")
        return
    filter_gatk(gvcf, f"{strain}.gatk_linear", work,
                work / f"{strain}.GATK_linear.shortvars.lt50bp.PASS.GT.decomposed.vcf.gz")


############################################################
# STEP 3: GATK-on-LOO-surject filtering (uses gvcf from run_loo)
############################################################
def run_gatk_surject_loo(strain: str):
    loo_work = LOO_BASE / f"loo_{strain}"
    work = OUTROOT / strain / "gatk_surject_loo"
    work.mkdir(parents=True, exist_ok=True)
    gvcf = loo_work / f"{strain}.loo.gatk_surject.g.vcf.gz"
    if not is_file(gvcf):
        print(f"WARNING: no LOO surject gvcf for {strain}, run_loo first")
        return
    filter_gatk(gvcf, f"{strain}.gatk_surject_loo", work,
                work / f"{strain}.GATK_surject_loo.shortvars.lt50bp.PASS.GT.decomposed.vcf.gz")


############################################################
# STEP 4: vg_hap_loo filtering (mirrors prepare_variants_3.sh)
############################################################
def run_vg_hap_loo_filter(strain: str):
    work = LOO_BASE / f"loo_{strain}"
    vcf = f"{strain}.loo.SV.call_haploid_.vcf.gz"
    sample = f"{strain}.loo"

    conda_sh("cactus", """bcftools norm -f %(ref)s -m -both %(vcf)s -Oz -o %(sample)s.pan.norm.unsorted.hap.vcf.gz""" % {
        "ref": REF, "vcf": vcf, "sample": sample}, cwd=work)
    conda_sh("cactus", """bcftools sort %(sample)s.pan.norm.unsorted.hap.vcf.gz -Oz -o %(sample)s.pan.norm.hap.vcf.gz""" % {
        "sample": sample}, cwd=work)
    conda_sh("cactus", """tabix -f -p vcf %(sample)s.pan.norm.hap.vcf.gz""" % {"sample": sample}, cwd=work)

    conda_sh("cactus", """bcftools view -f PASS -v snps %(sample)s.pan.norm.hap.vcf.gz -Oz -o %(sample)s.pan.snps.PASS.hap.vcf.gz""" % {
        "sample": sample}, cwd=work)
    conda_sh("cactus", """tabix -f -p vcf %(sample)s.pan.snps.PASS.hap.vcf.gz""" % {"sample": sample}, cwd=work)
    conda_sh("cactus", """bcftools view -f PASS -v indels -i 'abs(strlen(ALT)-strlen(REF))<50' \
        %(sample)s.pan.norm.hap.vcf.gz -Oz -o %(sample)s.pan.indels.lt50bp.PASS.hap.vcf.gz""" % {
        "sample": sample}, cwd=work)
    conda_sh("cactus", """tabix -f -p vcf %(sample)s.pan.indels.lt50bp.PASS.hap.vcf.gz""" % {"sample": sample}, cwd=work)
    conda_sh("cactus", """bcftools concat -a %(sample)s.pan.snps.PASS.hap.vcf.gz %(sample)s.pan.indels.lt50bp.PASS.hap.vcf.gz \
        -Oz -o %(sample)s.pan.shortvars.lt50bp.PASS.hap.vcf.gz""" % {"sample": sample}, cwd=work)
    conda_sh("cactus", """tabix -f -p vcf %(sample)s.pan.shortvars.lt50bp.PASS.hap.vcf.gz""" % {"sample": sample}, cwd=work)
    conda_sh("cactus", """bcftools view -f PASS -i 'abs(strlen(ALT)-strlen(REF))>=50' \
        %(sample)s.pan.norm.hap.vcf.gz -Oz -o %(sample)s.pan.SV.ge50bp.PASS.hap.vcf.gz""" % {"sample": sample}, cwd=work)
    conda_sh("cactus", """tabix -f -p vcf %(sample)s.pan.SV.ge50bp.PASS.hap.vcf.gz""" % {"sample": sample}, cwd=work)

    conda_sh("cactus", """bcftools +setGT %(sample)s.pan.shortvars.lt50bp.PASS.hap.vcf.gz -- -t q -n 'c:0/0' -i 'GT="0"' \
        | bcftools +setGT -- -t q -n 'c:1/1' -i 'GT="1"' | bgzip -c > tmp.vcf.gz""" % {"sample": sample}, cwd=work)
    (work / "tmp.vcf.gz").rename(work / f"{sample}.pan.shortvars.lt50bp.PASS.hap.vcf.gz")
    conda_sh("cactus", """tabix -f -p vcf %(sample)s.pan.shortvars.lt50bp.PASS.hap.vcf.gz""" % {"sample": sample}, cwd=work)

    conda_sh("cactus", """bcftools view -i 'FORMAT/DP>5' %(sample)s.pan.shortvars.lt50bp.PASS.hap.vcf.gz -Oz -o %(sample)s.pan.shortvars.lt50bp.PASS.hap.dp5.vcf.gz""" % {
        "sample": sample}, cwd=work)
    conda_sh("cactus", """tabix -f -p vcf %(sample)s.pan.shortvars.lt50bp.PASS.hap.dp5.vcf.gz""" % {"sample": sample}, cwd=work)
    conda_sh("cactus", """bcftools view -i 'GT="1/1" || GT="1"' %(sample)s.pan.shortvars.lt50bp.PASS.hap.dp5.vcf.gz \
        -Oz -o %(sample)s.pan.shortvars.lt50bp.PASS.hap.dp5.1_1_only.vcf.gz""" % {"sample": sample}, cwd=work)
    conda_sh("cactus", """tabix -f -p vcf %(sample)s.pan.shortvars.lt50bp.PASS.hap.dp5.1_1_only.vcf.gz""" % {"sample": sample}, cwd=work)

    final_loo = work / f"{sample}.pan.shortvars.lt50bp.PASS.hap.dp5.1_1_only.decomposed.vcf.gz"
    unlink_if_exists(final_loo)
    unlink_if_exists(Path(str(final_loo) + ".tbi"))
    conda_sh("truvari", """rtg vcfdecompose --break-mnps --break-indels \
        -i %(sample)s.pan.shortvars.lt50bp.PASS.hap.dp5.1_1_only.vcf.gz \
        -o %(out)s""" % {"sample": sample, "out": final_loo}, cwd=work)
    conda_sh("truvari", """tabix -f -p vcf %(out)s""" % {"out": final_loo}, cwd=work)


############################################################
# STEP 5: delly on LOO surjected BAM (SV caller)
############################################################
def run_delly_loo(strain: str):
    work = LOO_BASE / f"loo_{strain}"
    conda_sh("delly", """delly call -g %(ref)s %(strain)s.loo.sort.bam -o %(strain)s.loo.delly.bcf""" % {
        "ref": REF, "strain": strain}, cwd=work)
    conda_sh("fastq2matrix", """bcftools view %(strain)s.loo.delly.bcf -f PASS \
        -i 'INFO/SVTYPE!="" && abs(INFO/SVLEN)>=50' \
        -Oz -o %(strain)s.loo.delly.SV.ge50bp.PASS.vcf.gz""" % {"strain": strain}, cwd=work)
    conda_sh("fastq2matrix", """tabix -f -p vcf %(strain)s.loo.delly.SV.ge50bp.PASS.vcf.gz""" % {"strain": strain}, cwd=work)


############################################################
# STEP 5b: delly on the existing linear (BWA/BQSR) BAM -- true linear baseline
############################################################
def run_delly_linear(strain: str):
    work = OUTROOT / strain / "delly_linear"
    work.mkdir(parents=True, exist_ok=True)
    bam = READS_DIR / f"{strain}.bqsr.cram"
    if not is_file(bam):
        print(f"WARNING: no linear bqsr cram for {strain} at {bam}")
        return

    if not is_file(Path(str(bam) + ".crai")):
        conda_sh("fastq2matrix", """samtools index %(bam)s""" % {"bam": bam})

    bcf = work / f"{strain}.delly_linear.bcf"
    if not is_nonempty(bcf):
        conda_sh("delly", """delly call -g %(ref)s %(bam)s -o %(strain)s.delly_linear.bcf""" % {
            "ref": REF, "bam": bam, "strain": strain}, cwd=work)

    conda_sh("fastq2matrix", """bcftools view %(strain)s.delly_linear.bcf -f PASS \
        -i 'INFO/SVTYPE!="" && abs(INFO/SVLEN)>=50' \
        -Oz -o %(strain)s.delly_linear.SV.ge50bp.PASS.vcf.gz""" % {"strain": strain}, cwd=work)
    conda_sh("fastq2matrix", """tabix -f -p vcf %(strain)s.delly_linear.SV.ge50bp.PASS.vcf.gz""" % {"strain": strain}, cwd=work)


############################################################
# STEP 5c/5d: manta (loo + linear)
# NOTE: manta names its output diploidSV.vcf.gz regardless of organism
# ploidy -- that's just its fixed filename, not a ploidy assumption baked
# into variant *detection* (assembly/evidence-based). The reported GT does
# assume diploidy though, so downstream we treat calls as "present" via
# GT!="ref" rather than trusting het/hom calls (see run_truvari_all).
############################################################
def run_manta_loo(strain: str):
    work = LOO_BASE / f"loo_{strain}"
    bam = work / f"{strain}.loo.sort.bam"
    rundir = work / "manta_loo"
    if not is_file(bam):
        print(f"WARNING: no LOO surject bam for {strain}, run_loo first")
        return

    diploid = rundir / "results" / "variants" / "diploidSV.vcf.gz"
    if not is_nonempty(diploid):
        rmtree_if_exists(rundir)
        conda_sh("manta", """configManta.py --bam %(bam)s --referenceFasta %(ref)s --runDir %(rundir)s""" % {
            "bam": bam, "ref": REF, "rundir": rundir})
        conda_sh("manta", """%(rundir)s/runWorkflow.py -m local -j 16""" % {"rundir": rundir})

    out = work / f"{strain}.loo.manta.SV.ge50bp.PASS.vcf.gz"
    conda_sh("fastq2matrix", """bcftools view %(diploid)s -f PASS \
        -i 'INFO/SVTYPE!="" && (INFO/SVTYPE="BND" || abs(INFO/SVLEN)>=50)' \
        -Oz -o %(out)s""" % {"diploid": diploid, "out": out})
    conda_sh("fastq2matrix", """tabix -f -p vcf %(out)s""" % {"out": out})


def run_manta_linear(strain: str):
    work = OUTROOT / strain / "manta_linear"
    work.mkdir(parents=True, exist_ok=True)
    bam = READS_DIR / f"{strain}.bqsr.cram"
    rundir = work / "run"
    if not is_file(bam):
        print(f"WARNING: no linear bqsr cram for {strain} at {bam}")
        return

    conda_sh("fastq2matrix", """[[ -f %(bam)s.crai ]] || samtools index %(bam)s""" % {"bam": bam})

    diploid = rundir / "results" / "variants" / "diploidSV.vcf.gz"
    if not is_nonempty(diploid):
        rmtree_if_exists(rundir)
        conda_sh("manta", """configManta.py --bam %(bam)s --referenceFasta %(ref)s --runDir %(rundir)s""" % {
            "bam": bam, "ref": REF, "rundir": rundir})
        conda_sh("manta", """%(rundir)s/runWorkflow.py -m local -j 16""" % {"rundir": rundir})

    out = work / f"{strain}.manta_linear.SV.ge50bp.PASS.vcf.gz"
    conda_sh("fastq2matrix", """bcftools view %(diploid)s -f PASS \
        -i 'INFO/SVTYPE!="" && (INFO/SVTYPE="BND" || abs(INFO/SVLEN)>=50)' \
        -Oz -o %(out)s""" % {"diploid": diploid, "out": out})
    conda_sh("fastq2matrix", """tabix -f -p vcf %(out)s""" % {"out": out})


############################################################
# STEP 5g/5h: dysgu (loo + linear)
# NOTE: dysgu's ML component classifies alignment-signature features from
# this sample's own reads (split-read/discordant-pair/soft-clip patterns);
# it is not a model pretrained on a human (or any other) reference genome,
# so unlike DeepVariant it's not organism-locked by design.
############################################################
def run_dysgu_loo(strain: str):
    work = LOO_BASE / f"loo_{strain}"
    bam = work / f"{strain}.loo.sort.bam"
    tmpdir = work / "dysgu_loo_tmp"
    raw = work / f"{strain}.loo.dysgu.vcf"
    if not is_file(bam):
        print(f"WARNING: no LOO surject bam for {strain}, run_loo first")
        return

    if not is_nonempty(raw):
        rmtree_if_exists(tmpdir)
        conda_sh("dysgu", """dysgu run -p 4 -x %(ref)s %(tmp)s %(bam)s > %(raw)s""" % {
            "ref": REF, "tmp": tmpdir, "bam": bam, "raw": raw})

    out = work / f"{strain}.loo.dysgu.SV.ge50bp.PASS.vcf.gz"
    conda_sh("fastq2matrix", """bcftools view %(raw)s -f PASS \
        -i 'INFO/SVTYPE!="" && abs(INFO/SVLEN)>=50' \
        -Oz -o %(out)s""" % {"raw": raw, "out": out})
    conda_sh("fastq2matrix", """tabix -f -p vcf %(out)s""" % {"out": out})


############################################################
# STEP 5i: matched-linear BAM -- BWA-MEM on the SAME trimmed fastqs used by
# giraffe, with no BQSR, so the only variable left between this and the
# loo/surject BAM is the mapping approach itself (graph+surject vs direct
# linear alignment). The pre-existing READS_DIR/{strain}.bqsr.cram used
# elsewhere in this pipeline comes from a separate historical pipeline with
# its own trimming/aligner/BQSR choices, which confounds a loo-vs-linear
# comparison -- this function removes that confound for dysgu specifically.
############################################################
def run_bwa_matched_linear(strain: str):
    work = OUTROOT / strain / "matched_linear"
    work.mkdir(parents=True, exist_ok=True)
    fq1 = READS_DIR / f"{strain}_1.trimmed.fastq.gz"
    fq2 = READS_DIR / f"{strain}_2.trimmed.fastq.gz"
    bam = work / f"{strain}.matched_linear.sort.bam"
    bam_tmp = work / f"{strain}.matched_linear.sort.bam.tmp"

    if not is_file(fq1) or not is_file(fq2):
        print(f"WARNING: no trimmed fastqs for {strain} at {fq1}")
        return None

    if not is_file(Path(str(REF) + ".bwt")):
        conda_sh("fastq2matrix", """bwa index %(ref)s""" % {"ref": REF})

    if not is_nonempty(bam):
        # write + verify under a .tmp name first: a job killed mid-sort
        # (OOM, walltime) then just leaves a stale .tmp file behind rather
        # than a corrupt-but-nonempty bam that would fool is_nonempty()
        # into skipping the rebuild on the next run
        unlink_if_exists(bam_tmp)
        conda_sh("fastq2matrix", """bwa mem -t 16 -R "@RG\\tID:1\\tLB:lib1\\tSM:%(strain)s\\tPL:illumina\\tPU:unit1" \
            %(ref)s %(fq1)s %(fq2)s | samtools sort -O BAM -o %(bam_tmp)s --threads 8""" % {
            "strain": strain, "ref": REF, "fq1": fq1, "fq2": fq2, "bam_tmp": bam_tmp})
        conda_sh("fastq2matrix", """samtools quickcheck -v %(bam_tmp)s""" % {"bam_tmp": bam_tmp})
        bam_tmp.rename(bam)
        conda_sh("fastq2matrix", """samtools index %(bam)s -@ 8""" % {"bam": bam})

    return bam


def run_dysgu_linear(strain: str):
    work = OUTROOT / strain / "dysgu_linear"
    work.mkdir(parents=True, exist_ok=True)
    bam = run_bwa_matched_linear(strain)
    tmpdir = work / "tmp"
    raw = work / f"{strain}.dysgu_linear.vcf"
    if bam is None or not is_file(bam):
        print(f"WARNING: no matched-linear bam for {strain}")
        return

    if not is_nonempty(raw):
        rmtree_if_exists(tmpdir)
        conda_sh("dysgu", """dysgu run -p 4 -x %(ref)s %(tmp)s %(bam)s > %(raw)s""" % {
            "ref": REF, "tmp": tmpdir, "bam": bam, "raw": raw})

    out = work / f"{strain}.dysgu_linear.SV.ge50bp.PASS.vcf.gz"
    conda_sh("fastq2matrix", """bcftools view %(raw)s -f PASS \
        -i 'INFO/SVTYPE!="" && abs(INFO/SVLEN)>=50' \
        -Oz -o %(out)s""" % {"raw": raw, "out": out})
    conda_sh("fastq2matrix", """tabix -f -p vcf %(out)s""" % {"out": out})


############################################################
# HELPER: reheader raw paftools truth VCF to use the strain name as sample
############################################################
def reheader_truth(strain: str) -> Path:
    raw = COMPARISON_DIR / f"{strain}.paftools.snps_indels.vcf.gz"
    strain_dir = OUTROOT / strain
    strain_dir.mkdir(parents=True, exist_ok=True)
    renamed = strain_dir / f"{strain}.paftools.snps_indels.reheader.vcf.gz"

    if not is_file(renamed):
        current_name = conda_sh("truvari", """bcftools view -h %(raw)s | tail -1 | awk -F'\\t' '{print $NF}'""" % {
            "raw": raw}, capture=True).strip()
        rename_txt = strain_dir / "rename.txt"
        rename_txt.write_text(f"{current_name}\t{strain}\n")
        conda_sh("truvari", """bcftools reheader -s %(rename)s -o %(out)s %(raw)s""" % {
            "rename": rename_txt, "out": renamed, "raw": raw})
        conda_sh("truvari", """tabix -f -p vcf %(out)s""" % {"out": renamed})
        unlink_if_exists(rename_txt)
    return renamed


############################################################
# STEP 6: confidence tiers (base + VSA + vt-union [short var] / TRF [SV])
############################################################
def build_confidence(strain: str):
    bed_out = OUTROOT / strain / "beds"
    bed_out.mkdir(parents=True, exist_ok=True)

    truth_vcf = reheader_truth(strain)
    conf_raw = CONF_BED_DIR / f"{strain}.April2018" / f"{strain}.April2018.confident.bed"
    if not is_file(conf_raw):
        print(f"WARNING: no raw confidence bed for {strain}")
        return

    merged = bed_out / f"{strain}.confident.merged.bed"
    conda_sh("truvari", """sort -k1,1 -k2,2n %(raw)s | bedtools merge -i - > %(out)s""" % {
        "raw": conf_raw, "out": merged})

    trimmed = bed_out / f"{strain}.confident.trimmed.bed"
    conda_sh("truvari", """awk -v b=%(buf)s 'BEGIN{OFS="\\t"}{s=$2+b; e=$3-b; if(e>s) print $1,s,e}' \
        %(inp)s > %(out)s""" % {"buf": BOUNDARY_BUFFER, "inp": merged, "out": trimmed})

    windows = bed_out / "windows.bed"
    conda_sh("truvari", """bedtools makewindows -g %(genome)s -w %(win)s > %(out)s""" % {
        "genome": GENOME_FILE, "win": DENSITY_WINDOW, "out": windows})

    variants_bed = bed_out / f"{strain}.variants.bed"
    conda_sh("truvari", """bcftools view -H %(truth)s | awk 'BEGIN{OFS="\\t"}{print $1,$2-1,$2}' > %(out)s""" % {
        "truth": truth_vcf, "out": variants_bed})

    density = bed_out / f"{strain}.density.bed"
    conda_sh("truvari", """bedtools intersect -a %(windows)s -b %(variants)s -c > %(out)s""" % {
        "windows": windows, "variants": variants_bed, "out": density})

    thresh = conda_sh("truvari", """awk '{print $4}' %(density)s | sort -n | \
        awk -v p=%(pct)s '{a[NR]=$1} END{idx=int(NR*p/100); if(idx<1)idx=1; print a[idx]}'""" % {
        "density": density, "pct": DENSITY_PERCENTILE}, capture=True).strip()

    high_density = bed_out / f"{strain}.high_density.bed"
    conda_sh("truvari", """awk -v t=%(t)s '$4 > t {print $1"\\t"$2"\\t"$3}' %(density)s > %(out)s""" % {
        "t": thresh, "density": density, "out": high_density})

    strict = bed_out / f"{strain}.confident.strict.bed"
    conda_sh("truvari", """bedtools subtract -a %(trimmed)s -b %(hd)s > %(out)s""" % {
        "trimmed": trimmed, "hd": high_density, "out": strict})

    vsa_masked = bed_out / f"{strain}.confident.strict.vsa_masked.bed"
    conda_sh("truvari", """bedtools subtract -a %(strict)s -b %(vsa)s > %(out)s""" % {
        "strict": strict, "vsa": VSA_BED, "out": vsa_masked})

    print(f"=== [{strain}] vt on truth + gatk_linear + gatk_surject_loo + vg_hap_loo ===")
    vcfs = {
        "truth": truth_vcf,
        "gatk_linear": OUTROOT / strain / "gatk_linear" / f"{strain}.GATK_linear.shortvars.lt50bp.PASS.GT.decomposed.vcf.gz",
        "gatk_surject_loo": OUTROOT / strain / "gatk_surject_loo" / f"{strain}.GATK_surject_loo.shortvars.lt50bp.PASS.GT.decomposed.vcf.gz",
        "vg_hap_loo": LOO_BASE / f"loo_{strain}" / f"{strain}.loo.pan.shortvars.lt50bp.PASS.hap.dp5.1_1_only.decomposed.vcf.gz",
    }

    all_tracts = bed_out / f"{strain}.all_vt_repeat_tracts.bed"
    all_tracts.write_text("")
    for name, vcf in vcfs.items():
        if not is_file(vcf):
            print(f"  SKIP vt/{name}")
            continue
        vt_vcf = bed_out / f"{strain}.{name}.vt.vcf"
        conda_sh("truvari", """%(vt)s annotate_indels -r %(ref)s %(vcf)s -o %(out)s 2>&1 | tail -1""" % {
            "vt": VT_BIN, "ref": REF, "vcf": vcf, "out": vt_vcf})
        conda_sh("truvari", """bgzip -f %(out)s""" % {"out": vt_vcf})
        conda_sh("truvari", """tabix -f -p vcf %(out)s.gz""" % {"out": vt_vcf})
        conda_sh("truvari", """bcftools view -v indels %(vcf)s.gz | \
            bcftools query -f '%%CHROM\\t%%POS\\t%%INFO/FZ_REPEAT_TRACT\\t%%REF\\t%%ALT\\n' \
                -i "INFO/FZ_MOTIF!=\\"\\" && INFO/FZ_RL>=%(minrl)s" | \
            awk -F'\\t' 'BEGIN{OFS="\\t"}{
                split($3,a,","); ts=a[1]; te=a[2]; vs=$2-1; ve=$2-1+length($4)
                s=(vs<ts-1)?vs:ts-1; e=(ve>te)?ve:te; if(s<0)s=0; print $1,s,e
            }' >> %(tracts)s""" % {
            "vcf": vt_vcf, "minrl": VT_MIN_FZ_RL, "tracts": all_tracts})

    union = bed_out / f"{strain}.vt_repeat_tracts.union.bed"
    conda_sh("truvari", """sort -k1,1 -k2,2n %(tracts)s | bedtools merge -i - > %(out)s""" % {
        "tracts": all_tracts, "out": union})

    final_mask = bed_out / f"{strain}.final_repeat_mask.bed"
    conda_sh("truvari", """cat %(trf)s %(union)s | cut -f1-3 | sort -k1,1 -k2,2n | \
        bedtools merge -i - > %(out)s""" % {"trf": TRF_BED, "union": union, "out": final_mask})

    final_bed = bed_out / f"{strain}.confident.FINAL.bed"
    conda_sh("truvari", """bedtools subtract -a %(vsa_masked)s \
        -b %(mask)s > %(out)s""" % {"vsa_masked": vsa_masked, "mask": final_mask, "out": final_bed})

    for bed in sorted(BEDDIR.glob("*.bed")):
        bedname = bed.stem
        cat_out = bed_out / f"{bedname}.FINAL.bed"
        conda_sh("truvari", """bedtools intersect -a %(bed)s -b %(final)s > %(out)s""" % {
            "bed": bed, "final": final_bed, "out": cat_out})
    # NOTE: SVs deliberately use the single CORE_BED directly (no
    # confidence-tier masking) in run_truvari_all, not this per-category
    # FINAL bed loop.


############################################################
# STEP 7: vcfeval (short vars) -- 4 callers vs paftools truth, all BEDDIR categories
############################################################
def run_vcfeval_all(strain: str):
    bed_out = OUTROOT / strain / "beds"
    out = OUTROOT / strain / "vcfeval"

    truth_vcf = reheader_truth(strain)

    callers = {
        "pan_direct": COMPARISON_DIR / f"{strain}.pan.snps_indels.vcf.gz",   # CONFIRM: needs no further filtering?
        "gatk_linear": OUTROOT / strain / "gatk_linear" / f"{strain}.GATK_linear.shortvars.lt50bp.PASS.GT.decomposed.vcf.gz",
        "gatk_surject_loo": OUTROOT / strain / "gatk_surject_loo" / f"{strain}.GATK_surject_loo.shortvars.lt50bp.PASS.GT.decomposed.vcf.gz",
        "vg_hap_loo": LOO_BASE / f"loo_{strain}" / f"{strain}.loo.pan.shortvars.lt50bp.PASS.hap.dp5.1_1_only.decomposed.vcf.gz",
    }

    for caller, calls in callers.items():
        if not is_file(calls):
            print(f"  SKIP {caller}: not found")
            continue
        for bed in sorted(BEDDIR.glob("*.bed")):
            bedname = bed.stem
            cat_bed = bed_out / f"{bedname}.FINAL.bed"
            bp = sum_bed_bp(cat_bed)
            if bp < MIN_BP:
                continue
            sample_out = out / caller / bedname
            sample_out.parent.mkdir(parents=True, exist_ok=True)
            rmtree_if_exists(sample_out)
            print(f"  -> [{strain}] {caller} | {bedname}")
            conda_sh("truvari", """rtg vcfeval -m annotate --all-records --ref-overlap --squash-ploidy \
                --vcf-score-field=QUAL --bed-regions %(bed)s \
                -b %(truth)s -c %(calls)s -t %(sdf)s --sample %(strain)s \
                -o %(out)s 2>&1 | tail -2""" % {
                "bed": cat_bed, "truth": truth_vcf, "calls": calls, "sdf": REF_SDF,
                "strain": strain, "out": sample_out})


############################################################
# STEP 8: truvari (SVs) -- vs svim truth, core-genome bed only
############################################################
def run_truvari_all(strain: str):
    out = OUTROOT / strain / "truvari"

    sv_truth = COMPARISON_DIR / f"{strain}.svim.svs.vcf.gz"
    alt_sv_truth = OUTROOT / strain / f"{strain}.svim.svs.alt_only.mixed.vcf.gz"
    if not is_file(alt_sv_truth):
        conda_sh("truvari", """bcftools view -i 'GT!="ref"' %(truth)s -Oz -o %(out)s""" % {
            "truth": sv_truth, "out": alt_sv_truth})
        conda_sh("truvari", """bcftools index -f -t %(out)s""" % {"out": alt_sv_truth})

    # Which callers get (re)compared here is driven by the same RUN_SV /
    # RUN_SV_EXTRA / RUN_SV_DYSGU toggles that gate the calling steps above,
    # so re-running with only one toggle True reruns truvari for just that
    # caller set and leaves the already-computed results for the others alone.
    sv_callers = {}
    if RUN_SV:
        sv_callers["pan_direct_sv"] = COMPARISON_DIR / f"{strain}.pan.svs.vcf.gz"
        sv_callers["delly_loo"] = LOO_BASE / f"loo_{strain}" / f"{strain}.loo.delly.SV.ge50bp.PASS.vcf.gz"
        sv_callers["vg_hap_sv_loo"] = LOO_BASE / f"loo_{strain}" / f"{strain}.loo.pan.SV.ge50bp.PASS.hap.vcf.gz"
        sv_callers["delly_linear"] = OUTROOT / strain / "delly_linear" / f"{strain}.delly_linear.SV.ge50bp.PASS.vcf.gz"
    if RUN_SV_EXTRA:
        sv_callers["manta_loo"] = LOO_BASE / f"loo_{strain}" / f"{strain}.loo.manta.SV.ge50bp.PASS.vcf.gz"
        sv_callers["manta_linear"] = OUTROOT / strain / "manta_linear" / f"{strain}.manta_linear.SV.ge50bp.PASS.vcf.gz"
    if RUN_SV_DYSGU:
        sv_callers["dysgu_loo"] = LOO_BASE / f"loo_{strain}" / f"{strain}.loo.dysgu.SV.ge50bp.PASS.vcf.gz"
        sv_callers["dysgu_linear"] = OUTROOT / strain / "dysgu_linear" / f"{strain}.dysgu_linear.SV.ge50bp.PASS.vcf.gz"

    if not sv_callers:
        print("  [skip] run_truvari_all: RUN_SV, RUN_SV_EXTRA, and RUN_SV_DYSGU all false")
        return

    # SVs are evaluated against the core genome only (CONFIRM the exact
    # CORE_BED filename in the CONFIG section), not the full per-category
    # BEDDIR loop that vcfeval (short variants) still uses.
    if not is_file(CORE_BED):
        print(f"  ERROR: core genome bed not found at {CORE_BED} -- check CORE_BED path")
        return
    core_bp = sum_bed_bp(CORE_BED)
    if core_bp < MIN_BP:
        print(f"  ERROR: core genome bed at {CORE_BED} has only {core_bp}bp, below MIN_BP={MIN_BP}")
        return

    # Haploid-consistency filtering for callers that emit diploid-style
    # genotypes (manta, dysgu): a true call in a haploid organism should be
    # effectively fixed in the sample, so a weak 0/1 het at low support is
    # far more likely to be alignment noise than real biology (this is
    # exactly the pattern seen in dysgu's loo FPs -- low-AF hets
    # concentrated in tandem-repeat regions). vg_hap_sv_loo is exempt: it
    # was called with vg call --ploidy 1, so its GTs are already haploid
    # singletons ("0"/"1"), not diploid pairs, and delly is left as-is
    # since it wasn't part of this comparison. smoove was dropped entirely
    # (not performing well enough to be worth comparing).
    #
    # dysgu confirmed to emit FORMAT/AF, so it gets the AF>=0.7 escape
    # hatch (recovers true haploid calls a diploid genotyper miscalled as
    # 0/1 but with strong support). manta does NOT reliably emit an
    # AF-equivalent FORMAT field in this pipeline (it reports PR/SR
    # read-support pairs instead) -- using a tag that isn't declared in
    # the VCF header makes bcftools error out rather than silently skip
    # it, so it gets the stricter GT="1/1"-only filter.
    HAP_FILTER_EXPR = {
        "manta_loo": 'GT="1/1"',
        "manta_linear": 'GT="1/1"',
        "dysgu_loo": 'GT="1/1" || FORMAT/AF>=0.7',
        "dysgu_linear": 'GT="1/1" || FORMAT/AF>=0.7',
    }

    for caller, calls in sv_callers.items():
        if not is_file(calls):
            print(f"  SKIP {caller}")
            continue
        calls_str = str(calls)
        assert calls_str.endswith(".vcf.gz")
        alt_calls = Path(calls_str[: -len(".vcf.gz")] + ".alt_only.mixed.vcf.gz")
        filter_expr = HAP_FILTER_EXPR.get(caller, 'GT!="ref"')
        if not is_file(alt_calls):
            conda_sh("truvari", """bcftools view -i '%(expr)s' %(calls)s -Oz -o %(out)s""" % {
                "expr": filter_expr, "calls": calls, "out": alt_calls})
            conda_sh("truvari", """bcftools index -f -t %(out)s""" % {"out": alt_calls})

        sample_out = out / caller / CORE_BED.stem
        sample_out.parent.mkdir(parents=True, exist_ok=True)
        rmtree_if_exists(sample_out)
        print(f"  -> [{strain}] {caller} | {CORE_BED.stem}")
        conda_sh("truvari", """truvari bench -b %(truth)s -c %(calls)s -o %(out)s -f %(ref)s \
            -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 10000 --includebed %(bed)s""" % {
            "truth": alt_sv_truth, "calls": alt_calls, "out": sample_out, "ref": REF, "bed": CORE_BED})


############################################################
# MAIN
############################################################
def main():
    build_shared_beds()

    for strain in STRAINS:
        print()
        print("##########################################")
        print(f"# {strain}")
        print("##########################################")

        if RUN_MAPPING:
            run_loo(strain)
        else:
            print("[skip] mapping (RUN_MAPPING=false)")

        if RUN_SMALL_VAR:
            run_gatk_linear(strain)
            run_gatk_surject_loo(strain)
            run_vg_hap_loo_filter(strain)
            build_confidence(strain)      # builds FINAL (short-var) confidence tier
            run_vcfeval_all(strain)
        else:
            print("[skip] small-variant filtering + vcfeval (RUN_SMALL_VAR=false)")
            # no confidence-tier bed is needed for SVs (CORE_BED used directly
            # in run_truvari_all), so nothing to backfill here

        if RUN_SV:
            run_delly_loo(strain)
            run_delly_linear(strain)
        else:
            print("[skip] delly (RUN_SV=false)")

        if RUN_SV_EXTRA:
            run_manta_loo(strain)
            run_manta_linear(strain)
        else:
            print("[skip] manta (RUN_SV_EXTRA=false)")

        if RUN_SV_DYSGU:
            run_dysgu_loo(strain)
            run_dysgu_linear(strain)
        else:
            print("[skip] dysgu (RUN_SV_DYSGU=false)")

        if RUN_SV or RUN_SV_EXTRA or RUN_SV_DYSGU:
            run_truvari_all(strain)
        else:
            print("[skip] truvari SV comparison (RUN_SV=false, RUN_SV_EXTRA=false, RUN_SV_DYSGU=false)")

    print(f"DONE. Results in {OUTROOT}/<strain>/vcfeval and .../truvari")


############################################################
# ARGUMENT PARSING
# Every default below is exactly the hardcoded value from the CONFIG
# section above, so running with no arguments reproduces the original
# script exactly. Any path, threshold, or stage toggle can be overridden
# on the command line.
############################################################
def add_bool_flag(ap: argparse.ArgumentParser, name: str, default: bool, help_text: str):
    """--name / --no-name pair with a default, compatible with any Python 3
    version (argparse.BooleanOptionalAction needs Python 3.9+)."""
    dest = name.replace("-", "_")
    group = ap.add_mutually_exclusive_group()
    group.add_argument(f"--{name}", dest=dest, action="store_true", default=None,
                        help=f"{help_text} (default: {default})")
    group.add_argument(f"--no-{name}", dest=dest, action="store_false", default=None,
                        help=argparse.SUPPRESS)
    ap.set_defaults(**{dest: default})


def build_arg_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Pf pangenome LOO validation pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    ap.add_argument(
        "--strains", nargs="+",
        default=["PfDd2", "Pf7G8", "PfCD01", "PfGA01", "PfGB4", "PfGN01", "PfHB3",
                 "PfIT", "PfKE01", "PfKH01", "PfKH02", "PfML01", "PfSN01"],
        help="Strains to process")

    ap.add_argument("--seq-file", type=Path,
                     default=Path("/mnt/storage13/nbillows/pangenome/final_no_mix/pf3k_seq_v2.txt"))
    ap.add_argument("--reads-dir", type=Path,
                     default=Path("/mnt/storage13/nbillows/pangenome/analysis/pan_GT/og_reads/merged_illumina"))
    ap.add_argument("--linear-gvcf-dir", type=Path, default=None,
                     help="Raw linear GATK gvcfs ({strain}.g.vcf.gz); "
                          "defaults to --reads-dir if not given")
    ap.add_argument("--loo-base", type=Path,
                     default=Path("/mnt/storage13/nbillows/pangenome/analysis/pan_GT"))
    ap.add_argument("--comparison-dir", type=Path,
                     default=Path("/mnt/storage13/nbillows/pangenome/analysis/pan_GT/comparison_vcfs"),
                     help="Pre-existing pan_direct + truth files")

    ap.add_argument("--cactus-venv", type=Path,
                     default=Path("/mnt/storage13/nbillows/pangenome/cactus-bin-v2.9.3/venv-cactus-v2.9.3/bin/activate"))
    ap.add_argument("--ref", type=Path,
                     default=Path("/mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/Pfalciparum.genome.fasta"))
    ap.add_argument("--ref-gatk", type=Path,
                     default=Path("/mnt/storage13/nbillows/Pf_09_24/Pfalciparum_09_24_v2/reference/Pf3D7_v3/Pfalciparum.genome.fasta"),
                     help="CONFIRM: same file as --ref, different path used "
                          "historically for GATK steps")
    ap.add_argument("--ref-sdf", type=Path,
                     default=Path("/mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/Pfalciparum.genome.fasta.sdf"))
    ap.add_argument("--vt-bin", type=Path,
                     default=Path("/mnt/storage13/nbillows/pangenome/final_no_mix/review/reviewer2_5/vt/vt"))
    ap.add_argument("--gff", type=Path,
                     default=Path("/mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/Pfalciparum.genome.modified.new.gff3"),
                     help="SET THIS")

    ap.add_argument("--beddir", type=Path,
                     default=Path("/mnt/storage13/nbillows/pangenome/analysis/bed_files"))
    ap.add_argument("--vsa-bed", type=Path,
                     default=Path("/mnt/storage13/nbillows/pangenome/analysis/bed_files/Pfalciparum_variable_repetitive_regions.bed"))
    ap.add_argument("--trf-bed", type=Path,
                     default=Path("/mnt/storage13/nbillows/pangenome/analysis/bed_files/Pfalciparum_TRF_repeats_sensitive.merged.bed"))
    ap.add_argument("--conf-bed-dir", type=Path,
                     default=Path("/mnt/storage13/nbillows/pangenome/analysis/pan_GT/svim_asm_results"))
    ap.add_argument("--core-bed", type=Path, default=Path("/mnt/storage13/nbillows/pangenome/analysis/bed_files/Core_genome_Pf3D7_v3_ext.bed"))

    ap.add_argument("--outroot", type=Path,
                     default=Path("/mnt/storage13/nbillows/pangenome/final_no_mix/review/reviewer2_5/full_pipeline_sep26"))

    ap.add_argument("--vsa-pad", type=int, default=500)
    ap.add_argument("--vt-min-fz-rl", type=int, default=4)
    ap.add_argument("--boundary-buffer", type=int, default=1000)
    ap.add_argument("--density-window", type=int, default=200)
    ap.add_argument("--density-percentile", type=int, default=95)
    ap.add_argument("--min-bp", type=int, default=1000)

    add_bool_flag(ap, "run-mapping", False,
                  "run_loo (cactus build, giraffe, surject, pack, vg call, gatk-on-surject-gvcf)")
    add_bool_flag(ap, "run-small-var", False,
                  "gatk_linear/gatk_surject_loo/vg_hap_loo filtering + vcfeval")
    add_bool_flag(ap, "run-sv", False,
                  "delly_loo + delly_linear (+ their truvari comparison)")
    add_bool_flag(ap, "run-sv-extra", False,
                  "manta (loo + linear)")
    add_bool_flag(ap, "run-sv-dysgu", True,
                  "dysgu (loo + linear)")

    return ap


def apply_args(args: argparse.Namespace):
    """Overwrite the module-level config globals with parsed CLI values.
    Called once before main() so every function below (which reads these
    as module-level globals, unchanged from the original bash-variable
    style) sees the final, possibly-overridden values."""
    global STRAINS, SEQ_FILE, READS_DIR, LINEAR_GVCF_DIR, LOO_BASE, COMPARISON_DIR
    global CACTUS_VENV, REF, REF_GATK, REF_SDF, VT_BIN, GFF
    global BEDDIR, VSA_BED, TRF_BED, CONF_BED_DIR, CORE_BED, OUTROOT
    global VSA_PAD, VT_MIN_FZ_RL, BOUNDARY_BUFFER, DENSITY_WINDOW, DENSITY_PERCENTILE, MIN_BP
    global RUN_MAPPING, RUN_SMALL_VAR, RUN_SV, RUN_SV_EXTRA, RUN_SV_DYSGU

    STRAINS = args.strains
    SEQ_FILE = args.seq_file
    READS_DIR = args.reads_dir
    # preserves the original "defaults to reads-dir" aliasing even when
    # --reads-dir itself is overridden and --linear-gvcf-dir is not
    LINEAR_GVCF_DIR = args.linear_gvcf_dir if args.linear_gvcf_dir is not None else READS_DIR
    LOO_BASE = args.loo_base
    COMPARISON_DIR = args.comparison_dir

    CACTUS_VENV = args.cactus_venv
    REF = args.ref
    REF_GATK = args.ref_gatk
    REF_SDF = args.ref_sdf
    VT_BIN = args.vt_bin
    GFF = args.gff

    BEDDIR = args.beddir
    VSA_BED = args.vsa_bed
    TRF_BED = args.trf_bed
    CONF_BED_DIR = args.conf_bed_dir
    # preserves the "CORE_BED lives under BEDDIR" relationship even when
    # --beddir itself is overridden and --core-bed is not
    CORE_BED = args.core_bed if args.core_bed is not None else BEDDIR / "core_genome.bed"

    OUTROOT = args.outroot
    OUTROOT.mkdir(parents=True, exist_ok=True)

    VSA_PAD = args.vsa_pad
    VT_MIN_FZ_RL = args.vt_min_fz_rl
    BOUNDARY_BUFFER = args.boundary_buffer
    DENSITY_WINDOW = args.density_window
    DENSITY_PERCENTILE = args.density_percentile
    MIN_BP = args.min_bp

    RUN_MAPPING = args.run_mapping
    RUN_SMALL_VAR = args.run_small_var
    RUN_SV = args.run_sv
    RUN_SV_EXTRA = args.run_sv_extra
    RUN_SV_DYSGU = args.run_sv_dysgu


if __name__ == "__main__":
    _parser = build_arg_parser()
    _args = _parser.parse_args()
    apply_args(_args)
    main()
