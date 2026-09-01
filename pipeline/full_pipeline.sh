#!/usr/bin/env bash
set -euo pipefail

############################################################
# CONFIG
############################################################
STRAINS=(PfDd2 Pf7G8 PfCD01 PfGA01 PfGB4 PfGN01 PfHB3 PfIT PfKE01 PfKH01 PfKH02 PfML01 PfSN01 PfTG01)

SEQ_FILE="/mnt/storage13/nbillows/pangenome/final_no_mix/pf3k_seq_v2.txt"
READS_DIR="/mnt/storage13/nbillows/pangenome/analysis/pan_GT/og_reads/merged_illumina"
LINEAR_GVCF_DIR="${READS_DIR}"   # raw linear GATK gvcfs live here per strain: {strain}.g.vcf.gz
LOO_BASE="/mnt/storage13/nbillows/pangenome/analysis/pan_GT"
COMPARISON_DIR="/mnt/storage13/nbillows/pangenome/analysis/pan_GT/comparison_vcfs"   # pre-existing pan_direct + truth files

CACTUS_VENV="/mnt/storage13/nbillows/pangenome/cactus-bin-v2.9.3/venv-cactus-v2.9.3/bin/activate"
REF="/mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/Pfalciparum.genome.fasta"
REF_GATK="/mnt/storage13/nbillows/Pf_09_24/Pfalciparum_09_24_v2/reference/Pf3D7_v3/Pfalciparum.genome.fasta"  # CONFIRM: same file as REF, different path used historically for GATK steps
REF_SDF="/mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/Pfalciparum.genome.fasta.sdf"
VT_BIN="/mnt/storage13/nbillows/pangenome/final_no_mix/review/reviewer2_5/vt/vt"
GFF="/path/to/Pfalciparum_genome_modified_new.gff3"   # SET THIS

BEDDIR="/mnt/storage13/nbillows/pangenome/analysis/bed_files"
VSA_BED="/mnt/storage13/nbillows/pangenome/final_no_mix/review/reviewer2_5/Pfalciparum_variable_repetitive_regions.bed"
TRF_BED="/mnt/storage13/nbillows/pangenome/final_no_mix/review/reviewer2_5/Pfalciparum_TRF_repeats_sensitive.merged.bed"
CONF_BED_DIR="/mnt/storage13/nbillows/pangenome/analysis/pan_GT/svim_asm_results"

OUTROOT="/mnt/storage13/nbillows/pangenome/final_no_mix/review/reviewer2_5/full_pipeline"
mkdir -p "${OUTROOT}"

VSA_PAD=500
VT_MIN_FZ_RL=4
BOUNDARY_BUFFER=1000
DENSITY_WINDOW=200
DENSITY_PERCENTILE=95
MIN_BP=1000

############################################################
# STAGE TOGGLES
# Mapping (LOO graph build + giraffe) and small-variant validation
# are already complete for validated strains and should not be
# rerun by default. Set to true to re-enable for new strains.
############################################################
RUN_MAPPING=false        # run_loo (cactus build, giraffe, surject, pack, vg call, gatk-on-surject-gvcf)
RUN_SMALL_VAR=false      # gatk_linear/gatk_surject_loo/vg_hap_loo filtering + vcfeval
RUN_SV=true              # delly_loo + truvari SV comparison

############################################################
# ENV HELPERS
############################################################
conda_activate() {
    source "$(conda info --base)/etc/profile.d/conda.sh"
    set +u
    conda activate "$1"
    set -u
}

############################################################
# STEP 0: shared, genome-wide beds (build once)
############################################################
build_shared_beds() {
    if [[ ! -f "${VSA_BED}" ]]; then
        echo "=== building VSA bed ==="
        local tmp; tmp="$(dirname "${VSA_BED}")/vsa_tmp"; mkdir -p "${tmp}"
        awk -F'\t' '$3=="gene"{print}' "${GFF}" | \
            grep -P 'Name=(VAR|VAR-like|VAR2CSA|RIF|RIFA|MC-2TM|SURF)' | \
            awk -F'\t' 'BEGIN{OFS="\t"}{print $1,$4-1,$5}' > "${tmp}/vsa_genes.bed"
        awk -F'\t' '$3=="repeat_region"{print $1"\t"($4-1)"\t"$5}' "${GFF}" > "${tmp}/repeat_regions.bed"
        awk -F'\t' '$3=="centromere"{print $1"\t"($4-1)"\t"$5}' "${GFF}" > "${tmp}/centromeres.bed"
        cat "${tmp}"/*.bed | sort -k1,1 -k2,2n > "${tmp}/combined_raw.bed"
        awk -v p="${VSA_PAD}" 'BEGIN{OFS="\t"}{s=$2-p; if(s<0)s=0; print $1,s,$3+p}' "${tmp}/combined_raw.bed" | \
            sort -k1,1 -k2,2n | bedtools merge -i - > "${VSA_BED}"
    fi
    [[ -f "${TRF_BED}" ]] || { echo "ERROR: TRF bed not found"; exit 1; }

    GENOME_FILE="${OUTROOT}/genome.txt"
    if [[ ! -f "${GENOME_FILE}" ]]; then
        samtools faidx "${REF}"
        cut -f1,2 "${REF}.fai" > "${GENOME_FILE}"
    fi
}

############################################################
# STEP 1: LOO graph build + read mapping + vg call + GATK-on-surject
############################################################
run_loo() {
    local STRAIN="$1"
    local WORK="${LOO_BASE}/loo_${STRAIN}"
    mkdir -p "${WORK}"; cd "${WORK}"

    echo "=== [${STRAIN}] LOO seq file ==="
    grep -v "^${STRAIN}\b" "${SEQ_FILE}" > "pf3k_seq_no_${STRAIN}.txt"

    local GRAPH_NAME="PfPan_no_${STRAIN}"
    local OUT_DIR="${WORK}/${GRAPH_NAME}_pan"
    local GBZ="${OUT_DIR}/${GRAPH_NAME}.gbz"
    local SNARLS="${OUT_DIR}/${GRAPH_NAME}.snarls"

    echo "=== [${STRAIN}] cactus-pangenome (build if missing) ==="
    if [[ ! -f "${GBZ}" ]]; then
        bash -c "source ${CACTUS_VENV} && cactus-pangenome ${WORK}/${GRAPH_NAME}_js ${WORK}/pf3k_seq_no_${STRAIN}.txt \
            --outDir ${OUT_DIR} --outName ${GRAPH_NAME} --reference Pf3D7 \
            --filter 2 --haplo --giraffe clip filter \
            --gbz clip filter full --gfa clip full \
            --vcf --vcfReference Pf3D7 \
            --logFile ${WORK}/${GRAPH_NAME}.log \
            --workDir ${WORK} --consCores 8 --mgMemory 128Gi \
            --restart && deactivate"
    fi
    [[ -f "${GBZ}" ]] || { echo "ERROR: GBZ build failed for ${STRAIN}"; return 1; }

    conda_activate cactus

    echo "=== [${STRAIN}] vg giraffe + surject ==="
    local FQ1="${READS_DIR}/${STRAIN}_1.trimmed.fastq.gz"
    local FQ2="${READS_DIR}/${STRAIN}_2.trimmed.fastq.gz"
    local GAF="${WORK}/${STRAIN}.loo.gaf.gz"
    local SORT_BAM="${WORK}/${STRAIN}.loo.sort.bam"

    if [[ ! -s "${GAF}" ]]; then
        vg giraffe -p -t 16 -Z "${GBZ}" -f "${FQ1}" -f "${FQ2}" -o gaf | bgzip > "${GAF}"
    fi
    vg paths -x "${GBZ}" -L | grep "^Pf3D7" > "${WORK}/Pf3D7.paths.txt"

    if [[ ! -s "${SORT_BAM}" ]]; then
        vg surject -x "${GBZ}" -G "${GAF}" --interleaved \
            -F "${WORK}/Pf3D7.paths.txt" -b \
            -N "${STRAIN}" -R "ID:1 LB:lib1 SM:${STRAIN} PL:illumina PU:unit1" \
            | samtools reheader -c 'sed s/Pf3D7#0#//g' - > "${WORK}/${STRAIN}.loo.bam"
        samtools sort "${WORK}/${STRAIN}.loo.bam" -O BAM -o "${SORT_BAM}" --threads 8
        samtools index "${SORT_BAM}" --threads 8
    fi

    echo "=== [${STRAIN}] vg pack + vg call (haploid) ==="
    local PACK="${WORK}/${STRAIN}.loo.pack"
    local HAP_VCF="${WORK}/${STRAIN}.loo.SV.call_haploid_.vcf.gz"
    [[ -s "${PACK}" ]] || vg pack -x "${GBZ}" -Q5 -a "${GAF}" -o "${PACK}"
    if [[ ! -s "${HAP_VCF}" ]]; then
        vg call "${GBZ}" -r "${SNARLS}" -k "${PACK}" -t 16 --ploidy 1 \
            -s "${STRAIN}" -S Pf3D7 -az | bgzip > "${HAP_VCF}"
    fi

    echo "=== [${STRAIN}] fix contig naming on raw vg call output ==="
    fix_contigs "${HAP_VCF}" "loo_${STRAIN}"

    echo "=== [${STRAIN}] GATK on LOO surjected BAM ==="
    conda_activate fastq2matrix
    local sample="${STRAIN}.loo.gatk_surject"
    if [[ ! -s "${WORK}/${sample}.g.vcf.gz" ]]; then
        gatk HaplotypeCaller -I "${SORT_BAM}" -R "${REF_GATK}" -O "${WORK}/${sample}.g.vcf.gz" -ERC GVCF
    fi
}

############################################################
# fix_contigs: strip Pf3D7#0# prefix from vg call output (header + body)
############################################################
fix_contigs() {
    local vcf="$1" tag="$2"
    conda_activate cactus
    bcftools index -f -t "${vcf}"
    bcftools view -h "${vcf}" | grep "^##contig" | \
        sed -E 's/##contig=<ID=([^,]+),.*/\1/' | \
        awk '{new=$0; gsub("Pf3D7#0#","",new); print $0"\t"new}' > "rename_map_${tag}.txt"
    bcftools annotate --rename-chrs "rename_map_${tag}.txt" "${vcf}" -Oz -o "body_renamed_${tag}.vcf.gz"
    bcftools view -h "body_renamed_${tag}.vcf.gz" | sed 's/Pf3D7#0#//g' > "header_fixed_${tag}.txt"
    bcftools reheader -h "header_fixed_${tag}.txt" "body_renamed_${tag}.vcf.gz" -o "${vcf%.vcf.gz}.renamed.vcf.gz"
    bcftools index -f -t "${vcf%.vcf.gz}.renamed.vcf.gz"
    mv "${vcf%.vcf.gz}.renamed.vcf.gz" "${vcf}"
    mv "${vcf%.vcf.gz}.renamed.vcf.gz.tbi" "${vcf}.tbi"
    rm -f "body_renamed_${tag}.vcf.gz" "header_fixed_${tag}.txt" "rename_map_${tag}.txt"
}

############################################################
# GATK filter chain, reused for linear and for LOO-surject GVCFs
############################################################
filter_gatk() {
    local gvcf="$1" sample="$2" workdir="$3" final_out="$4"
    conda_activate fastq2matrix
    cd "${workdir}"

    gatk GenotypeGVCFs -R "${REF_GATK}" -V "${gvcf}" -O "${sample}.raw.vcf.gz"
    tabix -f -p vcf "${sample}.raw.vcf.gz"

    gatk SelectVariants -R "${REF_GATK}" -V "${sample}.raw.vcf.gz" --select-type-to-include SNP -O "${sample}.snps.vcf.gz"
    tabix -f -p vcf "${sample}.snps.vcf.gz"
    gatk SelectVariants -R "${REF_GATK}" -V "${sample}.raw.vcf.gz" --select-type-to-include INDEL -O "${sample}.indels.vcf.gz"
    tabix -f -p vcf "${sample}.indels.vcf.gz"

    gatk VariantFiltration -R "${REF_GATK}" -V "${sample}.snps.vcf.gz" \
        --filter-expression "QD < 2.0"             --filter-name "QD2" \
        --filter-expression "QUAL < 30.0"           --filter-name "QUAL30" \
        --filter-expression "SOR > 4.0"             --filter-name "SOR4" \
        --filter-expression "FS > 60.0"             --filter-name "FS60" \
        --filter-expression "MQ < 40.0"             --filter-name "MQ40" \
        --filter-expression "MQRankSum < -15.0"     --filter-name "MQRankSum-15" \
        --filter-expression "ReadPosRankSum < -5.0" --filter-name "ReadPosRankSum-5" \
        -O "${sample}.snps.filtered.vcf.gz"
    tabix -f -p vcf "${sample}.snps.filtered.vcf.gz"

    gatk VariantFiltration -R "${REF_GATK}" -V "${sample}.indels.vcf.gz" \
        --filter-expression "QD < 2.0"               --filter-name "QD2" \
        --filter-expression "FS > 200.0"             --filter-name "FS200" \
        --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
        -O "${sample}.indels.filtered.vcf.gz"
    tabix -f -p vcf "${sample}.indels.filtered.vcf.gz"

    bcftools concat -a "${sample}.snps.filtered.vcf.gz" "${sample}.indels.filtered.vcf.gz" -Oz -o "${sample}.filtered.vcf.gz"
    tabix -f -p vcf "${sample}.filtered.vcf.gz"
    gatk SelectVariants -R "${REF_GATK}" -V "${sample}.filtered.vcf.gz" --exclude-filtered -O "${sample}.PASS.vcf.gz"
    tabix -f -p vcf "${sample}.PASS.vcf.gz"
    bcftools norm -f "${REF_GATK}" -m -both "${sample}.PASS.vcf.gz" -Oz -o "${sample}.norm.vcf.gz"
    tabix -f -p vcf "${sample}.norm.vcf.gz"

    bcftools view -i 'TYPE="snp" || (TYPE="indel" && abs(strlen(ALT)-strlen(REF))<50)' \
        "${sample}.norm.vcf.gz" -Oz -o "${sample}.shortvars.lt50bp.PASS.vcf.gz"
    tabix -f -p vcf "${sample}.shortvars.lt50bp.PASS.vcf.gz"

    bcftools view "${sample}.shortvars.lt50bp.PASS.vcf.gz" | setGT.py --fraction 0.7 | \
        bcftools view -O z -c 1 -o "${sample}.shortvars.lt50bp.PASS.GT.vcf.gz"
    tabix -f -p vcf "${sample}.shortvars.lt50bp.PASS.GT.vcf.gz"

    conda_activate truvari
    rm -f "${final_out}" "${final_out}.tbi"
    rtg vcfdecompose --break-mnps --break-indels -i "${sample}.shortvars.lt50bp.PASS.GT.vcf.gz" -o "${final_out}"
    tabix -f -p vcf "${final_out}"
}

############################################################
# STEP 2: linear GATK (raw gvcfs -> filtered), reused across strains
############################################################
run_gatk_linear() {
    local STRAIN="$1"
    local WORK="${OUTROOT}/${STRAIN}/gatk_linear"; mkdir -p "${WORK}"
    local GVCF="${LINEAR_GVCF_DIR}/${STRAIN}.g.vcf.gz"
    [[ -f "${GVCF}" ]] || { echo "WARNING: no linear gvcf for ${STRAIN}"; return 1; }
    filter_gatk "${GVCF}" "${STRAIN}.gatk_linear" "${WORK}" \
        "${WORK}/${STRAIN}.GATK_linear.shortvars.lt50bp.PASS.GT.decomposed.vcf.gz"
}

############################################################
# STEP 3: GATK-on-LOO-surject filtering (uses gvcf from run_loo)
############################################################
run_gatk_surject_loo() {
    local STRAIN="$1"
    local LOO_WORK="${LOO_BASE}/loo_${STRAIN}"
    local WORK="${OUTROOT}/${STRAIN}/gatk_surject_loo"; mkdir -p "${WORK}"
    local GVCF="${LOO_WORK}/${STRAIN}.loo.gatk_surject.g.vcf.gz"
    [[ -f "${GVCF}" ]] || { echo "WARNING: no LOO surject gvcf for ${STRAIN}, run_loo first"; return 1; }
    filter_gatk "${GVCF}" "${STRAIN}.gatk_surject_loo" "${WORK}" \
        "${WORK}/${STRAIN}.GATK_surject_loo.shortvars.lt50bp.PASS.GT.decomposed.vcf.gz"
}

############################################################
# STEP 4: vg_hap_loo filtering (mirrors prepare_variants_3.sh)
############################################################
run_vg_hap_loo_filter() {
    local STRAIN="$1"
    local WORK="${LOO_BASE}/loo_${STRAIN}"; cd "${WORK}"
    conda_activate cactus
    local vcf="${STRAIN}.loo.SV.call_haploid_.vcf.gz"
    local sample="${STRAIN}.loo"

    bcftools norm -f "${REF}" -m -both "${vcf}" -Oz -o "${sample}.pan.norm.unsorted.hap.vcf.gz"
    bcftools sort "${sample}.pan.norm.unsorted.hap.vcf.gz" -Oz -o "${sample}.pan.norm.hap.vcf.gz"
    tabix -f -p vcf "${sample}.pan.norm.hap.vcf.gz"

    bcftools view -f PASS -v snps "${sample}.pan.norm.hap.vcf.gz" -Oz -o "${sample}.pan.snps.PASS.hap.vcf.gz"
    tabix -f -p vcf "${sample}.pan.snps.PASS.hap.vcf.gz"
    bcftools view -f PASS -v indels -i 'abs(strlen(ALT)-strlen(REF))<50' \
        "${sample}.pan.norm.hap.vcf.gz" -Oz -o "${sample}.pan.indels.lt50bp.PASS.hap.vcf.gz"
    tabix -f -p vcf "${sample}.pan.indels.lt50bp.PASS.hap.vcf.gz"
    bcftools concat -a "${sample}.pan.snps.PASS.hap.vcf.gz" "${sample}.pan.indels.lt50bp.PASS.hap.vcf.gz" \
        -Oz -o "${sample}.pan.shortvars.lt50bp.PASS.hap.vcf.gz"
    tabix -f -p vcf "${sample}.pan.shortvars.lt50bp.PASS.hap.vcf.gz"
    bcftools view -f PASS -i 'abs(strlen(ALT)-strlen(REF))>=50' \
        "${sample}.pan.norm.hap.vcf.gz" -Oz -o "${sample}.pan.SV.ge50bp.PASS.hap.vcf.gz"
    tabix -f -p vcf "${sample}.pan.SV.ge50bp.PASS.hap.vcf.gz"

    bcftools +setGT "${sample}.pan.shortvars.lt50bp.PASS.hap.vcf.gz" -- -t q -n 'c:0/0' -i 'GT="0"' \
        | bcftools +setGT -- -t q -n 'c:1/1' -i 'GT="1"' | bgzip -c > tmp.vcf.gz
    mv tmp.vcf.gz "${sample}.pan.shortvars.lt50bp.PASS.hap.vcf.gz"
    tabix -f -p vcf "${sample}.pan.shortvars.lt50bp.PASS.hap.vcf.gz"

    bcftools view -i 'FORMAT/DP>5' "${sample}.pan.shortvars.lt50bp.PASS.hap.vcf.gz" -Oz -o "${sample}.pan.shortvars.lt50bp.PASS.hap.dp5.vcf.gz"
    tabix -f -p vcf "${sample}.pan.shortvars.lt50bp.PASS.hap.dp5.vcf.gz"
    bcftools view -i 'GT="1/1" || GT="1"' "${sample}.pan.shortvars.lt50bp.PASS.hap.dp5.vcf.gz" \
        -Oz -o "${sample}.pan.shortvars.lt50bp.PASS.hap.dp5.1_1_only.vcf.gz"
    tabix -f -p vcf "${sample}.pan.shortvars.lt50bp.PASS.hap.dp5.1_1_only.vcf.gz"

    conda_activate truvari
    local FINAL_LOO="${sample}.pan.shortvars.lt50bp.PASS.hap.dp5.1_1_only.decomposed.vcf.gz"
    rm -f "${FINAL_LOO}" "${FINAL_LOO}.tbi"
    rtg vcfdecompose --break-mnps --break-indels \
        -i "${sample}.pan.shortvars.lt50bp.PASS.hap.dp5.1_1_only.vcf.gz" \
        -o "${FINAL_LOO}"
    tabix -f -p vcf "${FINAL_LOO}"
}

############################################################
# STEP 5: delly on LOO surjected BAM (SV caller, no linear-delly per your instruction)
############################################################
run_delly_loo() {
    local STRAIN="$1"
    local WORK="${LOO_BASE}/loo_${STRAIN}"; cd "${WORK}"
    conda_activate delly   # CONFIRM: env name
    delly call -g "${REF}" "${STRAIN}.loo.sort.bam" -o "${STRAIN}.loo.delly.bcf"
    conda_activate fastq2matrix
    bcftools view "${STRAIN}.loo.delly.bcf" -f PASS \
        -i 'INFO/SVTYPE!="" && abs(INFO/SVLEN)>=50' \
        -Oz -o "${STRAIN}.loo.delly.SV.ge50bp.PASS.vcf.gz"
    tabix -f -p vcf "${STRAIN}.loo.delly.SV.ge50bp.PASS.vcf.gz"
}

############################################################
# STEP 5b: delly on the existing linear (BWA/BQSR) BAM -- true linear baseline
############################################################
run_delly_linear() {
    local STRAIN="$1"
    local WORK="${OUTROOT}/${STRAIN}/delly_linear"; mkdir -p "${WORK}"; cd "${WORK}"
    local BAM="${READS_DIR}/${STRAIN}.bqsr.cram"
    [[ -f "${BAM}" ]] || { echo "WARNING: no linear bqsr cram for ${STRAIN} at ${BAM}"; return 1; }

    conda_activate fastq2matrix
    [[ -f "${BAM}.crai" ]] || samtools index "${BAM}"

    conda_activate delly   # CONFIRM: env name
    if [[ ! -s "${STRAIN}.delly_linear.bcf" ]]; then
        delly call -g "${REF}" "${BAM}" -o "${STRAIN}.delly_linear.bcf"
    fi

    conda_activate fastq2matrix
    bcftools view "${STRAIN}.delly_linear.bcf" -f PASS \
        -i 'INFO/SVTYPE!="" && abs(INFO/SVLEN)>=50' \
        -Oz -o "${STRAIN}.delly_linear.SV.ge50bp.PASS.vcf.gz"
    tabix -f -p vcf "${STRAIN}.delly_linear.SV.ge50bp.PASS.vcf.gz"
}

############################################################
# HELPER: reheader raw paftools truth VCF to use the strain name as sample
############################################################
reheader_truth() {
    local STRAIN="$1"
    local RAW="${COMPARISON_DIR}/${STRAIN}.paftools.snps_indels.vcf.gz"
    local RENAMED="${OUTROOT}/${STRAIN}/${STRAIN}.paftools.snps_indels.reheader.vcf.gz"
    mkdir -p "${OUTROOT}/${STRAIN}"
    conda_activate truvari

    if [[ ! -f "${RENAMED}" ]]; then
        local current_name
        current_name=$(bcftools view -h "${RAW}" | tail -1 | awk -F'\t' '{print $NF}')
        echo -e "${current_name}\t${STRAIN}" > "${OUTROOT}/${STRAIN}/rename.txt"
        bcftools reheader -s "${OUTROOT}/${STRAIN}/rename.txt" -o "${RENAMED}" "${RAW}"
        tabix -f -p vcf "${RENAMED}"
        rm -f "${OUTROOT}/${STRAIN}/rename.txt"
    fi
    TRUTH_VCF_RENAMED="${RENAMED}"
}

############################################################
# STEP 6: confidence tiers (base + VSA + vt-union [short var] / TRF [SV])
############################################################
build_confidence() {
    local STRAIN="$1"
    local BED_OUT="${OUTROOT}/${STRAIN}/beds"; mkdir -p "${BED_OUT}"
    conda_activate truvari

    reheader_truth "${STRAIN}"
    local TRUTH_VCF="${TRUTH_VCF_RENAMED}"
    local CONF_RAW="${CONF_BED_DIR}/${STRAIN}.April2018/${STRAIN}.April2018.confident.bed"
    [[ -f "${CONF_RAW}" ]] || { echo "WARNING: no raw confidence bed for ${STRAIN}"; return 1; }

    sort -k1,1 -k2,2n "${CONF_RAW}" | bedtools merge -i - > "${BED_OUT}/${STRAIN}.confident.merged.bed"
    awk -v b="${BOUNDARY_BUFFER}" 'BEGIN{OFS="\t"}{s=$2+b; e=$3-b; if(e>s) print $1,s,e}' \
        "${BED_OUT}/${STRAIN}.confident.merged.bed" > "${BED_OUT}/${STRAIN}.confident.trimmed.bed"

    bedtools makewindows -g "${GENOME_FILE}" -w "${DENSITY_WINDOW}" > "${BED_OUT}/windows.bed"
    bcftools view -H "${TRUTH_VCF}" | awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2}' > "${BED_OUT}/${STRAIN}.variants.bed"
    bedtools intersect -a "${BED_OUT}/windows.bed" -b "${BED_OUT}/${STRAIN}.variants.bed" -c > "${BED_OUT}/${STRAIN}.density.bed"
    local THRESH
    THRESH=$(awk '{print $4}' "${BED_OUT}/${STRAIN}.density.bed" | sort -n | \
        awk -v p="${DENSITY_PERCENTILE}" '{a[NR]=$1} END{idx=int(NR*p/100); if(idx<1)idx=1; print a[idx]}')
    awk -v t="${THRESH}" '$4 > t {print $1"\t"$2"\t"$3}' "${BED_OUT}/${STRAIN}.density.bed" > "${BED_OUT}/${STRAIN}.high_density.bed"
    bedtools subtract -a "${BED_OUT}/${STRAIN}.confident.trimmed.bed" -b "${BED_OUT}/${STRAIN}.high_density.bed" \
        > "${BED_OUT}/${STRAIN}.confident.strict.bed"
    bedtools subtract -a "${BED_OUT}/${STRAIN}.confident.strict.bed" -b "${VSA_BED}" \
        > "${BED_OUT}/${STRAIN}.confident.strict.vsa_masked.bed"

    echo "=== [${STRAIN}] vt on truth + gatk_linear + gatk_surject_loo + vg_hap_loo ==="
    declare -A VCFS
    VCFS["truth"]="${TRUTH_VCF}"
    VCFS["gatk_linear"]="${OUTROOT}/${STRAIN}/gatk_linear/${STRAIN}.GATK_linear.shortvars.lt50bp.PASS.GT.decomposed.vcf.gz"
    VCFS["gatk_surject_loo"]="${OUTROOT}/${STRAIN}/gatk_surject_loo/${STRAIN}.GATK_surject_loo.shortvars.lt50bp.PASS.GT.decomposed.vcf.gz"
    VCFS["vg_hap_loo"]="${LOO_BASE}/loo_${STRAIN}/${STRAIN}.loo.pan.shortvars.lt50bp.PASS.hap.dp5.1_1_only.decomposed.vcf.gz"

    > "${BED_OUT}/${STRAIN}.all_vt_repeat_tracts.bed"
    for NAME in "${!VCFS[@]}"; do
        local VCF="${VCFS[$NAME]}"
        [[ -f "${VCF}" ]] || { echo "  SKIP vt/${NAME}"; continue; }
        "${VT_BIN}" annotate_indels -r "${REF}" "${VCF}" -o "${BED_OUT}/${STRAIN}.${NAME}.vt.vcf" 2>&1 | tail -1
        bgzip -f "${BED_OUT}/${STRAIN}.${NAME}.vt.vcf"; tabix -f -p vcf "${BED_OUT}/${STRAIN}.${NAME}.vt.vcf.gz"
        bcftools view -v indels "${BED_OUT}/${STRAIN}.${NAME}.vt.vcf.gz" | \
            bcftools query -f '%CHROM\t%POS\t%INFO/FZ_REPEAT_TRACT\t%REF\t%ALT\n' \
                -i "INFO/FZ_MOTIF!=\"\" && INFO/FZ_RL>=${VT_MIN_FZ_RL}" | \
            awk -F'\t' 'BEGIN{OFS="\t"}{
                split($3,a,","); ts=a[1]; te=a[2]; vs=$2-1; ve=$2-1+length($4)
                s=(vs<ts-1)?vs:ts-1; e=(ve>te)?ve:te; if(s<0)s=0; print $1,s,e
            }' >> "${BED_OUT}/${STRAIN}.all_vt_repeat_tracts.bed"
    done
    unset VCFS
    sort -k1,1 -k2,2n "${BED_OUT}/${STRAIN}.all_vt_repeat_tracts.bed" | bedtools merge -i - \
        > "${BED_OUT}/${STRAIN}.vt_repeat_tracts.union.bed"

    cat "${TRF_BED}" "${BED_OUT}/${STRAIN}.vt_repeat_tracts.union.bed" | cut -f1-3 | sort -k1,1 -k2,2n | \
        bedtools merge -i - > "${BED_OUT}/${STRAIN}.final_repeat_mask.bed"
    bedtools subtract -a "${BED_OUT}/${STRAIN}.confident.strict.vsa_masked.bed" \
        -b "${BED_OUT}/${STRAIN}.final_repeat_mask.bed" > "${BED_OUT}/${STRAIN}.confident.FINAL.bed"

    for BED in "${BEDDIR}"/*.bed; do
        local BEDNAME; BEDNAME=$(basename "${BED}" .bed)
        bedtools intersect -a "${BED}" -b "${BED_OUT}/${STRAIN}.confident.FINAL.bed" > "${BED_OUT}/${BEDNAME}.FINAL.bed"
    done
    # NOTE: SVs deliberately use the plain category BEDs directly (no confidence-tier
    # masking) in run_truvari_all, matching the original (non-LOO) truvari pipeline's
    # stringency. No SV_TIER bed is built here.
}

############################################################
# STEP 7: vcfeval (short vars) -- 4 callers vs paftools truth
############################################################
run_vcfeval_all() {
    local STRAIN="$1"
    local BED_OUT="${OUTROOT}/${STRAIN}/beds"
    local OUT="${OUTROOT}/${STRAIN}/vcfeval"
    conda_activate truvari

    reheader_truth "${STRAIN}"
    local TRUTH_VCF="${TRUTH_VCF_RENAMED}"

    declare -A CALLERS
    CALLERS["pan_direct"]="${COMPARISON_DIR}/${STRAIN}.pan.snps_indels.vcf.gz"   # CONFIRM: needs no further filtering?
    CALLERS["gatk_linear"]="${OUTROOT}/${STRAIN}/gatk_linear/${STRAIN}.GATK_linear.shortvars.lt50bp.PASS.GT.decomposed.vcf.gz"
    CALLERS["gatk_surject_loo"]="${OUTROOT}/${STRAIN}/gatk_surject_loo/${STRAIN}.GATK_surject_loo.shortvars.lt50bp.PASS.GT.decomposed.vcf.gz"
    CALLERS["vg_hap_loo"]="${LOO_BASE}/loo_${STRAIN}/${STRAIN}.loo.pan.shortvars.lt50bp.PASS.hap.dp5.1_1_only.decomposed.vcf.gz"

    for CALLER in "${!CALLERS[@]}"; do
        local CALLS="${CALLERS[$CALLER]}"
        [[ -f "${CALLS}" ]] || { echo "  SKIP ${CALLER}: not found"; continue; }
        for BED in "${BEDDIR}"/*.bed; do
            local BEDNAME; BEDNAME=$(basename "${BED}" .bed)
            local CAT_BED="${BED_OUT}/${BEDNAME}.FINAL.bed"
            local BP; BP=$(awk '{s+=$3-$2}END{print s+0}' "${CAT_BED}")
            [[ "${BP}" -lt "${MIN_BP}" ]] && continue
            local SAMPLE_OUT="${OUT}/${CALLER}/${BEDNAME}"
            mkdir -p "$(dirname "${SAMPLE_OUT}")"; rm -rf "${SAMPLE_OUT}"
            echo "  -> [${STRAIN}] ${CALLER} | ${BEDNAME}"
            rtg vcfeval -m annotate --all-records --ref-overlap --squash-ploidy \
                --vcf-score-field=QUAL --bed-regions "${CAT_BED}" \
                -b "${TRUTH_VCF}" -c "${CALLS}" -t "${REF_SDF}" --sample "${STRAIN}" \
                -o "${SAMPLE_OUT}" 2>&1 | tail -2
        done
    done
    unset CALLERS
}

############################################################
# STEP 8: truvari (SVs) -- pan_direct, delly_loo, vg_hap_sv_loo vs svim truth
############################################################
run_truvari_all() {
    local STRAIN="$1"
    local OUT="${OUTROOT}/${STRAIN}/truvari"
    conda_activate truvari

    local SV_TRUTH="${COMPARISON_DIR}/${STRAIN}.svim.svs.vcf.gz"
    local ALT_SV_TRUTH="${OUTROOT}/${STRAIN}/${STRAIN}.svim.svs.alt_only.mixed.vcf.gz"
    if [[ ! -f "${ALT_SV_TRUTH}" ]]; then
        bcftools view -i 'GT!="ref"' "${SV_TRUTH}" -Oz -o "${ALT_SV_TRUTH}"
        bcftools index -f -t "${ALT_SV_TRUTH}"
    fi

    declare -A SV_CALLERS
    SV_CALLERS["pan_direct_sv"]="${COMPARISON_DIR}/${STRAIN}.pan.svs.vcf.gz"
    SV_CALLERS["delly_loo"]="${LOO_BASE}/loo_${STRAIN}/${STRAIN}.loo.delly.SV.ge50bp.PASS.vcf.gz"
    SV_CALLERS["vg_hap_sv_loo"]="${LOO_BASE}/loo_${STRAIN}/${STRAIN}.loo.pan.SV.ge50bp.PASS.hap.vcf.gz"
    SV_CALLERS["delly_linear"]="${OUTROOT}/${STRAIN}/delly_linear/${STRAIN}.delly_linear.SV.ge50bp.PASS.vcf.gz"

    for CALLER in "${!SV_CALLERS[@]}"; do
        local CALLS="${SV_CALLERS[$CALLER]}"
        [[ -f "${CALLS}" ]] || { echo "  SKIP ${CALLER}"; continue; }
        local ALT_CALLS="${CALLS%.vcf.gz}.alt_only.mixed.vcf.gz"
        if [[ ! -f "${ALT_CALLS}" ]]; then
            bcftools view -i 'GT!="ref"' "${CALLS}" -Oz -o "${ALT_CALLS}"
            bcftools index -f -t "${ALT_CALLS}"
        fi
        for BED in "${BEDDIR}"/*.bed; do
            local BEDNAME; BEDNAME=$(basename "${BED}" .bed)
            local BP; BP=$(awk '{s+=$3-$2}END{print s+0}' "${BED}")
            [[ "${BP}" -lt "${MIN_BP}" ]] && continue
            local SAMPLE_OUT="${OUT}/${CALLER}/${BEDNAME}"
            mkdir -p "$(dirname "${SAMPLE_OUT}")"; rm -rf "${SAMPLE_OUT}"
            echo "  -> [${STRAIN}] ${CALLER} | ${BEDNAME}"
            truvari bench -b "${ALT_SV_TRUTH}" -c "${ALT_CALLS}" -o "${SAMPLE_OUT}" -f "${REF}" \
                -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 10000 --includebed "${BED}"
        done
    done
    unset SV_CALLERS
}

############################################################
# MAIN
############################################################
build_shared_beds

for STRAIN in "${STRAINS[@]}"; do
    echo ""
    echo "##########################################"
    echo "# ${STRAIN}"
    echo "##########################################"

    if [[ "${RUN_MAPPING}" == "true" ]]; then
        run_loo "${STRAIN}"
    else
        echo "[skip] mapping (RUN_MAPPING=false)"
    fi

    if [[ "${RUN_SMALL_VAR}" == "true" ]]; then
        run_gatk_linear "${STRAIN}"
        run_gatk_surject_loo "${STRAIN}"
        run_vg_hap_loo_filter "${STRAIN}"
        build_confidence "${STRAIN}"      # builds FINAL (short-var) confidence tier
        run_vcfeval_all "${STRAIN}"
    else
        echo "[skip] small-variant filtering + vcfeval (RUN_SMALL_VAR=false)"
        # no confidence-tier bed is needed for SVs anymore (plain category beds used
        # directly in run_truvari_all), so nothing to backfill here
    fi

    if [[ "${RUN_SV}" == "true" ]]; then
        run_delly_loo "${STRAIN}"
        run_delly_linear "${STRAIN}"
        run_truvari_all "${STRAIN}"
    else
        echo "[skip] structural validation (RUN_SV=false)"
    fi
done

echo "DONE. Results in ${OUTROOT}/<strain>/vcfeval and .../truvari"
