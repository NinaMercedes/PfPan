#!/usr/bin/env bash
#
# QC filter for the merged pangenome SV callset (vg call, haploid samples).
#
#   usage:  ./qc_filter.sh in.vcf.gz out.vcf.gz
#           ./qc_filter.sh in.vcf.gz out.vcf.gz --dry-run    # counts only
#
# Why the ORDER matters
# ---------------------
# vg call emits a genotype even when no read supports either allele: the two
# genotype likelihoods come out exactly tied and the tie breaks toward the ALT,
# giving GT=1 with DP=0, AD=0,0 and GQ=0. Those calls look "called" to anything
# that only reads GT, so F_MISSING and AF computed on the raw merged file are
# wrong.
#
# So: blank the unsupported genotypes FIRST, then recompute F_MISSING/AC/AN/AF
# from what is left, then apply the site filters.
#
# Three traps in the obvious version of this filter
# -------------------------------------------------
#   * INFO/DP is the COHORT SUM, not per-sample depth. Across 878 samples a
#     threshold of 5 is met by 5 reads spread over the whole cohort, so it is
#     effectively inert. Use FMT/DP instead.
#   * FILTER is NOT usable on this file. It was produced by `bcftools merge`,
#     which defaults to --filter-logic '+' and applies every input's filters.
#     A site is PASS only if it passed in all 878 samples, so almost nothing is
#     (11 records genome-wide). Per-sample blanking below does this job properly.
#   * ILEN is empty on this callset (complex REF/ALT pairs, not simple indels)
#     and SVLEN is populated for only a subset, so size comes from the allele
#     strings via strlen().
#
set -euo pipefail

IN=${1:?usage: qc_filter.sh in.vcf.gz out.vcf.gz [--dry-run]}
OUT=${2:?usage: qc_filter.sh in.vcf.gz out.vcf.gz [--dry-run]}
DRY=${3:-}

# ---- thresholds -------------------------------------------------------------
MIN_DP=5         # per-sample read depth below this -> genotype set to missing
MIN_GQ=20        # per-sample genotype quality below this -> missing
MAX_MISS=0.2     # max fraction of missing genotypes AFTER blanking
MIN_QUAL=30      # site quality
MIN_AC=2         # drop singletons
MAX_LEN=10000    # drop events >= this; not reliably genotyped from short reads

command -v bcftools >/dev/null || { echo "bcftools not found" >&2; exit 1; }

SITE_EXPR="QUAL>=${MIN_QUAL} && F_MISSING<=${MAX_MISS} && AC>=${MIN_AC} && AN>0 \
 && strlen(REF)-strlen(ALT) < ${MAX_LEN} && strlen(ALT)-strlen(REF) < ${MAX_LEN}"

# ---- before ----------------------------------------------------------------
echo "=== input: ${IN}"
N_SAMP=$(bcftools query -l "$IN" | wc -l)
N_IN=$(bcftools view -H "$IN" | wc -l)
printf '  samples   %d\n' "$N_SAMP"
printf '  records   %d\n' "$N_IN"
echo
echo "  (FILTER column is a merge union and is deliberately NOT used - see notes)"
echo

if [[ "$DRY" == "--dry-run" ]]; then
  echo "=== would apply"
  echo "  blank genotypes with FMT/DP<${MIN_DP} or FMT/GQ<${MIN_GQ}"
  echo "  then keep sites matching:"
  echo "    ${SITE_EXPR}" | tr -s ' '
  echo
  echo "(dry run - no output written)"
  exit 0
fi

# ---- filter ----------------------------------------------------------------
echo "=== filtering"
echo "  blanking genotypes with FMT/DP<${MIN_DP} or FMT/GQ<${MIN_GQ}"
echo "  site expression:"
echo "    ${SITE_EXPR}" | tr -s ' '
echo

bcftools +setGT "$IN" -- -t q -i "FMT/DP<${MIN_DP} | FMT/GQ<${MIN_GQ}" -n . \
| bcftools +fill-tags -- -t F_MISSING,AC,AN,AF \
| bcftools filter -i "$SITE_EXPR" -Oz -o "$OUT"

bcftools index -t "$OUT"

# ---- after -----------------------------------------------------------------
N_OUT=$(bcftools view -H "$OUT" | wc -l)
echo
echo "=== output: ${OUT}"
printf '  records kept   %d / %d  (%.2f%%)\n' "$N_OUT" "$N_IN" \
       "$(awk -v a="$N_OUT" -v b="$N_IN" 'BEGIN{print (b?100*a/b:0)}')"
if [[ "$N_OUT" -gt 0 ]]; then
  echo
  echo "  size distribution of what survived (ALT-REF, bp):"
  bcftools view -H "$OUT" \
  | awk -F'\t' '{d=length($5)-length($4); if(d<0)d=-d;
      if(d==0)a++; else if(d<50)b++; else if(d<1000)c++; else if(d<10000)e++}
      END{printf "    SNV/MNP   %d\n    1-49bp    %d\n    50-999bp  %d\n    1-10kb    %d\n",a,b,c,e}'
fi
echo
echo "done."