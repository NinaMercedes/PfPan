# PfPan Short-Read Population Analysis

This section describes the population-level analysis of *P. falciparum* short-read sequencing data aligned to the PfPan pangenome graph and the Pf3D7 linear reference. The analysis covers variant filtering, principal component analysis, differentiation statistics, copy number variation, and structural variant characterisation.

---

## Overview

```
Stage 1 — VCF filtering and harmonisation (bcftools)
    Multi-sample SNP/indel VCFs (linear and pangenome) → quality-filtered biallelic SNP matrices

Stage 2 — PCA comparison (PfPan_pca_plink.R)
    Plink distance matrices → PCoA → compare linear vs pangenome SNP clustering

Stage 3 — Jost's D differentiation analysis (PfPan_josts_d.py)
    Population differentiation across genomic regions and watchlist genes

Stage 4 — CNV calling (PfPan_cnv_call.py)
    mosdepth coverage → flanking-region normalised amplification detection
    (MDR1, Plasmepsin II/III)

Stage 5 — SV characterisation (PfPan_sv_pipeline.sh)
    Core genome SVs ≥2kb → MAFFT sequence extraction → BLAST → gene annotation
    → whole/partial gene deletion detection → summary report
```

---

## Input Files

| File | Description |
|------|-------------|
| `pan_SR_GTK.bi.vcf.gz` | Biallelic SNP VCF — pangenome-aligned short reads (GATK) |
| `norm_SR_GTK.bi.vcf.gz` | Biallelic SNP VCF — linear reference-aligned short reads (GATK) |
| `pan_SR_GTK.bi.dist` / `.dist.id` | Plink SNP distance matrix — pangenome |
| `norm_SR_GTK.bi.dist` / `.dist.id` | Plink SNP distance matrix — linear reference |
| `merged_diploid.vcf.gz` | Merged diploid VCF for Jost's D analysis |
| `metadata_high_quality.csv` | Sample metadata with `sample`, `Region`, `Country` columns |
| `PfPan.norm.svinfo.vcf.gz` | Normalised pangenome VCF with svinfo annotations |
| `Core_genome_Pf3D7_v3_ext.bed` | Core genome BED file |
| `Pfalciparum.genome.fasta` | Pf3D7 v3 reference FASTA |
| `Pfalciparum.genome.modified.new.gff3` | Gene annotation GFF3 |
| `watchlist.txt` | Genes of interest for Jost's D prioritisation |
| `mosdepth_output/` | Per-sample mosdepth coverage files (`*.regions.bed.gz`) |

---

## Stage 1 — VCF Filtering and Harmonisation

Multi-sample SNPs and indels were called jointly with GATK across short-read samples aligned to both the linear reference (Pf3D7 v3) and the PfPan pangenome graph. The filtering history is recorded in the VCF headers and was applied as follows:

**Quality filters applied:**
- `QUAL < 20` — low-confidence variant calls removed
- `MQ < 40` — poor mapping quality sites removed
- `QD < 2` — low quality-by-depth removed
- `FS > 60` — strand bias filter
- `BaseQRankSum` and `MQRankSum` outside ±12.5 — rank sum outliers removed
- `AF == 0` — invariant sites removed
- `FMT/DP < 5` — individual genotypes with fewer than 5 reads set to missing
- `INFO/DP < 5` — sites with low total depth removed

**Post-filter steps:**
- Contig names harmonised between pangenome (e.g. `Pf3D7#0#Pf3D7_01_v3`) and linear reference (e.g. `Pf3D7_01_v3`) using `bcftools annotate --rename-chrs`
- SNPs and indels merged with `bcftools concat`
- Normalised with `bcftools norm -m -both`
- Subset to high-quality samples
- Filtered to biallelic SNPs only for PCA: `bcftools view -m2 -M2 -v snps`

**Plink distance matrix generation:**

The biallelic SNP VCF was converted to Plink format and a genome-wide SNP distance matrix generated for PCoA:

```bash
# Pangenome SNP distance matrix
plink --vcf pan_SR_GTK.bi.vcf.gz \
      --distance square \
      --double-id \
      --allow-extra-chr \
      --out pan_SR_GTK.bi

# Linear reference SNP distance matrix
plink --vcf norm_SR_GTK.bi.vcf.gz \
      --distance square \
      --double-id \
      --allow-extra-chr \
      --out norm_SR_GTK.bi
```

---

## Stage 2 — PCA Comparison (`PfPan_pca_plink.R`)

**Required R packages:** `dplyr`, `ggplot2`, `ape`, `patchwork`

Performs principal coordinates analysis (PCoA) on the Plink SNP distance matrices from both the pangenome and linear reference pipelines and produces a panelled comparison figure. Samples are coloured by geographic region from `metadata_high_quality.csv`.

```bash
Rscript PfPan_pca_plink.R
```

The script reads both distance matrices (prefixes `pan_SR_GTK.bi` and `norm_SR_GTK.bi`), computes PCoA with `cmdscale()`, and plots PC1 vs PC2 side by side (Panel A = pangenome, Panel B = linear) to compare how the two calling strategies resolve population structure.

**Key output:**

| File | Description |
|------|-------------|
| `PCA_panelled_AB.tiff` | Side-by-side PCoA — pangenome vs linear reference SNP calls, coloured by region |

---

## Stage 3 — Jost's D Differentiation Analysis (`PfPan_josts_d.py`)

**Required packages:** `pandas`, `numpy`, `scikit-allel`, `zarr`, `tqdm`

Calculates Jost's D, a true differentiation statistic (unlike FST it is not confounded by within-population diversity), across the genome and specifically within watchlist genes. Runs in parallel across chromosomes.

```bash
python PfPan_josts_d.py \
    --vcf merged_diploid.vcf.gz \
    --metadata metadata_high_quality.csv \
    --gff Pfalciparum.gff3 \
    --watchlist watchlist.txt \
    --output-dir josts_d_results \
    --threads 20 \
    --percentile 99
```

**Key arguments:**

| Argument | Description |
|----------|-------------|
| `--vcf` | Merged diploid VCF |
| `--metadata` | CSV with `sample` and `Region` columns for population assignment |
| `--gff` | GFF3 for gene-level annotation of outlier loci |
| `--watchlist` | Gene IDs to prioritise in output |
| `--percentile` | Top percentile threshold for flagging high-differentiation loci (default 99) |
| `--threads` | Parallel workers |

**Key outputs (`josts_d_results/`):**

| File | Description |
|------|-------------|
| `josts_d_genome.tsv` | Genome-wide Jost's D per variant |
| `josts_d_outliers.tsv` | Variants above the specified percentile threshold |
| `josts_d_watchlist.tsv` | Jost's D values for watchlist gene regions |
| `josts_d_summary.csv` | Per-gene summary statistics |
| `josts_d_plots/` | Genome-wide Manhattan-style plots |

---

## Stage 4 — CNV Calling (`PfPan_cnv_call.py`)

**Required packages:** `pandas`, `numpy`

Detects copy number amplifications at key drug resistance loci (MDR1) using mosdepth depth-of-coverage files. Coverage within each target region is normalised against flanking regions rather than the whole-genome mean, which avoids confounding from subtelomeric low-coverage artefacts. 

Results are expressed as concordance/discordance relative to MalariaGen linear mapping CNV calls, as the MalariaGen calls are used as a reference comparator rather than a gold standard.

**Setup — run mosdepth per sample before calling:**

```bash
for bam in ./output/*_sort.bam; do
    sample=$(basename "$bam" _sort.bam)
    mosdepth \
        --by 500 \
        --threads 4 \
        "${sample}" \
        "${bam}"
done
```

**Run CNV calling:**

```bash
python PfPan_cnv_call.py \
    --mosdepth-dir ./mosdepth_output \
    --output cnv_results.tsv
```

**Quality thresholds:**

| Parameter | Value | Description |
|-----------|-------|-------------|
| `MIN_FLANK_COVERAGE` | 5 | Minimum mean flank depth to trust a ratio |
| `MAX_FLANK_RATIO` | 1.5 | Maximum left/right flank imbalance before flagging a gradient |
| `MIN_LOW_CONF_PATTERNS` | 2 | Patterns required to call amplification at low-confidence loci |

**Key outputs:**

| File | Description |
|------|-------------|
| `cnv_results.tsv` | Per-sample CNV calls with coverage ratios and confidence flags |
| `cnv_concordance.tsv` | Concordance with MalariaGen reference calls |

---

## Stage 5 — SV Characterisation (`PfPan_sv_pipeline.sh`)

**Required tools:** `bcftools`, `mafft`, `blastn`, `makeblastdb`, `bedtools`, `python3`

A 12-step pipeline that characterises large structural variants (≥2kb by default) from the pangenome VCF within the core genome. Because the pangenome VCF uses sequence-resolved alleles rather than symbolic `<DEL>`/`<INS>` tags, SVs are identified by REF/ALT length difference and sequences extracted via MAFFT alignment.

```bash
# Standard run
bash PfPan_sv_pipeline.sh \
    -v PfPan.norm.svinfo.vcf.gz \
    -b Core_genome_Pf3D7_v3_ext.bed \
    -r Pfalciparum.genome.fasta \
    -g Pfalciparum.genome.modified.new.gff3 \
    -o sv_results \
    -t 16 \
    -m 2000

# Or use the provided run script
bash PfPan_sv_pipeline_run.sh

# Restart from after BLAST (skip steps 1-7)
bash PfPan_sv_pipeline.sh [args] -s 8
```

**Pipeline flags:**

| Flag | Default | Description |
|------|---------|-------------|
| `-v` | — | Input VCF (required) |
| `-b` | — | Core genome BED (required) |
| `-r` | — | Reference FASTA (required) |
| `-g` | — | GFF3 annotation (required) |
| `-o` | — | Output directory (required) |
| `-t` | 8 | Threads |
| `-m` | 2000 | Minimum SV length (bp) |
| `-p` | 85 | Minimum BLAST percent identity |
| `-q` | 20 | Minimum BLAST query coverage |
| `-s` | 1 | Start from this step (for restarts) |

**Pipeline steps:**

| Step | Description |
|------|-------------|
| 1 | Filter VCF to core genome using BED file |
| 2 | Classify SVs by REF/ALT length difference into insertions and deletions |
| 3 | Build BLAST nucleotide database from reference FASTA |
| 4 | Extract insertion sequences — align REF vs ALT with MAFFT, extract gap blocks present in ALT |
| 5 | Extract deletion sequences — align REF vs ALT with MAFFT, extract gap blocks present in REF |
| 6 | BLAST insertion sequences against reference genome |
| 7 | BLAST deletion sequences against reference genome |
| 8 | Parse GFF3 to gene BED file |
| 9 | Intersect SV coordinates with gene annotations using bedtools |
| 10 | Classify deletions as whole gene (>90% overlap) or partial (20–90% overlap) |
| 11 | Annotate BLAST results — assign overlapping genes, estimate copy number from bitscore ratio, flag putative tandem duplications (same-chromosome hits) |
| 12 | Generate summary report — flags resistance gene hits (e.g. CRT, MDR1, DHFR, DHPS, GCH1, Kelch13, Plasmepsin II/III) |

**Copy number estimation:** For insertion sequences with BLAST hits, the ratio of total bitscore to maximum single-hit bitscore is used to estimate copy number (`cn = round(total / max)`). An insertion with `cn > 1` and a same-chromosome hit is flagged as a putative tandem duplication.

**Gene flank window:** BLAST hits within 2kb of a gene body are included in the annotation, tagged with `~` to indicate promoter/downstream rather than coding sequence overlap.

**Examples of Resistance genes flagged automatically (Other Targets included):**

| Gene ID | Drug |
|---------|------|
| PF3D7_0709000 | CRT (chloroquine) |
| PF3D7_0523000 | MDR1 (multidrug) |
| PF3D7_0417200 | DHFR (antifolate) |
| PF3D7_0810800 | DHPS (antifolate) |
| PF3D7_1224000 | GCH1 (antifolate) |
| PF3D7_1343700 | Kelch13 (artemisinin) |
| PF3D7_1408000 | Plasmepsin II (piperaquine) |
| PF3D7_1408100 | Plasmepsin III (piperaquine) |

**Key outputs (`sv_results/`):**

| File | Description |
|------|-------------|
| `vcf/core_genome.vcf.gz` | VCF filtered to core genome |
| `vcf/core_insertions.vcf.gz` | Classified insertion VCF |
| `vcf/core_deletions.vcf.gz` | Classified deletion VCF |
| `sequences/insertion_sequences.fa` | Extracted insertion sequences (MAFFT) |
| `sequences/deletion_sequences.fa` | Extracted deletion sequences (MAFFT) |
| `sequences/insertion_summary.tsv` | Per-SV extraction summary |
| `sequences/deletion_summary.tsv` | Per-SV extraction summary |
| `blast/insertions_blast_raw.txt` | Raw BLAST output — insertions |
| `blast/deletions_blast_raw.txt` | Raw BLAST output — deletions |
| `results/insertions_annotated.tsv` | Annotated insertions with gene overlaps and copy number estimates |
| `results/deletions_annotated.tsv` | Annotated deletions with gene overlaps |
| `results/whole_gene_deletions.tsv` | SVs with >90% gene body overlap |
| `results/partial_gene_deletions.tsv` | SVs with 20–90% gene body overlap |
| `results/summary_report.txt` | Plain-text summary with resistance gene flags |
| `logs/pipeline.log` | Full run log |

**Standalone scripts:** `PfPan_extract_sv_sequences.py` and `PfPan_annotate_svs.py` can be run independently outside the pipeline if individual steps need to be re-run or applied to a different VCF.

---

## Script Summary

| Script | Language | Description |
|--------|----------|-------------|
| `PfPan_sv_pipeline.sh` | Bash | Main SV characterisation pipeline (12 steps) |
| `PfPan_sv_pipeline_run.sh` | Bash | Example run command for `PfPan_sv_pipeline.sh` |
| `PfPan_extract_sv_sequences.py` | Python | Standalone SV sequence extraction via MAFFT |
| `PfPan_annotate_svs.py` | Python | Standalone BLAST result annotation |
| `PfPan_josts_d.py` | Python | Parallel Jost's D genome-wide differentiation analysis |
| `PfPan_cnv_call.py` | Python | Flanking-region normalised CNV detection |
| `PfPan_pca_plink.R` | R | PCoA comparison of pangenome vs linear SNP calls |

