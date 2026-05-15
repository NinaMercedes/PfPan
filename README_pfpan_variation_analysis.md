# PfPan Variation Analysis

This section describes the analysis of genetic variation captured in the PfPan pangenome graph, covering variant extraction and annotation from the Cactus VCF, genomic characterisation, identification of genes of interest, and visualisation of haplotype relationships. The analysis is divided into two stages: command-line preprocessing and R/Python-based statistical analysis and plotting.

---

## Overview

```
Stage 1 — VCF preprocessing (bash/bcftools)
    Normalise → annotate SVs → filter → functional annotation → GOI extraction → graph visualisation

Stage 2 — Statistical analysis and plotting
    analysis2.R     — variant processing, regional analysis, GOI analysis, complex loci
    analysis.R      — SV genomic distribution, fold enrichment, BED export
    PfPan_pca_network.py — PCA, haplotype similarity network, summary figure
```

---

## Input Files

| File | Description |
|------|-------------|
| `PfPan.vcf.gz` | Raw Cactus pangenome VCF (multi-sample, projected onto Pf3D7) |
| `Pfalciparum.genome.fasta` | Pf3D7 v3 linear reference FASTA |
| `Pfalciparum.genome.modified.new.gff3` | Pf3D7 v3 gene annotation GFF3 |
| `watchlist.txt` | List of genes of interest (one gene ID per line) |
| `samples.txt` | Sample names (one per line, matches VCF sample order) |
| `*.bed` | Genomic region BED files (core genome, coding, centromere, hypervariable, subtelomeric) |

---

## Stage 1 — VCF Preprocessing

### 1.1 Normalise and annotate SVs

```bash
# Normalise multiallelic sites
bcftools norm -m -any -Oz -o PfPan.norm.vcf.gz PfPan.vcf.gz
tabix PfPan.norm.vcf.gz

# Annotate with SV type and length using truvari svinfo
# This adds SVTYPE and SVLEN INFO fields to sequence-resolved variants
# from vg call which otherwise lack symbolic allele annotations
truvari anno svinfo -o PfPan.norm.svinfo.vcf.gz PfPan.norm.vcf.gz
tabix PfPan.norm.svinfo.vcf.gz
```

---

### 1.2 Inspect headers and extract variant table

```bash
# Check available INFO and FORMAT fields
bcftools view -h PfPan.norm.svinfo.vcf.gz | grep "##INFO"   | cut -d',' -f1 | cut -d'=' -f3
bcftools view -h PfPan.norm.svinfo.vcf.gz | grep "##FORMAT" | cut -d',' -f1 | cut -d'=' -f3

# Extract variant info + per-sample genotypes to TSV
bcftools query \
  -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%AC\t%AN\t%AF\t%NS\t%SVTYPE\t%SVLEN[\t%GT]\n' \
  PfPan.norm.svinfo.vcf.gz > PfPan.all_variants_info_GT.tsv
```

This TSV is the primary input for all downstream R and Python analysis.

---

### 1.3 Quality filter

```bash
# Keep variants present in >1 sample, polymorphic (0 < AF < 1),
# and with absolute allele length <= 10kb
bcftools view \
  -i 'INFO/NS>1 && INFO/AF>0 && INFO/AF<1 && abs(strlen(ALT)-strlen(REF)) <= 10000' \
  PfPan.norm.svinfo.vcf.gz \
  -O z -o PfPan.filtered.vcf.gz
tabix -p vcf PfPan.filtered.vcf.gz
```

---

### 1.4 Functional annotation with bcftools csq

```bash
bcftools csq \
  -f Pfalciparum.genome.fasta \
  -g Pfalciparum.genome.modified.new.gff3 \
  -o PfPan.filtered.ann.vcf.gz -O z \
  PfPan.filtered.vcf.gz
```

Adds `BCSQ` FORMAT field with consequence annotations (missense, synonymous, etc.) for variants overlapping coding sequences.

---

### 1.5 Genes of interest (GOI) extraction

Extracts a subset of genes from the annotation using a watchlist, builds extended gene body BED files (with promoter and downstream flanks), and extracts overlapping variants.

```bash
# Extract GOI features from GFF
grep -F -f watchlist.txt Pfalciparum.genome.modified.new.gff3 > GOI.gff3

# Build BED with 2kb upstream promoter and 500bp downstream flank
awk -v prom=2000 -v down=500 'BEGIN{OFS="\t"} $3=="gene" {
    match($9, /ID=([^;]+)/, a); id=a[1];
    if ($7=="+") { start=$4-prom; if (start<0) start=0; end=$5+down; }
    else         { start=$4-down; if (start<0) start=0; end=$5+prom; }
    print $1, start-1, end, id, ".", $7;
}' GOI.gff3 > genes_of_interest.bed

# Index and extract GOI variants
bcftools index -t PfPan.filtered.ann.vcf.gz
bcftools view -R genes_of_interest.bed -Oz -o PfPan.filtered.ann.GOI.vcf.gz PfPan.filtered.ann.vcf.gz
bcftools index -t PfPan.filtered.ann.GOI.vcf.gz

# Export annotated tables
bcftools query \
  -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%AC\t%AN\t%AF\t%NS\t%SVTYPE\t%SVLEN\t%BCSQ[\t%GT]\n' \
  PfPan.filtered.ann.vcf.gz > annotated_table.tsv

bcftools query \
  -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%AC\t%AN\t%AF\t%NS\t%SVTYPE\t%SVLEN\t%BCSQ[\t%GT]\n' \
  PfPan.filtered.ann.GOI.vcf.gz > GOI_table.tsv

bcftools query -l PfPan.filtered.ann.GOI.vcf.gz > haplos.txt
```

---

### 1.6 Graph region extraction and visualisation

Regions of interest are extracted from the pangenome graph for visualisation in the sequence tube map (https://graph-genome.github.io/) and as ODGI images. The example shown is the *pfcrt* region on chromosome 7.

```bash
# Extract region from full ODGI graph
odgi extract \
  -i PfPan.full.og \
  -r Pf3D7#0#Pf3D7_07_v3:401221-406817 \
  -o crt_region.og

# Generate ODGI linear visualisation
odgi viz -i crt_region.og -o crt_region.png

# Export to GFA for sequence tube map
odgi view -i crt_region.og -g > crt_region.gfa

# Alternatively extract as VG subgraph for sequence tube map
conda activate cactus
vg index -x chr07.xg pf_pan_2024.chroms/Pf3D7_07_v3.vg
vg snarls chr07.xg > chr07.snarls
vg chunk \
  -x chr07.xg \
  -p Pf3D7#0#Pf3D7_07_v3:401221-406817 \
  -T -O vg -S chr07.snarls > crt_region.vg
```

The resulting `.gfa` or `.vg` file can be loaded into the web-based sequence tube map to visualise haplotype paths through the region, as shown in the *pfcrt* example figure.

---

## Stage 2 — Statistical Analysis and Plotting

### 2.1 Variant processing and analysis (`PfPan_variation.R`)

**Required R packages:** `dplyr`, `ggplot2`, `GenomicRanges`, `tidyr`, `stringr`, `readr`, `purrr`, `scales`

**Input files:**
- `PfPan.all_variants_info_GT.tsv` — extracted from Stage 1.2
- `samples.txt` — sample names
- `genes_of_interest.bed` — from Stage 1.5
- `hrp.bed` — HRP2/3 gene regions
- `Pfalciparum.genome.modified.new.gff3` — gene annotation
- `*.bed` — genomic region BED files (core genome, coding, centromere, hypervariable, subtelomeric)

Input paths are set in the `CONFIGURATION` block at the top of the script.

The script runs 15 sections in order:

**Section 2 — Variant processing:** The `process_variants()` function filters and classifies the raw TSV:
- Removes variants with >90% missing genotypes across samples
- Removes fixed variants (AF = 0)
- Removes variants with absolute allele length >10kb
- Classifies variants by type: SNP (LEN=0), Small Insertion/Deletion (1–50bp), Structural Insertion/Deletion (>50bp)

**Sections 7–8 — GOI and HRP analysis:** `analyse_GOI()` overlaps variants with extended gene bodies (2kb upstream + 500bp downstream) for genes in the watchlist and the HRP locus separately.

**Section 9 — Complex loci:** `analyse_complex_loci()` clusters nearby variant sites into windows, identifies loci with high allelic complexity spanning >1kb, and annotates overlapping genes from the GFF.

**Sections 10–12 — Regional analysis:** Variants are mapped to genomic regions (core genome, coding regions, centromere, internal hypervariable, sub-telomeric hypervariable, sub-telomeric repeat) using `GenomicRanges::findOverlaps()`.

**Sections 13–15 — SV genomic distribution:** Bins SVs into 10kb and 100kb windows, calculates fold enrichment over the per-chromosome mean, and exports a BED file for integration with genome browsers.

**Key outputs:**

| File | Description |
|------|-------------|
| `all_variants_variant_type_counts.csv` | Count of each variant type |
| `all_variants_AF_summary.csv` | Mean and median AF per variant type |
| `all_variants_SV_counts.csv` | SV insertion and deletion counts |
| `all_variants_singletons.csv` | Variant counts by type and allele count |
| `all_variants_<type>_length_stats.csv` | Median/min/max length per variant type |
| `variant_type_counts_by_region_summary.csv` | Variant type counts and percentages per genomic region |
| `GOI_summary.csv` | Variant counts per GOI gene and type |
| `HRP_summary.csv` | Variant counts per HRP gene and type |
| `pf_complex_summary.csv` | Complex loci with allele counts, span, and overlapping genes |
| `variant_types.tiff` | Bar chart of variant type counts |
| `sv_length_distribution.tiff` | SV length histogram by type |
| `AF_distribution.tiff` | Allele frequency distribution by variant type |
| `multipanel_signed_SV_by_region.tiff` | Signed SV length distributions faceted by genomic region |
| `multipanel_variant_types_by_region.tiff` | Variant type counts faceted by genomic region |
| `multipanel_AF_by_region.tiff` | AF distributions faceted by genomic region |
| `multipanel_variant_types_counts_by_region.tiff` | Stacked bar chart of variant type counts per region |
| `sv_density_by_chr.pdf/.png` | SV density histogram faceted by chromosome and SV type |
| `sv_length_dist.png` | SV length distribution on log scale |
| `sv_fold_enrichment.pdf/.png` | SV fold enrichment per 10kb window (line plot) |
| `sv_fold_enrichment_coverage.pdf/.png` | SV fold enrichment per 10kb window (area plot, per chromosome) |
| `sv_enrichment.bed` | 100kb window BED with fold enrichment classification (low/mid/high/hotspot) |
| `sv_enrichment.idmap` | Colour map for BED enrichment classes |

---

### 2.3 PCA, network and summary figure (`PfPan_pca_network.py`)

**Required Python packages:** `pandas`, `numpy`, `matplotlib`, `sklearn`, `networkx`, `adjustText`

**Input:** `variants` dataframe loaded from `PfPan.all_variants_info_GT.tsv` (as processed in Stage 2.1)

Performs PCA and Jaccard similarity network analysis on the SV genotype matrix (SVs >50bp, present in ≥2 samples). Missing genotypes are imputed as 0.5. Produces a five-panel summary figure:

| Panel | Description |
|-------|-------------|
| A | Bar chart of variant type counts |
| B | SV length histogram (insertions vs deletions, >50bp) |
| C | Allele frequency distribution by variant type |
| D | PCA of SV genotype matrix (PC1 vs PC2) with inset scree plot and SE Asian cluster annotation |
| E | Haplotype similarity network — nodes sized by alt SV count, edges weighted by Jaccard similarity |

**Key outputs:**

| File | Description |
|------|-------------|
| `sv_summary_figure.png` | Five-panel summary figure (200 dpi) |
| `sv_summary_figure.tiff` | Five-panel summary figure (300 dpi) |

---

## Required Software

| Tool | Used for |
|------|---------|
| `bcftools` | Normalisation, filtering, annotation, querying |
| `truvari` | `anno svinfo` — SV type/length annotation |
| `odgi` | Graph region extraction and visualisation |
| `vg` | Graph indexing and chunk extraction |
| R (≥4.0) | Statistical analysis and plotting |
| Python (≥3.8) | PCA, network, summary figure |

