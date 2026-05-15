# Pf3k Variant Calling and Benchmarking Pipeline

This section covers the full variant calling and benchmarking workflow for *Plasmodium falciparum* Illumina short-read data, comparing linear reference-based and pangenome-based approaches for both short variants (SNPs/indels) and structural variants (SVs). The pipeline is split across four Python scripts which should be run in the order described below.

---

## Overview

```
PfPan_SVIM_ASM.py                  ← assembly-based truth set (run once)
PfPan_linear_map_and_call.py       ← linear BWA mapping + DELLY (per cohort)
PfPan_map_and_call.py              ← pangenome mapping + vg call + DELLY (per sample)
PfPan_variant_evaluation.py        ← filtering, benchmarking, plotting
```

The linear and pangenome calling pipelines are independent and can be run in parallel. Both feed into `PfPan_variant_evaluation.py` for joint benchmarking.

---

## Samples

14 clinical isolates and lab strains from the Pf3k dataset:

| Sample | Origin |
|--------|--------|
| Pf7G8 | Lab strain |
| PfCD01 | Cambodia |
| PfDd2 | Lab strain (Indochina) |
| PfGA01 | Gabon |
| PfGB4 | Lab strain |
| PfGN01 | Guinea |
| PfHB3 | Lab strain |
| PfIT | Lab strain |
| PfKE01 | Kenya |
| PfKH01 | Cambodia |
| PfKH02 | Cambodia |
| PfML01 | Mali |
| PfSN01 | Senegal |
| PfTG01 | Togo |

Note: only 12 were included in PfPan.

---

## Pipeline 0 — Assembly-based Truth Set (`PfPan_SVIM_ASM.py`)

**Conda environment:** requires minimap2, svim-asm, paftools, samtools, bcftools

Run once to generate the assembly-based truth callsets used for benchmarking. For each strain, the PacBio assembly is aligned to Pf3D7 and SVs and SNPs are called independently of Illumina reads.

```bash
python ./pipeline/PfPan_SVIM_ASM.py

# Single strain only
python ./pipeline/PfPan_SVIM_ASM.py --strain PfDd2

# Skip FASTA copy/decompress if already present locally
python ./pipeline/PfPan_SVIM_ASM.py --skip-copy
```

**Steps:**
1. Align each assembly to Pf3D7 with `minimap2` (asm5 preset) → sorted BAM
2. Call SVs from the BAM using `svim-asm` (haploid mode)
3. Align each assembly to Pf3D7 with `minimap2` (PAF output + cs string)
4. Call SNPs/indels from the PAF using `paftools.js call`, normalised with `bcftools`

**Key outputs per strain (`svim_asm_results/<strain>/`):**

| File | Description |
|------|-------------|
| `<strain>.vs.Pf3D7.sorted.bam` | Assembly-to-reference BAM |
| `svim/variants.vcf` | Raw SVIM-asm SV calls |
| `<strain>.paftools.snps.vcf.gz` | Normalised paftools SNP/indel calls — used as short variant baseline in vcfeval |

The SVIM-asm VCFs are used as the SV truth baseline in Truvari benchmarking. The paftools VCFs are used as the short variant truth baseline in RTG vcfeval.

---

## Pipeline 1a — Linear Reference Mapping (`PfPan_linear_map_and_call.py`)

**Conda environments:**
- Steps 1: fastq2matrix environment (for `fastq2vcf.py`)
- Steps 2–5: `delly` environment

Maps trimmed Illumina reads to the Pf3D7 linear reference using BWA, performs BQSR, and calls variants with GATK. Runs a multi-sample DELLY SV calling workflow across the whole cohort.

```bash
# Step 1 only (fastq2matrix env)
conda activate <fastq2matrix_env>
python ./pipeline/PfPan_linear_map_and_call.py --samples-file fastqs.txt --skip-delly

# Steps 2-5 only (delly env)
conda activate delly
python ./pipeline/PfPan_linear_map_and_call.py --samples-file fastqs.txt --skip-fastq2vcf
```

**Steps:**
1. `fastq2vcf.py all` — BWA mapping, duplicate marking, BQSR, GATK HaplotypeCaller per sample (sequential, 1 sample at a time)
2. `delly call` — per-sample SV calling against Pf3D7 (parallel, 10 samples at a time)
3. `delly merge` — merge all per-sample BCFs into a unified sites file
4. `delly call -v sites.bcf` — re-genotype all samples at merged sites (parallel) — ensures all samples are called at the same loci for joint analysis
5. `bcftools merge` — merge all re-genotyped BCFs into a single multi-sample VCF

**Key outputs (`./linear/`):**

| File | Description |
|------|-------------|
| `<sample>.bqsr.bam` | Mapped, recalibrated BAM |
| `<sample>_gatk.g.vcf.gz` | GATK HaplotypeCaller GVCF |
| `<sample>_delly_sites.bcf` | Per-sample DELLY calls re-genotyped at merged sites |
| `sites.bcf` | Merged DELLY sites across all samples |
| `pan_delly.vcf.gz` | Final multi-sample DELLY SV VCF |

---

## Pipeline 1b — Pangenome Mapping (`PfPan_map_and_call.py`)

**Conda environment:** `cactus` (or whichever env has vg, samtools, gatk)

Maps trimmed Illumina reads to the pangenome graph using `vg giraffe`, surjects back to the Pf3D7 linear reference, and calls variants using both GATK (short variants) and vg call (graph-based variants, diploid and haploid).

```bash
python ./pipeline/PfPan_map_and_call.py --sample PfDd2
python ./pipeline/PfPan_map_and_call.py --sample PfDd2 --skip-delly
python ./pipeline/PfPan_map_and_call.py --sample PfDd2 --skip-haploid
```

**Steps:**
1. `vg giraffe` — map trimmed reads to pangenome graph → GAF
2. `vg surject` — project graph alignments onto Pf3D7 paths → BAM
3. `samtools sort/index/flagstat` — sort, index, QC the BAM
4. `gatk HaplotypeCaller` — short variant calling in GVCF mode from the linear BAM
5. `delly call` — SV calling from the surjected BAM
6. `vg pack` — build coverage index from graph alignments
7. `vg call` (diploid) — call variants directly from the graph (default ploidy=2)
8. `vg call --ploidy 1` (haploid) — call variants from the graph in haploid mode, appropriate for the *P. falciparum* blood stage

**Key outputs per sample (`./output/`):**

| File | Description |
|------|-------------|
| `<sample>.gaf.gz` | Graph alignments (vg giraffe) |
| `<sample>_sort.bam` | Sorted, indexed surjected BAM |
| `<sample>_stat.txt` | samtools flagstat summary |
| `<sample>_gatk.g.vcf.gz` | GATK HaplotypeCaller GVCF |
| `<sample>_delly_sites.bcf` | DELLY SV calls |
| `<sample>.pack` | vg pack coverage index |
| `<sample>.SV.call.vcf.gz` | vg call diploid variant calls |
| `<sample>.SV.call_haploid_.vcf.gz` | vg call haploid variant calls |

---

## Pipeline 2 — Variant Filtering, Benchmarking and Evaluation (`PfPan_variant_evaluation.py`)

**Conda environment:** requires bcftools, gatk, rtg-tools, truvari, R

Filters and prepares all variant callsets, runs benchmarking against assembly-based truth sets, and generates summary plots and tables.

```bash
# Full pipeline
python ./pipeline/PfPan_variant_evaluation.py

# Individual steps
python ./pipeline/PfPan_variant_evaluation.py --steps filter_delly_gatk filter_pangenome_vcfs
python ./pipeline/PfPan_variant_evaluation.py --steps prebench decompose vcfeval truvari plot

# Skip steps
python ./pipeline/PfPan_variant_evaluation.py --skip recode_haploid_genotypes plot
```

**Steps:**

| Step | Description |
|------|-------------|
| `filter_delly_gatk` | Filter DELLY BCFs to PASS SVs ≥50bp; run full GATK hard-filter pipeline on gVCFs → normalised SNPs + indels <50bp with genotype recoding via `setGT.py` (AF threshold 0.7) |
| `filter_pangenome_vcfs` | Filter diploid vg call VCFs into short variants (<50bp) and SVs (≥50bp) |
| `filter_pangenome_vcfs_haploid` | Same as above for haploid vg call VCFs |
| `recode_haploid_genotypes` | Recode haploid GTs (0→0/0, 1→1/1) for compatibility with vcfeval and truvari |
| `prebench` | Rename SVIM-asm VCFs, filter to SVs ≥50bp; split pangenome SV VCFs into SNP/indel and SV subsets (hom-alt only); copy paftools baselines |
| `decompose` | Decompose MNPs and complex indels using `rtg vcfdecompose` |
| `vcfeval` | RTG vcfeval short variant benchmarking against paftools baselines, across all BED regions and callsets. Applies DP>5 filter to pangenome and GATK callsets before filtering to hom-alt (1/1) genotypes |
| `truvari` | Truvari bench SV benchmarking against SVIM-asm truth, across all BED regions and callsets. Filters to ALT-only calls (keeping mixed genotypes) |
| `plot` | Calls `plot_evaluation.R` to generate F-measure boxplot figures and summary CSVs |

**Callsets benchmarked:**

| Callset | Variant type | Tool |
|---------|-------------|------|
| `gatk_shortvars` | SNPs + indels <50bp | GATK (linear mapping) |
| `pan_snps` | SNPs + indels | Pangenome graph (Cactus VCF) |
| `pan_shortvars` | SNPs + indels <50bp | vg call diploid |
| `pan_shortvars_hap` | SNPs + indels <50bp | vg call haploid |
| `delly` | SVs ≥50bp | DELLY (linear mapping) |
| `pan_SV` | SVs ≥50bp | vg call diploid |
| `pan_SV_hap` | SVs ≥50bp | vg call haploid |
| `pan_svs` | SVs ≥50bp | Pangenome graph (Cactus VCF) |

**Genomic regions evaluated (BED files):**

| Region | Description |
|--------|-------------|
| Core Genome | High-confidence mappable regions |
| Coding Regions | Annotated coding sequences |
| Centromeric | Centromere regions |
| Internal Hypervariable | Internal hypervariable gene families (e.g. var, rifin) |
| Sub-telomeric Hypervariable | Subtelomeric hypervariable regions |
| Sub-telomeric Repeat | Subtelomeric repeat regions |

**Key outputs:**

| File | Description |
|------|-------------|
| `vcfeval_raw_strain_data.csv` | Per-strain precision, sensitivity, F-measure for all callsets and regions |
| `vcfeval_summary_shortvars.csv` | Aggregated summary statistics for short variants |
| `vcfeval_summary_svs.csv` | Aggregated summary statistics for SVs |
| `vcfeval_ab_all_regions.pdf/.png` | Combined panel figure — short variants (panel a) and SVs (panel b) across all regions |
| `vcfeval_a_shortvars.pdf/.png` | Short variant F-measure panel only |
| `vcfeval_b_svs.pdf/.png` | SV F-measure panel only |
| `vcfeval_ab_core_genome.pdf/.png` | Core genome only |
| `vcfeval_ab_coding_regions.pdf/.png` | Coding regions only |

---

## Required Reference Files

| File | Description |
|------|-------------|
| `Pfalciparum.genome.fasta` | Pf3D7 v3 linear reference FASTA — hardcoded at `/mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/` |
| `Pfalciparum.genome.fasta.sdf` | RTG sequence dictionary (see generation command below) |
| `PfPan_Pf3D7_pan/PfPan_test.gbz` | Pangenome graph (from pangenome construction pipeline) |
| `PfPan_Pf3D7_pan/PfPan_test.snarls` | Snarl index for vg call |
| `Pf3D7.paths.txt` | List of Pf3D7 reference paths for vg surject (see generation command below) |
| `bed_files/*.bed` | Genomic region BED files — hardcoded at `/mnt/storage13/nbillows/pangenome/analysis/bed_files/` |

The following files must be generated once before running the pangenome calling pipeline:

```bash
# Extract Pf3D7 path names from the graph for vg surject
vg paths -x ./PfPan_Pf3D7_pan/PfPan_test.gbz -L | grep Pf3D7 > Pf3D7.paths.txt

# Build RTG sequence dictionary for vcfeval
rtg format -o Pfalciparum.genome.fasta.sdf Pfalciparum.genome.fasta
```

