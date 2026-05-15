# PfPan — *Plasmodium falciparum* Pangenome Analysis Pipeline

This repository contains the full analysis pipeline for the PfPan project, covering the construction of a *P. falciparum* pangenome graph from long-read assemblies, characterisation of genetic variation, population-level short-read analysis, and benchmarking of pangenome versus linear reference variant calling approaches.

---

## Pipeline Overview

```
1. Pangenome Construction
   └── Minigraph-Cactus graph from Pf3k PacBio assemblies → GBZ/GFA/VCF outputs

2. Variation Analysis
   └── VCF preprocessing → functional annotation → SV characterisation → PCA/network

3. Short-Read Population Analysis
   └── Pangenome + linear SNP calling → PCoA → Jost's D → CNV → SV pipeline

4. Mapping & Calling Benchmarking
   └── Assembly truth sets → vg giraffe + GATK → Truvari + vcfeval evaluation
```

---

## Documentation

| Module | Description |
|--------|-------------|
| [Pangenome Construction](README_pangenome_construction.md) | Automated construction of the PfPan pangenome graph using Minigraph-Cactus, including assembly download, seq file generation, graph building, and pangenome statistics (panacus, vg stats) |
| [Variation Analysis](README_pfpan_variation_analysis.md) | Variant extraction and annotation from the Cactus VCF, genomic characterisation, SV analysis, gene-of-interest extraction, PCA, and haplotype network |
| [Short-Read Population Analysis](README_short_read_population_analysis.md) | Population-level analysis of short-read data aligned to the pangenome graph and linear reference, including PCoA comparison, Jost's D differentiation, CNV calling, and SV characterisation |
| [Mapping & Calling Benchmarking](README_mapping_calling_benchmarking.md) | Benchmarking of pangenome (vg giraffe + vg call) versus linear reference (BWA + GATK + DELLY) approaches for SNPs, indels, and SVs against assembly-based truth sets |

---

## Repository Structure

```
Pfpan/
├── README.md                                   ← this file
├── README_pangenome_construction.md
├── README_pfpan_variation_analysis.md
├── README_short_read_population_analysis.md
├── README_mapping_calling_benchmarking.md
│
├── pipeline/
│   ├── construct_PfPan.py                      ← pangenome construction (Step 1)
│   ├── wget_pf3k_fasta.py
│   ├── make_seqtxt.py
│   ├── replot_panacus.py
│   ├── PfPan_SVIM_ASM.py                       ← assembly truth set (Step 2)
│   ├── PfPan_linear_map_and_call.py            ← linear mapping pipeline (Step 3a)
│   ├── PfPan_map_and_call.py                   ← pangenome mapping pipeline (Step 3b)
│   └── PfPan_variant_evaluation.py             ← benchmarking (Step 4)
│
├── analysis/
│   ├── PfPan_variation.R                       ← variation analysis
│   ├── PfPan_pca_network.py                    ← PCA and haplotype network
│   ├── PfPan_pca_plink.R                       ← PCoA comparison
│   ├── PfPan_josts_d.py                        ← Jost's D differentiation
│   ├── PfPan_cnv_call.py                       ← CNV detection
│   ├── PfPan_sv_pipeline.sh                    ← SV characterisation pipeline
│   ├── PfPan_extract_sv_sequences.py           ← standalone SV sequence extraction
│   └── PfPan_annotate_svs.py                   ← standalone SV annotation
│
├── PfPan_Pf3D7_pan/                            ← dir contaning PfPan files (can be found here: https://genomics.lshtm.ac.uk/data/)

```

---

## Quick Start

### 1. Build the pangenome graph

```bash
conda activate cactus
python ./pipeline/construct_PfPan.py --download   # includes assembly download
```

See [Pangenome Construction](README_pangenome_construction.md) for full usage and step-skipping options.

### 2. Characterise variation

Run Stage 1 (VCF preprocessing) from [Variation Analysis](README_pfpan_variation_analysis.md), then the R and Python analysis scripts.

### 3. Generate the assembly-based truth set (once)

```bash
python ./pipeline/PfPan_SVIM_ASM.py
```

### 4. Run variant calling pipelines

```bash
# Linear reference
conda activate <fastq2matrix_env>
python ./pipeline/PfPan_linear_map_and_call.py --samples-file fastqs.txt --skip-delly
conda activate delly
python ./pipeline/PfPan_linear_map_and_call.py --samples-file fastqs.txt --skip-fastq2vcf

# Pangenome
conda activate cactus
python ./pipeline/PfPan_map_and_call.py --sample <SAMPLE>
```

### 5. Benchmark and evaluate

```bash
python ./pipeline/PfPan_variant_evaluation.py
```

See [Mapping & Calling Benchmarking](README_mapping_calling_benchmarking.md) for individual step options.

---

## Dependencies

| Tool | Used in |
|------|---------|
| Minigraph-Cactus (`cactus-pangenome`) | Pangenome construction |
| `vg` (giraffe, surject, call, pack, stats) | Graph mapping and variant calling |
| `panacus` | Pangenome growth statistics |
| `odgi` | Graph region extraction and visualisation |
| `bcftools` | VCF normalisation, filtering, annotation, querying |
| `truvari` | SV annotation and benchmarking |
| `GATK` | Short variant calling |
| `DELLY` | Linear-mapping SV calling |
| `minimap2` + `svim-asm` | Assembly-based truth set generation |
| `paftools.js` | Assembly-based SNP/indel calling |
| `RTG Tools` (vcfeval) | Short variant benchmarking |
| `MAFFT` + `BLAST` | SV sequence extraction and annotation |
| `mosdepth` | Coverage profiling for CNV |
| `plink` | SNP distance matrix generation |
| R (≥4.0) | Statistical analysis and plotting |
| Python (≥3.8) | PCA, network, differentiation, CNV, SV pipelines |

---

## Reference

All analyses use the *Plasmodium falciparum* 3D7 reference genome version 3 (`Pfalciparum.genome.fasta`) and the Pf3k long-read assembly dataset. Pf3D7 (West African origin) is used as the pangenome reference backbone for variant projection and interpretation.
