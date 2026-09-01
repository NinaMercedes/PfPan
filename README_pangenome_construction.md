# PfPan Pangenome Construction Pipeline

This pipeline automates the construction of a *Plasmodium falciparum* pangenome graph from Pf3k long-read assemblies using Minigraph-Cactus, and computes downstream pangenome statistics using panacus and vg. The pipeline is orchestrated by a single Python script (`construct_PfPan.py`) which manages environment activation, working directories, and all intermediate steps automatically.

---

## Requirements

- **Conda environment:** `cactus` (Python 3.8) — used to run the pipeline script itself
- **Cactus virtualenv:** sourced automatically at runtime for `cactus-pangenome`
- **Tools on PATH:** `panacus`, `vg`, `wget`
- **Pipeline scripts** (all located in the same directory as `construct_PfPan.py`):
  - `wget_pf3k_fasta.py` — downloads assemblies from Sanger FTP
  - `make_seqtxt.py` — generates the seq file for cactus
  - `plot_figure1_node.py` — plots panacus growth curves
  - `plot_figure1_bp.py` — plots panacus growth curves

---

## Directory Structure

All paths are resolved relative to the location of `construct_PfPan.py` using `Path(__file__).resolve()`, ensuring the pipeline works correctly regardless of which directory it is invoked from.

```
custom_scripts/
├── pipeline/
│   ├── construct_PfPan.py       ← main pipeline script
│   ├── wget_pf3k_fasta.py
│   ├── make_seqtxt.py
│   ├── mplot_figure1_bp.py
│   └── plot_figure1_node.py
├── genomes/                     ← assembly FASTA files and seq file (auto-created)
├── PfPan_Pf3D7_pan/             ← pangenome graph outputs (auto-created)
├── PfPan_Pf3D7_js/              ← Cactus job store, can be deleted after run (auto-created)
├── stats_panacus/               ← pangenome statistics outputs (auto-created)
└── pipeline.log                 ← timestamped run log
```

---

## Usage

```bash
conda activate cactus

# Full run including assembly download
python ./pipeline/construct_PfPan.py --download

# Note post-download replaced the 3D7 sequence with the standard Plasmodium falciparum 3D7 reference version 3

# Full run, assemblies already present in genomes/
python ./pipeline/construct_PfPan.py

# Run individual steps selectively
python ./pipeline/construct_PfPan.py --skip-seq            # skip seq file generation
python ./pipeline/construct_PfPan.py --skip-graph          # skip graph construction
python ./pipeline/construct_PfPan.py --skip-stats          # skip statistics
python ./pipeline/construct_PfPan.py --skip-seq --skip-graph   # statistics only
```

---

## Pipeline Steps

### Optional: Download Assemblies (`--download`)

**Script:** `wget_pf3k_fasta.py`

Fetches *P. falciparum* PacBio long-read assemblies from the Pf3k Sanger FTP server. Rather than using recursive wget (which only retrieves the index page), the script explicitly parses the FTP directory listing for `.fasta.gz` hrefs and downloads each file individually, skipping `.fai` index files. Files are downloaded directly into `genomes/` by running wget with `cwd=GENOME_DIR`.

**Source:** `https://ftp.sanger.ac.uk/pub/project/pathogens/Plasmodium/falciparum/PF3K/ReferenceGenomes_Version1/GENOMES/`

**Output — `genomes/`:**

| File | Size | Description |
|------|------|-------------|
| `Pfalciparum.genome.fasta.gz` | 6.1M | 3D7 Reference strain version 3 — West African origin |
| `Pf7G8.April2018.fasta.gz` | 5.9M | Lab strain |
| `PfCD01.April2018.fasta.gz` | 6.1M | Clinical isolate — Cambodia |
| `PfDd2.April2018.fasta.gz` | 5.9M | Lab strain — Indochina origin |
| `PfGA01.April2018.fasta.gz` | 6.0M | Clinical isolate — Gabon |
| `PfGB4.April2018.fasta.gz` | 6.1M | Lab strain |
| `PfGN01.April2018.fasta.gz` | 6.2M | Clinical isolate — Guinea |
| `PfHB3.April2018.fasta.gz` | 5.9M | Lab strain |
| `PfIT.April2018.fasta.gz` | 6.0M | Lab strain — Italian origin |
| `PfKE01.April2018.fasta.gz` | 5.9M | Clinical isolate — Kenya |
| `PfKH01.April2018.fasta.gz` | 6.1M | Clinical isolate — Cambodia |
| `PfKH02.April2018.fasta.gz` | 6.0M | Clinical isolate — Cambodia |
| `PfML01.April2018.fasta.gz` | 6.7M | Clinical isolate (not included) — Mali |
| `PfSD01.April2018.fasta.gz` | 5.9M | Clinical isolate (not included) — Sudan |
| `PfSN01.April2018.fasta.gz` | 6.1M | Clinical isolate — Senegal |
| `PfTG01.April2018.fasta.gz` | 6.7M | Clinical isolate (not included) — Togo |


---

### Step 1: Generate Seq File

**Script:** `make_seqtxt.py`

Generates the two-column sequence file required by `cactus-pangenome`. The script lists all `.fasta.gz` files in `genomes/`, extracts the sample name from the filename (everything before the first `.`), and writes the sample name and its absolute path as space-separated columns. The script runs with `cwd=GENOME_DIR` to ensure it finds the assemblies and writes the seq file to the correct location.

**Output — `genomes/`:**

| File | Description |
|------|-------------|
| `pf3k_seq_v2.txt` | Two-column seq file mapping sample names to absolute FASTA paths |
| `fasta_list.txt` | Intermediate file listing all `.fasta.gz` filenames — can be ignored |

Example `pf3k_seq_v2.txt` content:
```
Pf3D7 /mnt/storage13/nbillows/custom_scripts/genomes/Pf3D7.April2018.fasta.gz
PfDd2 /mnt/storage13/nbillows/custom_scripts/genomes/PfDd2.April2018.fasta.gz
...
```

---

### Step 2: Build Pangenome Graph

**Tool:** `cactus-pangenome` (run inside Cactus virtualenv, sourced automatically)

Constructs the pangenome graph using Minigraph-Cactus. Pf3D7 (West African origin) is used as the reference backbone sequence. Although the pangenome is largely reference-free in structure, a reference is required for variant projection and interpretation. The Cactus virtualenv is sourced and deactivated automatically within a single bash session, so no manual environment management is needed.

**Cactus flags used:**

| Flag | Description |
|------|-------------|
| `--reference Pf3D7` | Reference sequence for variant projection and graph backbone |
| `--filter 2` | Filter out alignments supported by fewer than 2 haplotypes |
| `--haplo` | Generate haplotype index for personalised graph construction with `vg haplotypes` |
| `--giraffe clip filter` | Build Giraffe short-read mapping indices on clipped and filtered graphs |
| `--viz` | Generate visual representation of the graph |
| `--odgi` | Output full graph in ODGI format |
| `--chrom-vg clip filter` | Per-chromosome vg graphs at clipped and filtered levels |
| `--chrom-og` | Per-chromosome ODGI graphs |
| `--gbz clip filter full` | GBZ format graphs at clipped, filtered, and full levels |
| `--gfa clip full` | GFA format graphs at clipped and full levels |
| `--vcf --vcfReference Pf3D7` | Call and output variants projected onto Pf3D7 |
| `--consCores 8` | Number of cores for the consensus step |
| `--mgMemory 128Gi` | Memory allocated to Minigraph |

**Output — `PfPan_Pf3D7_pan/`:**

| File | Size | Format | Description |
|------|------|--------|-------------|
| `PfPan.gbz` | 42M | GBZ | Clipped pangenome graph — primary file used for most downstream tools |
| `PfPan.full.gbz` | 52M | GBZ | Full (unclipped) pangenome graph |
| `PfPan.d2.gbz` | 38M | GBZ | Distance-indexed graph for Giraffe short-read mapping |
| `PfPan.gfa.gz` | 37M | GFA | Clipped pangenome graph in GFA format |
| `PfPan.full.gfa.gz` | 44M | GFA | Full pangenome graph in GFA format |
| `PfPan.sv.gfa.gz` | 6.9M | GFA | Structural variant subgraph |
| `PfPan.sv.gfa.fa.gz` | 6.9M | FASTA | Sequences associated with the SV subgraph |
| `PfPan.full.hal` | 152M | HAL | Full whole-genome alignment in Cactus native format |
| `PfPan.full.og` | 392M | OG | Full graph in ODGI format |
| `PfPan.full.snarls` | 3.1M | Snarls | Snarl decomposition of the full graph |
| `PfPan.snarls` | 3.0M | Snarls | Snarl decomposition of the clipped graph |
| `PfPan.dist` | 100M | Dist | Distance index for the clipped graph |
| `PfPan.d2.dist` | 40M | Dist | Distance index for the Giraffe-ready graph |
| `PfPan.d2.min` | 290M | Min | Minimiser index for Giraffe read mapping |
| `PfPan.d2.snarls` | 3.6M | Snarls | Snarl index for the Giraffe-ready graph |
| `PfPan.min` | 325M | Min | Minimiser index for the clipped graph |
| `PfPan.hapl` | 21M | Hapl | Haplotype index for personalised graph construction |
| `PfPan.gaf.gz` | 17M | GAF | Read-to-graph alignments in graph alignment format |
| `PfPan.paf` | 33M | PAF | Pairwise alignments of input assemblies to the graph |
| `PfPan.paf.unfiltered.gz` | 7.9M | PAF | Unfiltered pairwise alignments |
| `PfPan.paf.filter.log` | 34K | Log | Log of filtering applied to PAF alignments |
| `PfPan.vcf.gz` | 12M | VCF | Filtered variant calls projected onto Pf3D7 |
| `PfPan.vcf.gz.tbi` | 12K | TBI | Tabix index for filtered VCF |
| `PfPan.raw.vcf.gz` | 15M | VCF | Raw (unfiltered) variant calls |
| `PfPan.raw.vcf.gz.tbi` | 13K | TBI | Tabix index for raw VCF |
| `PfPan.stats.tgz` | 67K | TGZ | Archive of graph statistics |
| `PfPan.viz/` | — | Dir | Visual representations of the pangenome graph |
| `PfPan.chroms/` | — | Dir | Per-chromosome graph files |
| `chrom-alignments/` | — | Dir | Per-chromosome alignment files |
| `chrom-subproblems/` | — | Dir | Intermediate files from per-chromosome Cactus alignment — can be deleted after a successful run |
| `pf3k_seq_v2.txt` | 1.5K | TXT | Copy of the input seq file |

---

### Step 3: Pangenome Statistics

Two complementary tools are used to characterise the pangenome: **panacus** for growth curve analysis and **vg stats** for basic graph metrics. The GFA graph is copied and decompressed into `stats_panacus/` for use by panacus. If the haplotypes file (`paths.haplotypes.txt`) does not already exist, it is generated automatically from the first column of the seq file.

#### 3a. Panacus — Base Pair Coverage

Runs `panacus hist` to compute a histogram of base pair coverage across haplotypes, then `panacus histgrowth` to model how much sequence is core, soft-core, shell, or private as haplotypes are incrementally added to the pangenome. Results are plotted as a growth curve PDF.

Growth curve coverage thresholds (`-l 1,2,1,1,1 -q 0,0,1,0.5,0.1`):

| Category | Threshold | Description |
|----------|-----------|-------------|
| Core | 100% | Present in all haplotypes |
| Soft-core | ≥50% | Present in at least half of haplotypes |
| Shell | ≥10% | Present in a subset of haplotypes |
| Private | unique | Found in only one haplotype |

**Output — `stats_panacus/bp/`:**

| File | Size | Description |
|------|------|-------------|
| `bp.hist` | 2.0K | Histogram of base pair coverage counts across haplotypes |
| `bp.growth` | 1.2K | Pangenome growth table for base pairs — core/soft-core/shell/private breakdown |
| `panacus_growth_bp.pdf` | 20K | Growth curve plot for base pairs |

#### 3b. Panacus — Node Coverage

Runs `panacus histgrowth` at the graph node (segment) level, using the same coverage thresholds as above, and additionally computes pairwise node-level similarity between all haplotypes.

**Output — `stats_panacus/node/`:**

| File | Size | Description |
|------|------|-------------|
| `node.growth` | 1.1K | Pangenome growth table at node level |
| `panacus_growth_node.pdf` | 20K | Growth curve plot for nodes |

#### 3c. vg Stats — Basic Graph Metrics

Runs `vg stats` on the clipped GBZ graph to report basic graph-level summary statistics including total sequence length, number of nodes, number of edges, and snarl counts.

**Output — `stats_panacus/`:**

| File | Size | Description |
|------|------|-------------|
| `PfPan.gfa` | 205M | Decompressed GFA used as input for panacus |
| `paths.haplotypes.txt` | 105B | List of haplotype names — auto-generated from seq file |
| `vg_stats.txt` | 296B | Basic graph metrics from `vg stats` |

---

## Configuration

Key parameters are set at the top of `construct_PfPan.py` and can be adjusted as needed:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `CACTUS_VENV` | `/mnt/storage13/.../bin/activate` | Path to the Cactus virtualenv activate script |
| `REFERENCE` | `Pf3D7` | Reference strain for variant projection |
| `GRAPH_NAME` | `PfPan` | Prefix for all output files |
| `CONS_CORES` | `8` | CPU cores for the consensus step |
| `MG_MEMORY` | `128Gi` | Memory for Minigraph |
| `HIST_LEVELS` | `1,2,1,1,1` | Panacus coverage level thresholds |
| `HIST_QUANTILES` | `0,0,1,0.5,0.1` | Panacus core/shell/private quantile thresholds |

