# PfPan leave-one-out (LOO) benchmarking: pipeline, Truvari sensitivity sweep and plot

Three scripts run in order for the benchmarking. Two upstream mapping-and-calling scripts (`PfPan_map_and_call.py`, `PfPan_linear_map_and_call.py`) make the read-derived inputs and are described under "Upstream scripts" below.

| Order | Script | Purpose |
|---|---|---|
| 1 | `full_pipeline.py` | Builds each leave-one-out graph, maps reads, calls variants with every caller, and benchmarks against the assembly-derived truth sets (vcfeval for small variants, Truvari for SVs). |
| 2 | `truvari_param_sweep.py` | Re-runs only the Truvari SV benchmark under alternative matching parameters, using the VCFs already built by step 1. |
| 3 | `plot_truvari_sweep.py` | Plots median precision, recall and F1 per caller per parameter setting (Supplementary Figure 16). |

`truvari_param_sweep.py` imports `full_pipeline.py` (`import full_pipeline as fp`), so both must be in the same directory. Run everything from that directory.

---

## Requirements and software versions

- Python 3 with `pandas` and `matplotlib` (plot script only).
- All external tools are called as shell commands (`bash -c`, `set -euo pipefail`). The Python code handles orchestration only. Each command activates a named conda environment, so the environment names must match those below.
- Exact tool versions are recorded in the environment YAML files in the `envs/` folder of the GitHub pipeline directory.

| Environment (name used in the script) | Used for | Version file |
|---|---|---|
| `cactus` | `vg` (giraffe, surject, pack, call, paths), samtools, bcftools, bgzip | `envs/cactus.yml` |
| `fastq2matrix` | GATK, bcftools, samtools, bwa, tabix | `envs/fastq2matrix.yml` |
| `truvari` | truvari, bedtools, rtg-tools (vcfeval, vcfdecompose) | `envs/truvari.yml` |
| `delly` | Delly | `envs/delly.yml` (to add) |
| `manta` | Manta | `envs/manta.yml` (to add) |
| `dysgu` | Dysgu | `envs/dysgu.yml` (to add) |

Tool versions used (from the Methods; the YAMLs should agree with these):

| Tool | Version | Tool | Version |
|---|---|---|---|
| fastp | v1.0.1 | Delly | v1.5.0 |
| Hostile | v2.0.0 | Manta | v1.6.0 |
| vg (giraffe, pack, call, surject) | v1.61.0 | Dysgu | v1.9.0 |
| GATK | v4.1.4.1 | bwa | v0.7.17 |
| bcftools | v1.12 | Truvari | v5.3.0  |
| Cactus-pangenome | v2.9.3 | | |

Not covered by a conda YAML, so record these separately (for example in the same folder or in Methods):

- **Cactus-pangenome v2.9.3**: a Python virtualenv (`--cactus-venv`, `cactus-bin-v2.9.3/venv-cactus-v2.9.3`), not a conda environment.
- **`vt`**: a standalone binary (`--vt-bin`), used only for the small-variant repeat-tract masks. Record the commit or version.
- **paftools / minimap2 and SVIM-asm**: used to make the truth sets before this pipeline. Record their versions too.

To export each environment (run with the environment active, or use `-n`):

```bash
mkdir -p envs
for e in cactus fastq2matrix truvari delly manta dysgu; do
    conda env export -n "$e" --no-builds > envs/"$e".yml
done
```

`--no-builds` keeps the YAML portable across machines. Drop it if exact build strings are needed to reproduce the environment. Recreate an environment with `conda env create -f envs/<name>.yml`.

---

## How the input files were generated

`full_pipeline.py` does not make its truth sets or region files. They are built beforehand. The commands below are the standard way to make files.

| Input file | Made with | Status |
|---|---|---|
| Assemblies (`pf3k_seq_v2.txt` paths) | Near-complete long-read assemblies, one per strain, plus Pf3D7 | 
| `{strain}.paftools.snps_indels.vcf.gz` (small-variant truth) | minimap2 + `paftools.js call`, strain assembly vs Pf3D7 | 
| `{strain}.svim.svs.vcf.gz` (SV truth) | minimap2 + SVIM-asm, strain assembly vs Pf3D7 | 
| `{strain}.pan.snps_indels.vcf.gz`, `{strain}.pan.svs.vcf.gz` (pan_direct) | Variants for that strain taken directly from the pangenome VCF, with no reads involved - using VG Giraffe| 
| `{strain}.April2018.confident.bed` | Confident regions for each strain's truth set |
| `Core_genome_Pf3D7_v3_ext.bed` | Core-genome regions of Pf3D7 v3 | 
| `Pfalciparum_TRF_repeats_sensitive.merged.bed` | Tandem Repeats Finder run on Pf3D7, hits merged | 
| `Pfalciparum_variable_repetitive_regions.bed` | Built from the GFF by `build_shared_beds()` (var, rif, SURF, etc. genes, repeat regions and centromeres, padded 500 bp and merged) | 
| Category BEDs in `--beddir` | Genome-category regions used by vcfeval | 
| Reads (`{strain}_1/2.trimmed.fastq.gz`) | Raw Illumina reads, host reads removed with Hostile v2.0.0 and trimmed/filtered with fastp v1.0.1 | 
| `{strain}.bqsr.cram`, `{strain}.g.vcf.gz` (linear baseline) | `PfPan_linear_map_and_call.py` (fastq2vcf.py: bwa v0.7.17, BQSR, GATK HaplotypeCaller `-ERC GVCF`) | 
| `Pfalciparum.genome.fasta` and its `.sdf` | Pf3D7 v3 reference, then `rtg format` | 

### Upstream scripts: read processing, mapping and calling

**Read processing.** Raw Illumina reads (strains from public data) were filtered before any mapping:

```bash
# host (human) read removal
hostile clean --fastq1 STRAIN_1.fastq.gz --fastq2 STRAIN_2.fastq.gz ...     # Hostile v2.0.0
# adapter/quality trimming and filtering
fastp -i STRAIN_1.clean.fastq.gz -I STRAIN_2.clean.fastq.gz \
      -o STRAIN_1.trimmed.fastq.gz -O STRAIN_2.trimmed.fastq.gz ...         # fastp v1.0.1 
```

The trimmed files `{strain}_1.trimmed.fastq.gz` and `{strain}_2.trimmed.fastq.gz` are the read inputs to every script below. 

**`PfPan_map_and_call.py` (graph mapping and calling).** Run per sample with `python PfPan_map_and_call.py --sample <name>`. Reads in `./reads/`, outputs in `./output/`. Inputs are the pangenome `.gbz` and `.snarls`, a list of the Pf3D7 paths (`vg paths -x <gbz> -L | grep Pf3D7 > Pf3D7.paths.txt`) and the Pf3D7 FASTA.

| Step | Command | Output |
|---|---|---|
| 1 | `vg giraffe -Z <gbz> -f R1 -f R2 -o gaf` | `<sample>.gaf.gz` |
| 2 | `vg surject -x <gbz> -G <gaf> --interleaved -F Pf3D7.paths.txt -b` then `samtools reheader` to strip `Pf3D7#0#` | `<sample>.bam` |
| 3 | `samtools sort`, `index`, `flagstat` | `<sample>_sort.bam`, `.bai`, `_stat.txt` |
| 4 | `gatk HaplotypeCaller -ERC GVCF` | `<sample>_gatk.g.vcf.gz` |
| 5 | `delly call` (skip with `--skip-delly`) | `<sample>_delly_sites.bcf` |
| 6 | `vg pack -Q5 -a <gaf>` | `<sample>.pack` |
| 7 | `vg call ... --ploidy 1 -S Pf3D7 -az` (haploid) | `<sample>.SV.call_haploid_.vcf.gz` |

**`PfPan_linear_map_and_call.py` (linear baseline).** Run with `python PfPan_linear_map_and_call.py --samples-file fastqs.txt` (one sample name per line). Reads in `../reads/`, outputs in `../linear/`.

| Step | What | Output |
|---|---|---|
| 1 | `fastq2vcf.py all --mapper bwa` per sample: bwa mapping, BQSR, GATK HaplotypeCaller (`fastq2matrix` environment) | `<sample>.bqsr.bam`, gVCF |
| 2 | `delly call` per sample (`delly` environment) | `<sample>_delly.bcf` |
| 3 | `delly merge` across samples | `sites.bcf` |
| 4 | `delly call -v sites.bcf` per sample (forced genotyping) | `<sample>_delly_sites.bcf` |
| 5 | `bcftools merge -m id` | `pan_delly.vcf.gz` |

BQSR uses a high-quality *P. falciparum* genetic-crosses dataset as known sites. VQSR is deliberately not applied. This script is for the multi-sample population analysis. The LOO benchmark uses only its per-sample BAM and gVCF, and `full_pipeline.py` runs Delly, Manta and Dysgu itself per sample, so the merged Delly steps (2 to 5) are not part of the benchmark.

**Filtering:**

- The population Methods filter GATK calls with QUAL < 20, MQ < 40, QD < 2, FS > 60, rank-sum outside +/-12.5 and AF == 0. The benchmark in `full_pipeline.py` (`filter_gatk`) uses different thresholds (SNPs: QD < 2, QUAL < 30, SOR > 4, FS > 60, MQ < 40, MQRankSum < -15, ReadPosRankSum < -5; indels: QD < 2, FS > 200, ReadPosRankSum < -20).


### Small-variant truth (paftools)

```bash
minimap2 -cx asm5 --cs -t 16 Pf3D7.fa STRAIN.fa \
  | sort -k6,6 -k8,8n \
  | paftools.js call -f Pf3D7.fa -s STRAIN - \
  | bgzip > STRAIN.paftools.snps_indels.vcf.gz
tabix -p vcf STRAIN.paftools.snps_indels.vcf.gz
```

`reheader_truth()` in the pipeline later renames the VCF sample to the strain name.

### SV truth (SVIM-asm)

```bash
minimap2 -a -x asm5 --cs -r2k -t 16 Pf3D7.fa STRAIN.fa \
  | samtools sort -@ 4 -m 4G -o STRAIN.asm.bam
samtools index STRAIN.asm.bam
svim-asm haploid svim_STRAIN/ STRAIN.asm.bam Pf3D7.fa
# keep SVs only (>= 50 bp), then compress and index
bcftools view -i 'abs(INFO/SVLEN)>=50' svim_STRAIN/variants.vcf -Oz -o STRAIN.svim.svs.vcf.gz
tabix -p vcf STRAIN.svim.svs.vcf.gz
```

### pan_direct

These are the strain's own variants read directly out of the full pangenome (Cactus `--vcf`, reference Pf3D7), so no reads are used. A typical approach is to take that strain's column, drop reference calls, and split into small variants (< 50 bp) and SVs (>= 50 bp). Replace this with the exact commands used:

```bash
bcftools view -s STRAIN PfPan.vcf.gz | bcftools norm -f Pf3D7.fa -m -both ... 
# small:  abs(strlen(ALT)-strlen(REF)) < 50   -> STRAIN.pan.snps_indels.vcf.gz
# SV:     abs(strlen(ALT)-strlen(REF)) >= 50  -> STRAIN.pan.svs.vcf.gz
```

### Reference and RTG SDF

```bash
samtools faidx Pfalciparum.genome.fasta
bwa index Pfalciparum.genome.fasta
rtg format -o Pfalciparum.genome.fasta.sdf Pfalciparum.genome.fasta
```

### Tandem-repeat BED 

```bash
trf Pfalciparum.genome.fasta 2 7 7 80 10 50 500 -f -d -h -m
# convert the .dat output to BED, then
sort -k1,1 -k2,2n trf.bed | bedtools merge -i - > Pfalciparum_TRF_repeats_sensitive.merged.bed
```

### BED files
The BED files used in the analysis can be found in the pipeline directory.

---

## 1. `full_pipeline.py`

### What it does, per strain

Each strain is held out of the graph, then re-genotyped from its own reads. Truth is the strain's own assembly.

| Stage (toggle) | What happens |
|---|---|
| `run_loo` (`--run-mapping`) | Removes the strain from the seq file, builds `PfPan_no_<strain>` with `cactus-pangenome` (`--filter 2 --haplo --giraffe clip filter --gbz clip filter full --gfa clip full --vcf`, reference Pf3D7), maps reads with `vg giraffe` (GAF), surjects to Pf3D7 (BAM), runs `vg pack -Q5` and `vg call --ploidy 1`, and runs GATK HaplotypeCaller on the surjected BAM. |
| `run_gatk_linear`, `run_gatk_surject_loo`, `run_vg_hap_loo_filter`, `build_confidence`, `run_vcfeval_all` (`--run-small-var`) | Filters small variants (< 50 bp) from linear GATK, surject+GATK and `vg call`. Builds confidence-region BEDs and benchmarks with `rtg vcfeval` against the paftools truth. |
| `run_delly_*` (`--run-sv`) | Delly on the LOO surjected BAM and on the linear BAM. |
| `run_manta_*` (`--run-sv-extra`) | Manta on the LOO surjected BAM and on the linear BAM. |
| `run_dysgu_*` (`--run-sv-dysgu`, default on) | Dysgu on the LOO surjected BAM and on a matched BWA-MEM BAM made from the same trimmed FASTQs, so mapping approach is the only difference. |
| `run_truvari_all` | Benchmarks every enabled SV caller against the SVIM-asm truth in the core genome only. |

The stage toggles are all off by default except Dysgu. Pass `--run-mapping`, `--run-small-var`, `--run-sv`, `--run-sv-extra` or `--run-sv-dysgu` to turn stages on, and `--no-<name>` to turn them off.

### Inputs

| Argument | What it is | Expected content |
|---|---|---|
| `--strains` | Strains to process | Default: the 12 held-out strains (PfDd2, Pf7G8, PfCD01, PfGA01, PfGB4, PfGN01, PfHB3, PfIT, PfKE01, PfKH01, PfKH02, PfSN01). Together with the Pf3D7 reference these make up the 13 pangenome samples. |
| `--seq-file` | Cactus seq file (`pf3k_seq_v2.txt`) | Sample name and assembly path per line, including Pf3D7. The strain is removed with `grep -v`. |
| `--reads-dir` | Illumina reads and linear alignments | `{strain}_1.trimmed.fastq.gz`, `{strain}_2.trimmed.fastq.gz`, `{strain}.bqsr.cram`, and raw linear GATK gVCF `{strain}.g.vcf.gz`. |
| `--linear-gvcf-dir` | Linear gVCFs | Defaults to `--reads-dir`. |
| `--loo-base` | Where `loo_<strain>/` working directories are created | Each holds the LOO graph, BAM, GAF, pack and VCFs. |
| `--comparison-dir` | Pre-existing truth and direct-from-assembly calls | `{strain}.paftools.snps_indels.vcf.gz` (small-variant truth), `{strain}.svim.svs.vcf.gz` (SV truth), `{strain}.pan.snps_indels.vcf.gz` and `{strain}.pan.svs.vcf.gz` (pan_direct calls). |
| `--ref`, `--ref-gatk`, `--ref-sdf` | Pf3D7 v3 FASTA, the same FASTA at the path used for GATK, and the RTG SDF | `--ref-gatk` and `--ref` should be the same sequence. |
| `--gff` | Pf3D7 annotation | Used to build the variable-surface-antigen BED. |
| `--beddir` | Category BEDs used by vcfeval | One BED per genome category. |
| `--vsa-bed`, `--trf-bed` | VSA and tandem-repeat masks | Built from the GFF if the VSA BED is missing. The TRF BED must already exist. |
| `--conf-bed-dir` | Per-strain confident regions | `{strain}.April2018/{strain}.April2018.confident.bed`. |
| `--core-bed` | Core-genome BED | **Region used for all SV benchmarking** (`Core_genome_Pf3D7_v3_ext.bed`). |
| `--outroot` | Output root | See outputs below. |
| `--vsa-pad`, `--vt-min-fz-rl`, `--boundary-buffer`, `--density-window`, `--density-percentile`, `--min-bp` | Confidence-region parameters | Defaults 500, 4, 1000, 200, 95, 1000. |

### Outputs (under `--outroot/<strain>/` unless stated)

| Path | Content |
|---|---|
| `<loo-base>/loo_<strain>/` | LOO graph (`PfPan_no_<strain>_pan/`: `.gbz`, `.snarls`, `.gfa`), GAF, `<strain>.loo.sort.bam`, `.pack`, vg call VCF, GATK gVCF, and Delly/Manta/Dysgu VCFs. |
| `gatk_linear/`, `gatk_surject_loo/` | Filtered, decomposed small-variant VCFs. |
| `beds/` | Confidence-region BEDs (`*.FINAL.bed`). |
| `vcfeval/<caller>/<category>/` | vcfeval results for each caller and genome category. |
| `delly_linear/`, `manta_linear/`, `dysgu_linear/`, `matched_linear/` | Linear-baseline calls and BAMs. |
| `<strain>.svim.svs.alt_only.mixed.vcf.gz` | SV truth with `GT!="ref"` kept. |
| `*.alt_only.mixed.vcf.gz` next to each caller's SV VCF | **Filtered caller VCFs that the sweep re-uses.** |
| `truvari/<caller>/<core bed name>/` | Truvari bench output (including `summary.json`) at the chosen settings. |

### SV benchmarking details (`run_truvari_all`)

Chosen Truvari command:

```
truvari bench -b <truth> -c <calls> -o <out> -f <ref> \
    -r 1000 -C 1000 -O 0.0 -p 0.0 -P 0.3 -s 50 -S 15 --sizemax 10000 \
    --includebed <core.bed>
```

Genotype filtering before Truvari, to keep only calls the caller considers present:

| Caller | Filter |
|---|---|
| truth (SVIM-asm), pan_direct_sv, vg_hap_sv_loo, Delly | `GT!="ref"` |
| Manta (LOO and linear) | `GT="1/1"` (Manta emits no AF field here) |
| Dysgu (LOO and linear) | `GT="1/1" \|\| FORMAT/AF>=0.7` |

The Manta and Dysgu filters are applied because they report diploid-style genotypes. `vg call` was run with `--ploidy 1`, so it is already haploid.

---

## 2. `truvari_param_sweep.py`

Responds to the reviewer request to justify the Truvari thresholds and show robustness to alternatives. It re-runs `truvari bench` for each caller and strain under each parameter setting. **It does not rebuild or re-filter any VCFs.** It reads the `*.alt_only.mixed.vcf.gz` files and the SV truth already made by `run_truvari_all`, so only the matching parameters change.

### Prerequisite

Run `full_pipeline.py` SV stages first (with the relevant toggles) so that these exist:

- Truth: `<outroot>/<strain>/<strain>.svim.svs.alt_only.mixed.vcf.gz`
- Calls, per caller (`caller_alt_vcf()`):

| Caller | Source VCF (before `.alt_only.mixed`) |
|---|---|
| `pan_direct_sv` | `<comparison-dir>/<strain>.pan.svs.vcf.gz` |
| `delly_loo`, `manta_loo`, `dysgu_loo` | `<loo-base>/loo_<strain>/<strain>.loo.<caller>.SV.ge50bp.PASS.vcf.gz` |
| `vg_hap_sv_loo` | `<loo-base>/loo_<strain>/<strain>.loo.pan.SV.ge50bp.PASS.hap.vcf.gz` |
| `delly_linear`, `manta_linear`, `dysgu_linear` | `<outroot>/<strain>/<caller>/<strain>.<caller>.SV.ge50bp.PASS.vcf.gz` |

Missing files are reported and skipped (the row is written with an `error` value).

### Settings tested

All other options are fixed: `-s 50 -S 15 --sizemax 10000 --includebed <core.bed>`. `--chunksize` is set to `max(1000, refdist)` because Truvari requires it to be at least `refdist`.

| Setting | refdist | pctseq | pctsize | pctovl |
|---|---|---|---|---|
| `current` (**chosen**) | 1000 | 0.0 | 0.3 | 0.0 |
| `truvari_default` | 500 | 0.7 | 0.7 | 0.0 |
| `refdist_200` / `_500` / `_2000` | 200 / 500 / 2000 | 0.0 | 0.3 | 0.0 |
| `pctseq_0.3` / `_0.7` | 1000 | 0.3 / 0.7 | 0.3 | 0.0 |
| `pctsize_0.5` / `_0.7` | 1000 | 0.0 | 0.5 / 0.7 | 0.0 |
| `pctovl_0.3` | 1000 | 0.0 | 0.3 | 0.3 |
| `moderate_all` | 500 | 0.3 | 0.5 | 0.0 |

Truvari's own default for reciprocal overlap (`pctovl`) is 0, so "no reciprocal overlap" is not a loosening relative to the default.

### Usage

```bash
# default: all 11 settings x all 8 callers x all strains
python3 truvari_param_sweep.py --csv-out truvari_param_sweep_results.csv

# a single strain / subset of callers
python3 truvari_param_sweep.py --strains PfDd2 \
    --callers dysgu_loo dysgu_linear delly_loo pan_direct_sv \
    --csv-out truvari_sweep_PfDd2.csv

# named presets only
python3 truvari_param_sweep.py --settings current truvari_default pctovl_0.3

# optional full-factorial grid (interactions between parameters)
python3 truvari_param_sweep.py --grid-refdist 200 1000 2000 --grid-pctseq 0 0.3 0.7
```

- All path arguments from `full_pipeline.py` (`--outroot`, `--loo-base`, `--ref`, `--core-bed`, ...) are accepted with the same defaults.
- Giving any `--grid-*` argument turns off the named presets unless `--settings` is also given. Any grid dimension not given uses the `current` value.
- Run size: settings x strains x callers (default 11 x 12 x 8 = 1,056 Truvari runs).
- `--strains` defaults to the 12 held-out strains in `full_pipeline.py`, the same 12 plotted below.

### Output

- `--csv-out` (default `truvari_param_sweep_results.csv`), one row per strain, caller and setting:
  `strain, caller, setting, refdist, chunksize, pctseq, pctsize, pctovl, TP_base, TP_comp, FP, FN, precision, recall, f1, base_cnt, comp_cnt, error`
- Per-run Truvari output in `<outroot>/truvari_param_sweep/<strain>/<caller>/<setting>/` (override with `--sweep-out`).

---

## 3. `plot_truvari_sweep.py`

```bash
python3 plot_truvari_sweep.py truvari_param_sweep_results.csv [output.png]
```

- **Input:** the CSV from the sweep. Rows with an `error` value are dropped.
- **Strain filter:** only the 12 strains in `V2_STRAINS` (Pf7G8, PfCD01, PfDd2, PfGA01, PfGB4, PfGN01, PfHB3, PfIT, PfKE01, PfKH01, PfKH02, PfSN01) are plotted, matching the pipeline default. Edit the list if the panel changes.
- **Summary:** median across strains for each caller and setting (`AGG_FUNC`, can be set to `"mean"`).
- **Layout:** three stacked panels (precision, recall, F1), one line per caller, shared x-axis of parameter settings. The chosen setting is labelled `CHOSEN` and shaded. A single legend sits in the right margin. The header gives the number of strains plotted.
- **Output:** `truvari_sweep_plot.png` by default (150 dpi), or the second argument.
- **Configuration at the top of the script:** `SETTING_ORDER`, `CHOSEN_SETTING`, `PARAMS` (x-axis labels), `CALLERS`, `PALETTE` and font sizes. If a setting is added to the sweep, add it to `SETTING_ORDER` and `PARAMS` as well.

---
