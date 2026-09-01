# Pf Pangenome Leave-One-Out (LOO) Validation Pipeline

Bash pipeline that, for each of 14 *P. falciparum* strains, builds a
strain-excluded ("leave-one-out") pangenome graph, maps that strain's reads
back to it, calls short variants and structural variants (SVs) via several
methods, and benchmarks each caller against assembly-derived truth sets,
stratified by genome region.

Ran using:

```
bash full_pipeline.sh 
```
This is currently hardcoded. Will make more flexible in the future.

Four calling strategies are compared per strain:

| Short variants        | Structural variants |
|------------------------|----------------------|
| `pan_direct`  (external, pre-computed) | `pan_direct_sv` (external, pre-computed) |
| `gatk_linear` (GATK on linear BWA mapping) | `delly_linear` (Delly on linear BWA mapping) |
| `gatk_surject_loo` (GATK on LOO-graph-surjected mapping) | `delly_loo` (Delly on LOO-graph-surjected mapping) |
| `vg_hap_loo` (`vg call` on the LOO graph directly) | `vg_hap_sv_loo` (`vg call` on the LOO graph directly) |

---

## Requirements

- conda/mamba, with environments named: `cactus`, `fastq2matrix`, `truvari`,
  and a Delly environment (name currently unconfirmed — see **Open items**).
- Tools invoked: `cactus-pangenome`, `vg` (giraffe/surject/pack/call/paths),
  `samtools`, `bcftools`, `tabix`, `gatk` (HaplotypeCaller,
  GenotypeGVCFs, SelectVariants, VariantFiltration), `rtg` (vcfeval,
  vcfdecompose), `truvari`, `delly`, `bedtools`, `vt` (`annotate_indels`),
  `setGT.py` (custom script, must be on `PATH`).
- `set -euo pipefail` is active — the script aborts on the first error or
  unset variable.

---

## Inputs

### Configured paths (edit at top of script)

| Variable | What it is | Required? |
|---|---|---|
| `STRAINS` | Array of 14 strain names to process | fixed list |
| `SEQ_FILE` | Pf3k sequence manifest used by `cactus-pangenome` (list of genome FASTA paths, one per strain, used to build each LOO graph) | yes |
| `READS_DIR` | Directory of per-strain trimmed Illumina reads (`{strain}_1/2.trimmed.fastq.gz`), linear BQSR CRAMs (`{strain}.bqsr.cram`), and raw linear GATK GVCFs (`{strain}.g.vcf.gz`) | yes |
| `LINEAR_GVCF_DIR` | = `READS_DIR`; where raw linear GVCFs live | yes |
| `LOO_BASE` | Base directory for per-strain LOO graph-build/mapping/calling working directories (`loo_{STRAIN}/`) | yes |
| `COMPARISON_DIR` | Pre-existing directory of external "pan_direct" calls and truth VCFs (`{strain}.paftools.snps_indels.vcf.gz`, `{strain}.pan.snps_indels.vcf.gz`, `{strain}.pan.svs.vcf.gz`, `{strain}.svim.svs.vcf.gz`) | yes |
| `CACTUS_VENV` | Path to the Cactus Python venv `activate` script | yes |
| `REF` | Pf3D7 v3 reference FASTA (used for graph building, `bcftools norm`, Delly, `samtools faidx`) | yes |
| `REF_GATK` | Reference FASTA used in GATK steps — **currently believed identical to `REF`, needs confirming** | yes |
| `REF_SDF` | RTG `.sdf` index of the reference (for `rtg vcfeval`) | yes |
| `VT_BIN` | Path to the `vt` binary (for `annotate_indels`) | yes |
| `GFF` | Genome annotation GFF3 — **placeholder path, must be set** | yes |
| `BEDDIR` | Directory of genome-partition/category BED files used to stratify all benchmarking (e.g. genic, subtelomeric, etc.) | yes |
| `VSA_BED` | Output path for the variable-surface-antigen (VAR/RIF/STEVOR/etc.) gene mask; **built automatically from `GFF` if it doesn't already exist** | auto-built |
| `TRF_BED` | Pre-computed sensitive tandem-repeat (TRF) BED — **must already exist**, script exits if missing | yes |
| `CONF_BED_DIR` | Directory of per-strain SVIM-asm "confident regions" BEDs (`{strain}.April2018/{strain}.April2018.confident.bed`) | yes |
| `OUTROOT` | Root output directory for all pipeline results (created if missing) | yes |

### Tunable parameters

| Variable | Meaning | Default |
|---|---|---|
| `VSA_PAD` | bp padding added around VSA gene annotations before merging into `VSA_BED` | 500 |
| `VT_MIN_FZ_RL` | Minimum repeat-tract length (`vt`'s `FZ_RL`) for an indel to be flagged as sitting in an ambiguous repeat tract | 4 |
| `BOUNDARY_BUFFER` | bp trimmed off both ends of each confident region, to avoid edge artifacts | 1000 |
| `DENSITY_WINDOW` | Window size (bp) used to compute local truth-variant density | 200 |
| `DENSITY_PERCENTILE` | Density percentile above which a window is excluded as "high density" | 95 |
| `MIN_BP` | Minimum total bp a BED category must cover to be benchmarked; smaller categories are skipped | 1000 |

### Stage toggles

| Variable | Controls | Default |
|---|---|---|
| `RUN_MAPPING` | LOO graph build (`cactus-pangenome`), read mapping (`giraffe`), surject, `vg pack`/`vg call`, GATK-on-surjected-BAM | `false` (assumed already done) |
| `RUN_SMALL_VAR` | GATK filtering (linear + LOO-surject) + `vg_hap_loo` filtering + confidence-tier building + `vcfeval` benchmarking | `false` (assumed already done) |
| `RUN_SV` | Delly (LOO + linear) + `truvari` SV benchmarking | `true` |

Set the relevant toggle(s) to `true` to (re)run that stage for the strains
in `STRAINS`; each stage's steps check for existing output files and skip
recomputation where they can (e.g. `[[ ! -f "${GBZ}" ]]`).

---

## Pipeline stages

### 0. `build_shared_beds` (once, before the strain loop)
- **Builds `VSA_BED`** (if absent) from `GFF`: extracts `gene` features named
  `VAR`/`VAR-like`/`VAR2CSA`/`RIF`/`RIFA`/`MC-2TM`/`SURF`, plus
  `repeat_region` and `centromere` features, pads by `VSA_PAD` bp, merges.
- **Checks `TRF_BED` exists** (exits with error if not).
- **Builds `${OUTROOT}/genome.txt`**, a `samtools faidx`-derived chrom-size
  file used by `bedtools makewindows` later.

### 1. `run_loo` (per strain, if `RUN_MAPPING=true`)
Builds the strain-excluded pangenome graph and calls variants from that
strain's own reads mapped back to it.
- **In:** `SEQ_FILE` (strain line removed), `CACTUS_VENV`, `READS_DIR/{strain}_{1,2}.trimmed.fastq.gz`, `REF_GATK`.
- **Out** (all under `${LOO_BASE}/loo_${STRAIN}/`):
  - `pf3k_seq_no_${STRAIN}.txt` — filtered manifest
  - `PfPan_no_${STRAIN}_pan/` — Cactus graph outputs (`.gbz`, `.snarls`, etc.)
  - `${STRAIN}.loo.gaf.gz` — Giraffe alignment
  - `Pf3D7.paths.txt` — Pf3D7 path names extracted from the graph
  - `${STRAIN}.loo.sort.bam` (+`.bai`) — surjected, sorted, indexed BAM in Pf3D7 linear coordinates
  - `${STRAIN}.loo.pack` — `vg pack` coverage index
  - `${STRAIN}.loo.SV.call_haploid_.vcf.gz` — `vg call` haploid VCF, contig-renamed via `fix_contigs`
  - `${STRAIN}.loo.gatk_surject.g.vcf.gz` — GATK HaplotypeCaller GVCF on the surjected BAM

### `fix_contigs` (helper, called from `run_loo`)
Strips the `Pf3D7#0#` path-name prefix from a VCF's header and body.
- **In:** VCF path, a tag used for temp filenames.
- **Out:** rewrites the VCF in place (+ `.tbi`); temp files removed.

### `filter_gatk` (shared helper, called from steps 2 and 3)
Standard GATK short-variant filter chain: genotype GVCF → split
SNPs/indels → hard-filter → concat → PASS-only → normalize → restrict to
variants `<50 bp` → genotype-refine (`setGT.py --fraction 0.7`) → decompose
MNPs/indels (`rtg vcfdecompose`).
- **In:** GVCF, sample name, working directory, final output path.
- **Out:** many intermediates in the working directory, ending in the
  caller's `*.decomposed.vcf.gz` (+ `.tbi`).

### 2. `run_gatk_linear` (per strain, if `RUN_SMALL_VAR=true`)
- **In:** `${LINEAR_GVCF_DIR}/${STRAIN}.g.vcf.gz`
- **Out:** `${OUTROOT}/${STRAIN}/gatk_linear/${STRAIN}.GATK_linear.shortvars.lt50bp.PASS.GT.decomposed.vcf.gz`

### 3. `run_gatk_surject_loo` (per strain, if `RUN_SMALL_VAR=true`)
- **In:** `${LOO_BASE}/loo_${STRAIN}/${STRAIN}.loo.gatk_surject.g.vcf.gz` (from step 1)
- **Out:** `${OUTROOT}/${STRAIN}/gatk_surject_loo/${STRAIN}.GATK_surject_loo.shortvars.lt50bp.PASS.GT.decomposed.vcf.gz`

### 4. `run_vg_hap_loo_filter` (per strain, if `RUN_SMALL_VAR=true`)
Filters the `vg call` pangenome VCF from step 1 into separate short-variant
and SV callsets.
- **In:** `${LOO_BASE}/loo_${STRAIN}/${STRAIN}.loo.SV.call_haploid_.vcf.gz`, `REF`
- **Out** (in `loo_${STRAIN}/`):
  - `${STRAIN}.loo.pan.shortvars.lt50bp.PASS.hap.dp5.1_1_only.decomposed.vcf.gz` — **final short-variant** pangenome calls (normalized, PASS, depth >5, homozygous-alt only, decomposed)
  - `${STRAIN}.loo.pan.SV.ge50bp.PASS.hap.vcf.gz` — **final SV** pangenome calls (≥50 bp, PASS)

### 5 / 5b. `run_delly_loo` / `run_delly_linear` (per strain, if `RUN_SV=true`)
SV calling with Delly on the LOO-surjected BAM and on the true linear
(BWA/BQSR) baseline, respectively.
- **In:** `${STRAIN}.loo.sort.bam` (from step 1) *or* `${READS_DIR}/${STRAIN}.bqsr.cram`, `REF`
- **Out:** `${STRAIN}.loo.delly.SV.ge50bp.PASS.vcf.gz` (in `loo_${STRAIN}/`) *or* `${OUTROOT}/${STRAIN}/delly_linear/${STRAIN}.delly_linear.SV.ge50bp.PASS.vcf.gz`

### `reheader_truth` (helper, called from steps 6 and 7)
Renames the sample column of the paftools truth VCF to the strain name.
- **In:** `${COMPARISON_DIR}/${STRAIN}.paftools.snps_indels.vcf.gz`
- **Out:** `${OUTROOT}/${STRAIN}/${STRAIN}.paftools.snps_indels.reheader.vcf.gz` (+`.tbi`); sets `TRUTH_VCF_RENAMED`

### 6. `build_confidence` (per strain, if `RUN_SMALL_VAR=true`)
Builds the strict, multi-layer confident-regions BED used for short-variant
benchmarking, and intersects every category BED in `BEDDIR` with it.
- **In:** reheadered truth VCF, `${CONF_BED_DIR}/${STRAIN}.April2018/${STRAIN}.April2018.confident.bed` (skips strain if missing), `genome.txt`, `VSA_BED`, `TRF_BED`, the four decomposed VCFs from steps 2–4 (+truth), `BEDDIR/*.bed`
- **Out** (in `${OUTROOT}/${STRAIN}/beds/`):
  - progressively refined confident BEDs (`*.confident.merged` → `*.trimmed` → `*.strict` → `*.strict.vsa_masked` → **`*.confident.FINAL.bed`**), removing SVIM-asm-flagged regions, boundary edges, high truth-variant-density windows, VSA genes, and a repeat-tract mask (TRF ∪ ambiguous-indel regions from all 4 callers)
  - `${BEDNAME}.FINAL.bed` for every category BED — used by `run_vcfeval_all`
  - *(SVs deliberately skip this confidence tier — `run_truvari_all` uses the plain `BEDDIR` BEDs directly, matching the original non-LOO pipeline's stringency)*

### 7. `run_vcfeval_all` (per strain, if `RUN_SMALL_VAR=true`)
Benchmarks all 4 short-variant callers against the paftools truth with
`rtg vcfeval`, per genome category.
- **In:** truth VCF, the 4 caller VCFs (`pan_direct` from `COMPARISON_DIR`; `gatk_linear`, `gatk_surject_loo`, `vg_hap_loo` from steps 2–4), `${BEDNAME}.FINAL.bed` per category (categories below `MIN_BP` skipped), `REF_SDF`
- **Out:** `${OUTROOT}/${STRAIN}/vcfeval/${CALLER}/${BEDNAME}/` — standard `rtg vcfeval` output directory (tp/fp/fn VCFs, `summary.txt`, ROC data) for every caller × category combination

### 8. `run_truvari_all` (per strain, if `RUN_SV=true`)
Benchmarks all 4 SV callers against the SVIM-asm truth with `truvari
bench`, per genome category.
- **In:** SV truth (`${COMPARISON_DIR}/${STRAIN}.svim.svs.vcf.gz`, filtered to non-ref genotypes once), the 4 SV caller VCFs (`pan_direct_sv` from `COMPARISON_DIR`; `delly_loo`, `vg_hap_sv_loo` from steps 4–5; `delly_linear` from step 5b), each filtered to non-ref genotypes, plain `${BEDDIR}/*.bed` categories (below `MIN_BP` skipped), `REF`
- **Out:** `${OUTROOT}/${STRAIN}/truvari/${CALLER}/${BEDNAME}/` — standard `truvari bench` output directory (`tp-base.vcf.gz`, `tp-comp.vcf.gz`, `fp.vcf.gz`, `fn.vcf.gz`, `summary.json`) for every caller × category combination

---

## Output directory structure

```
${OUTROOT}/
├── genome.txt
├── <STRAIN>/
│   ├── <STRAIN>.paftools.snps_indels.reheader.vcf.gz(.tbi)
│   ├── gatk_linear/         # step 2 intermediates + final decomposed VCF
│   ├── gatk_surject_loo/    # step 3 intermediates + final decomposed VCF
│   ├── delly_linear/        # step 5b BCF + filtered SV VCF
│   ├── beds/                # step 6: confidence tiers + *.FINAL.bed per category
│   ├── vcfeval/<CALLER>/<CATEGORY>/   # step 7 benchmarking output
│   └── truvari/<CALLER>/<CATEGORY>/   # step 8 benchmarking output
└── ...

${LOO_BASE}/loo_<STRAIN>/     # step 1 graph build, mapping, vg call, delly_loo,
                              # and step 4 vg_hap_loo filtered VCFs
```

---

## Open items (flagged `CONFIRM` in the script)

- `GFF` is a placeholder path (`/path/to/...`) — must be set before running `build_shared_beds`.
- `REF_GATK` is assumed identical to `REF`, just historically referenced via a different path — worth verifying they're the same file.
- The Delly conda environment name is unconfirmed in both `run_delly_loo` and `run_delly_linear`.
- `pan_direct` (external short-variant calls used in `run_vcfeval_all`) is used as-is — unclear whether it needs additional filtering before comparison.
