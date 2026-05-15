#!/usr/bin/env python3
"""
PfPan Pangenome Pipeline
========================
Automates the full pangenome construction and statistics workflow:
  - (Optional) Retrieve Pf3k assemblies
  - Step 1: Generate seq file
  - Step 2: Run Cactus pangenome
  - Step 3: Run pangenome statistics (panacus + vg stats)

Usage:
    python construct_PfPan.py [--download] [--skip-seq] [--skip-graph] [--skip-stats]

Requirements:
    - cactus virtualenv at CACTUS_VENV (see config below)
    - panacus, vg available on PATH (or loaded via module)
    - Existing pipeline scripts: wget_pf3k_fasta.py, make_seqtxt.py, replot_panacus.py
"""

import argparse
import logging
import os
import subprocess
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuration — edit these paths to match your environment
# ---------------------------------------------------------------------------

CACTUS_VENV        = "/mnt/storage13/nbillows/pangenome/cactus-bin-v2.9.3/venv-cactus-v2.9.3/bin/activate"
PIPELINE_DIR       = str(Path(__file__).resolve().parent)
WORK_DIR           = str(Path(__file__).resolve().parent.parent)  # cactus --workDir
PANGENOME_BASE_DIR = str(Path(__file__).resolve().parent.parent)  # root directory for outputs

# Genome / graph settings
GENOME_DIR         = str(Path(__file__).resolve().parent.parent / "genomes")
SEQ_FILE           = "pf3k_seq_v2.txt"
REFERENCE          = "Pf3D7"
GRAPH_NAME         = "PfPan_test"
JOB_STORE          = str(Path(__file__).resolve().parent.parent / f"PfPan_{REFERENCE}_js")
OUT_DIR            = str(Path(__file__).resolve().parent.parent / f"PfPan_{REFERENCE}_pan")
LOG_FILE           = str(Path(__file__).resolve().parent.parent / f"PfPan_{REFERENCE}_log.log")

# Cactus resource settings
CONS_CORES         = 8
MG_MEMORY          = "128Gi"

# Pf3k FTP path for assembly download
PF3K_FTP_URL       = (
    "https://ftp.sanger.ac.uk/pub/project/pathogens/Plasmodium/falciparum/"
    "PF3K/ReferenceGenomes_Version1/GENOMES/"
)

# Panacus settings
PANACUS_STATS_DIR  = str(Path(__file__).resolve().parent.parent / "stats_panacus")
HAPLOTYPES_FILE    = str(Path(__file__).resolve().parent.parent / "stats_panacus" / "paths.haplotypes.txt")  # must exist before running stats
HIST_LEVELS        = "1,2,1,1,1"
HIST_QUANTILES     = "0,0,1,0.5,0.1"

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler("pipeline.log"),
    ],
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run(cmd: str, shell: bool = True, check: bool = True, env=None, cwd=None) -> subprocess.CompletedProcess:
    """Run a shell command, stream output, and raise on failure."""
    log.info(f"Running: {cmd}")
    result = subprocess.run(
        cmd,
        shell=shell,
        check=check,
        text=True,
        env=env,
        cwd=cwd,
    )
    return result


def run_in_cactus_venv(cmd: str) -> subprocess.CompletedProcess:
    """
    Source the Cactus virtualenv and run a command in a single shell session.
    This is the equivalent of:
        source <venv>/bin/activate && <cmd> && deactivate
    """
    wrapped = f'bash -c "source {CACTUS_VENV} && {cmd} && deactivate"'
    log.info(f"[cactus-venv] {cmd}")
    return run(wrapped)


def ensure_dir(path: str) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


# ---------------------------------------------------------------------------
# Optional: Download assemblies
# ---------------------------------------------------------------------------

def download_assemblies():
    """Retrieve Pf3k assemblies from the Sanger FTP via wget_pf3k_fasta.py."""
    log.info("=== [Optional] Downloading Pf3k assemblies ===")
    ensure_dir(GENOME_DIR)
    run(
        f'python "{PIPELINE_DIR}/wget_pf3k_fasta.py" '
        f"--file_path {PF3K_FTP_URL}",
        cwd=GENOME_DIR,
    )
    log.info("Assemblies downloaded to: %s", GENOME_DIR)


# ---------------------------------------------------------------------------
# Step 1: Make seq file
# ---------------------------------------------------------------------------

def make_seq_file():
    """Generate the two-column seq file required by cactus-pangenome."""
    log.info("=== Step 1: Generating seq file (%s) ===", SEQ_FILE)
    run(
        f'python "{PIPELINE_DIR}/make_seqtxt.py" '
        f"--file_name {SEQ_FILE}",
        cwd=GENOME_DIR,
    )
    seq_path = Path(GENOME_DIR) / SEQ_FILE
    if not seq_path.exists():
        raise FileNotFoundError(f"Seq file was not created: {seq_path}")
    log.info("Seq file created: %s", seq_path)


# ---------------------------------------------------------------------------
# Step 2: Build pangenome graph
# ---------------------------------------------------------------------------

def build_pangenome():
    """
    Run cactus-pangenome inside the Cactus virtualenv.

    Key flags:
      --filter 2            : filter low-quality alignments
      --haplo               : generate haplotype index
      --giraffe clip filter : build Giraffe mapping indexes
      --viz --odgi          : visualisation and ODGI outputs
      --gbz clip filter full: GBZ outputs at multiple levels
      --gfa clip full       : GFA outputs
      --vcf --vcfReference  : call variants against the reference
    """
    log.info("=== Step 2: Building pangenome graph ===")

    seq_path = Path(GENOME_DIR) / SEQ_FILE

    cactus_cmd = (
        f"cactus-pangenome {JOB_STORE} {seq_path} "
        f"--outDir {OUT_DIR} "
        f"--outName {GRAPH_NAME} "
        f"--reference {REFERENCE} "
        f"--filter 2 "
        f"--haplo "
        f"--giraffe clip filter "
        f"--viz --odgi "
        f"--chrom-vg clip filter "
        f"--chrom-og "
        f"--gbz clip filter full "
        f"--gfa clip full "
        f"--vcf --vcfReference {REFERENCE} "
        f"--logFile {LOG_FILE} "
        f"--workDir {WORK_DIR} "
        f"--consCores {CONS_CORES} "
        f"--mgMemory {MG_MEMORY} "
        f"chrom-alignments"
    )

    run_in_cactus_venv(cactus_cmd)
    log.info("Pangenome graph written to: %s", OUT_DIR)


# ---------------------------------------------------------------------------
# Step 3a: Panacus — Base Pair coverage
# ---------------------------------------------------------------------------

def panacus_bp(gfa_path: str):
    """
    Run panacus hist + histgrowth for base-pair coverage, then plot.

    Growth curve parameters:
      -l 1,2,1,1,1      : coverage level thresholds
      -q 0,0,1,0.5,0.1  : quantile thresholds (core / soft-core / shell / private)
      -S                 : stratified output
      -a                 : cumulative (all-haplotype) output
    """
    bp_dir = Path(PANACUS_STATS_DIR) / "bp"
    ensure_dir(bp_dir)
    log.info("--- Panacus: base-pair coverage (output: %s) ---", bp_dir)

    run(f"panacus hist --count bp {gfa_path} > {bp_dir}/bp.hist")

    run(
        f"panacus histgrowth --count bp {gfa_path} "
        f"-l {HIST_LEVELS} "
        f"-q {HIST_QUANTILES} "
        f"-S -a "
        f"-s {HAPLOTYPES_FILE} "
        f"> {bp_dir}/bp.growth"
    )

    run(
        f"python {PIPELINE_DIR}/replot_panacus.py {bp_dir}/bp.growth "
        f"> {bp_dir}/panacus_growth_bp.pdf"
    )
    log.info("BP growth curve: %s/panacus_growth_bp.pdf", bp_dir)


# ---------------------------------------------------------------------------
# Step 3b: Panacus — Node coverage
# ---------------------------------------------------------------------------

def panacus_node(gfa_path: str):
    """Run panacus histgrowth + similarity at the graph-node level."""
    node_dir = Path(PANACUS_STATS_DIR) / "node"
    ensure_dir(node_dir)
    log.info("--- Panacus: node coverage (output: %s) ---", node_dir)

    run(
        f"panacus histgrowth --count node {gfa_path} "
        f"-l {HIST_LEVELS} "
        f"-q {HIST_QUANTILES} "
        f"-S -a "
        f"-s {HAPLOTYPES_FILE} "
        f"> {node_dir}/node.growth"
    )

    run(
        f"python {PIPELINE_DIR}/replot_panacus.py {node_dir}/node.growth "
        f"> {node_dir}/panacus_growth_node.pdf"
    )


# ---------------------------------------------------------------------------
# Step 3c: vg stats — Basic graph metrics
# ---------------------------------------------------------------------------

def vg_stats():
    """Run vg stats on the clipped GBZ graph for basic graph metrics."""
    gbz = Path(OUT_DIR) / f"{GRAPH_NAME}.gbz"
    out = Path(PANACUS_STATS_DIR) / "vg_stats.txt"

    if not gbz.exists():
        log.warning("GBZ not found, skipping vg stats: %s", gbz)
        return

    log.info("--- vg stats (graph: %s) ---", gbz)
    run(f"vg stats -z -l -N -E -s {gbz} > {out}")
    log.info("vg stats written to: %s", out)


# ---------------------------------------------------------------------------
# Step 3 (orchestrator): Pangenome statistics
# ---------------------------------------------------------------------------

def run_statistics():
    """
    Orchestrate all pangenome statistics steps:
      3a. Panacus base-pair coverage & growth curve
      3b. Panacus node coverage & growth curve + similarity
      3c. vg stats
    """
    log.info("=== Step 3: Pangenome statistics ===")
    ensure_dir(PANACUS_STATS_DIR)

    # Copy and decompress GFA into stats directory
    gfa_gz  = Path(OUT_DIR) / f"{GRAPH_NAME}.gfa.gz"
    gfa_out = Path(PANACUS_STATS_DIR) / f"{GRAPH_NAME}.gfa"

    if not gfa_out.exists():
        if not gfa_gz.exists():
            raise FileNotFoundError(f"GFA not found: {gfa_gz}")
        log.info("Decompressing GFA to %s ...", gfa_out)
        run(f"cp {gfa_gz} {PANACUS_STATS_DIR}/ && gunzip {PANACUS_STATS_DIR}/{GRAPH_NAME}.gfa.gz")
    else:
        log.info("GFA already decompressed: %s", gfa_out)

    # Generate haplotypes file from seq file if not present
    if not Path(HAPLOTYPES_FILE).exists():
        seq_path = Path(GENOME_DIR) / SEQ_FILE
        if not seq_path.exists():
            raise FileNotFoundError(
                f"Cannot generate haplotypes file: seq file not found at {seq_path}"
            )
        log.info("Generating haplotypes file from seq file: %s", seq_path)
        haplotypes = [line.split()[0] for line in seq_path.read_text().splitlines() if line.strip()]
        Path(HAPLOTYPES_FILE).write_text("\n".join(haplotypes) + "\n")
        log.info("Haplotypes file written to: %s", HAPLOTYPES_FILE)
    else:
        log.info("Haplotypes file already exists: %s", HAPLOTYPES_FILE)

    panacus_bp(str(gfa_out))
    panacus_node(str(gfa_out))
    vg_stats()

    log.info("Statistics complete. Outputs in: %s", PANACUS_STATS_DIR)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Pf3k Pangenome Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--download", action="store_true",
        help="(Optional) Download Pf3k assemblies from Sanger FTP first",
    )
    parser.add_argument(
        "--skip-seq",   action="store_true", help="Skip Step 1 (seq file creation)",
    )
    parser.add_argument(
        "--skip-graph", action="store_true", help="Skip Step 2 (cactus-pangenome)",
    )
    parser.add_argument(
        "--skip-stats", action="store_true", help="Skip Step 3 (panacus + vg stats)",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    if args.download:
        download_assemblies()

    if not args.skip_seq:
        make_seq_file()
    else:
        log.info("Skipping Step 1 (seq file creation)")

    if not args.skip_graph:
        build_pangenome()
    else:
        log.info("Skipping Step 2 (cactus-pangenome)")

    if not args.skip_stats:
        run_statistics()
    else:
        log.info("Skipping Step 3 (panacus + vg stats)")

    log.info("Pipeline complete.")


if __name__ == "__main__":
    main()