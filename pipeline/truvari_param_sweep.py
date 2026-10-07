#!/usr/bin/env python3
"""
truvari_param_sweep.py -- SV benchmarking parameter sensitivity analysis

Addresses reviewer comment: "The SV benchmarking parameters require
justification. The Truvari settings appear permissive, including 1-kb
positional tolerance and no required reciprocal overlap or sequence
similarity. The authors should explain how these thresholds were selected
and ideally demonstrate robustness to alternative settings."

What this does: reruns `truvari bench` for a chosen set of SV callers
against the same truth set and core-genome bed already used by
full_pipeline.py's run_truvari_all(), at several parameter settings --
the paper's chosen settings, truvari's own out-of-the-box defaults, and
settings that isolate the effect of relaxing each individual parameter
one at a time (holding the others at the paper's chosen values). This
lets you report, per parameter, how sensitive precision/recall/F1 are to
that specific choice, rather than just asserting the settings were
reasonable.

INPUTS: this script does NOT rebuild or re-filter any caller VCFs -- it
reuses the exact same already-built, already-haploid-filtered
*.alt_only.mixed.vcf.gz files that run_truvari_all() produces, so the
sweep is a clean comparison of *truvari settings only*, not a mix of
different filtering and different settings. Run full_pipeline.py's SV
stages first for the strains/callers you want to sweep.

USAGE (run from the same directory as full_pipeline.py):
    python3 truvari_param_sweep.py --strains PfDd2 \
        --callers dysgu_loo dysgu_linear delly_loo pan_direct_sv \
        --csv-out truvari_sweep_PfDd2.csv

All of full_pipeline.py's path arguments (--outroot, --loo-base, --ref,
--core-bed, etc.) are available here too, with the same defaults, so this
only needs overriding if you're pointing at a non-default location.
"""

import argparse
import csv
import json
from pathlib import Path

import full_pipeline as fp  # reuse path config, conda_sh, is_file, etc.

############################################################
# PARAMETER SETTINGS TO SWEEP
# "current" matches exactly what full_pipeline.py's run_truvari_all uses.
# "truvari_default" is truvari's own out-of-the-box defaults (refdist=500,
# pctseq=0.7, pctsize=0.7; pctovl actually defaults to 0.0 in truvari
# itself too -- worth noting in the rebuttal that "no required reciprocal
# overlap" is truvari's own default behaviour, not something loosened
# specifically for this paper). The rest each change ONE parameter at a
# time away from "current", to isolate that parameter's individual effect.
############################################################
SETTINGS = {
    "current":         dict(refdist=1000, pctseq=0.0, pctsize=0.3, pctovl=0.0),
    "truvari_default": dict(refdist=500,  pctseq=0.7, pctsize=0.7, pctovl=0.0),
    "refdist_200":      dict(refdist=200,  pctseq=0.0, pctsize=0.3, pctovl=0.0),
    "refdist_500":     dict(refdist=500,  pctseq=0.0, pctsize=0.3, pctovl=0.0),
    "refdist_2000":    dict(refdist=2000, pctseq=0.0, pctsize=0.3, pctovl=0.0),
    "pctseq_0.3":      dict(refdist=1000, pctseq=0.3, pctsize=0.3, pctovl=0.0),
    "pctseq_0.7":      dict(refdist=1000, pctseq=0.7, pctsize=0.3, pctovl=0.0),
    "pctsize_0.5":     dict(refdist=1000, pctseq=0.0, pctsize=0.5, pctovl=0.0),
    "pctsize_0.7":     dict(refdist=1000, pctseq=0.0, pctsize=0.7, pctovl=0.0),
    "pctovl_0.3":      dict(refdist=1000, pctseq=0.0, pctsize=0.3, pctovl=0.3),
    "moderate_all":    dict(refdist=500,  pctseq=0.3, pctsize=0.5, pctovl=0.0),
}

DEFAULT_CALLERS = [
    "pan_direct_sv", "delly_loo", "delly_linear",
    "manta_loo", "manta_linear", "dysgu_loo", "dysgu_linear",
    "vg_hap_sv_loo",
]


def caller_alt_vcf(strain: str, caller: str) -> Path:
    """Path to the already-built, already-haploid-filtered alt_only.mixed
    VCF for this caller -- the exact same input run_truvari_all() uses in
    full_pipeline.py, so this sweep compares truvari settings on identical
    inputs rather than re-deriving its own filtering."""
    base = {
        "pan_direct_sv": fp.COMPARISON_DIR / f"{strain}.pan.svs.vcf.gz",
        "delly_loo": fp.LOO_BASE / f"loo_{strain}" / f"{strain}.loo.delly.SV.ge50bp.PASS.vcf.gz",
        "delly_linear": fp.OUTROOT / strain / "delly_linear" / f"{strain}.delly_linear.SV.ge50bp.PASS.vcf.gz",
        "manta_loo": fp.LOO_BASE / f"loo_{strain}" / f"{strain}.loo.manta.SV.ge50bp.PASS.vcf.gz",
        "manta_linear": fp.OUTROOT / strain / "manta_linear" / f"{strain}.manta_linear.SV.ge50bp.PASS.vcf.gz",
        "dysgu_loo": fp.LOO_BASE / f"loo_{strain}" / f"{strain}.loo.dysgu.SV.ge50bp.PASS.vcf.gz",
        "dysgu_linear": fp.OUTROOT / strain / "dysgu_linear" / f"{strain}.dysgu_linear.SV.ge50bp.PASS.vcf.gz",
        "vg_hap_sv_loo": fp.LOO_BASE / f"loo_{strain}" / f"{strain}.loo.pan.SV.ge50bp.PASS.hap.vcf.gz",
    }.get(caller)
    if base is None:
        raise ValueError(f"unknown caller '{caller}'")
    base_str = str(base)
    assert base_str.endswith(".vcf.gz")
    return Path(base_str[: -len(".vcf.gz")] + ".alt_only.mixed.vcf.gz")


def run_one(strain: str, caller: str, setting_name: str, params: dict, sweep_out: Path) -> dict:
    truth = fp.OUTROOT / strain / f"{strain}.svim.svs.alt_only.mixed.vcf.gz"
    calls = caller_alt_vcf(strain, caller)

    if not fp.is_file(truth):
        return {"error": f"missing truth vcf: {truth}"}
    if not fp.is_file(calls):
        return {"error": f"missing calls vcf: {calls}"}

    out_dir = sweep_out / strain / caller / setting_name
    out_dir.parent.mkdir(parents=True, exist_ok=True)
    fp.rmtree_if_exists(out_dir)

    # truvari requires --chunksize >= --refdist; full_pipeline.py's
    # run_truvari_all() hardcodes both at 1000 since refdist is fixed
    # there, but this sweep varies refdist, so chunksize has to track it
    chunksize = max(1000, params["refdist"])

    fp.conda_sh("truvari", """truvari bench -b %(truth)s -c %(calls)s -o %(out)s -f %(ref)s \
        -r %(refdist)s -C %(chunksize)s -O %(pctovl)s -p %(pctseq)s -P %(pctsize)s -s 50 -S 15 \
        --sizemax 10000 --includebed %(bed)s""" % {
        "truth": truth, "calls": calls, "out": out_dir, "ref": fp.REF,
        "refdist": params["refdist"], "chunksize": chunksize, "pctovl": params["pctovl"],
        "pctseq": params["pctseq"], "pctsize": params["pctsize"],
        "bed": fp.CORE_BED,
    })

    summary_path = out_dir / "summary.json"
    if not fp.is_file(summary_path):
        return {"error": "truvari did not produce summary.json"}
    return json.loads(summary_path.read_text())


def build_grid_settings(refdists, pctseqs, pctsizes, pctovls) -> dict:
    """Full factorial combination of the given parameter values, so
    interaction effects (e.g. does pctseq=0.3 behave differently at
    refdist=200 vs refdist=2000?) can be seen, not just one-parameter-at-
    a-time isolation like the named SETTINGS above."""
    import itertools
    grid = {}
    for r, seq, size, ovl in itertools.product(refdists, pctseqs, pctsizes, pctovls):
        name = f"grid_r{r}_seq{seq}_size{size}_ovl{ovl}"
        grid[name] = dict(refdist=r, pctseq=seq, pctsize=size, pctovl=ovl)
    return grid


def main():
    ap = fp.build_arg_parser()  # inherits all of full_pipeline's path args/defaults
    ap.add_argument("--callers", nargs="+", default=DEFAULT_CALLERS,
                     help="Which SV callers to sweep (must already have alt_only VCFs built)")
    ap.add_argument("--settings", nargs="+", default=None,
                     help="Named presets to run (see SETTINGS dict for the list). "
                          "Default: all named presets, UNLESS any --grid-* argument "
                          "is given, in which case the default becomes none (grid "
                          "combinations only) -- pass --settings explicitly to run "
                          "both named presets and a grid together.")
    ap.add_argument("--grid-refdist", nargs="+", type=int, default=None,
                     help="Full factorial sweep: refdist values to combine with "
                          "--grid-pctseq/--grid-pctsize/--grid-pctovl. Any dimension "
                          "not given defaults to current's single value for that "
                          "parameter.")
    ap.add_argument("--grid-pctseq", nargs="+", type=float, default=None)
    ap.add_argument("--grid-pctsize", nargs="+", type=float, default=None)
    ap.add_argument("--grid-pctovl", nargs="+", type=float, default=None)
    ap.add_argument("--sweep-out", type=Path, default=None,
                     help="Where to write per-run truvari outputs (default: <outroot>/truvari_param_sweep)")
    ap.add_argument("--csv-out", type=Path, default=Path("truvari_param_sweep_results.csv"))
    args = ap.parse_args()
    fp.apply_args(args)

    grid_dims_given = any([args.grid_refdist, args.grid_pctseq, args.grid_pctsize, args.grid_pctovl])

    run_settings = {}
    if args.settings is not None:
        for name in args.settings:
            if name not in SETTINGS:
                print(f"WARNING: unknown preset '{name}', skipping (known: {list(SETTINGS.keys())})")
                continue
            run_settings[name] = SETTINGS[name]
    elif not grid_dims_given:
        run_settings.update(SETTINGS)  # default when nothing else specified: all named presets

    if grid_dims_given:
        current = SETTINGS["current"]
        refdists = args.grid_refdist or [current["refdist"]]
        pctseqs = args.grid_pctseq or [current["pctseq"]]
        pctsizes = args.grid_pctsize or [current["pctsize"]]
        pctovls = args.grid_pctovl or [current["pctovl"]]
        grid = build_grid_settings(refdists, pctseqs, pctsizes, pctovls)
        print(f"Grid: {len(refdists)} refdist x {len(pctseqs)} pctseq x "
              f"{len(pctsizes)} pctsize x {len(pctovls)} pctovl = {len(grid)} combinations")
        run_settings.update(grid)

    total_runs = len(run_settings) * len(fp.STRAINS) * len(args.callers)
    print(f"{len(run_settings)} settings x {len(fp.STRAINS)} strains x "
          f"{len(args.callers)} callers = {total_runs} truvari runs planned\n")

    sweep_out = args.sweep_out or (fp.OUTROOT / "truvari_param_sweep")
    sweep_out.mkdir(parents=True, exist_ok=True)

    rows = []
    for strain in fp.STRAINS:
        for caller in args.callers:
            for setting_name, params in run_settings.items():
                print(f"-> {strain} | {caller} | {setting_name} "
                      f"(r={params['refdist']} pctseq={params['pctseq']} "
                      f"pctsize={params['pctsize']} pctovl={params['pctovl']})")
                result = run_one(strain, caller, setting_name, params, sweep_out)

                row = {
                    "strain": strain, "caller": caller, "setting": setting_name,
                    "refdist": params["refdist"], "chunksize": max(1000, params["refdist"]),
                    "pctseq": params["pctseq"], "pctsize": params["pctsize"], "pctovl": params["pctovl"],
                }
                if "error" in result:
                    print(f"   SKIP: {result['error']}")
                    row["error"] = result["error"]
                else:
                    row.update({
                        "TP_base": result.get("TP-base"),
                        "TP_comp": result.get("TP-comp"),
                        "FP": result.get("FP"),
                        "FN": result.get("FN"),
                        "precision": result.get("precision"),
                        "recall": result.get("recall"),
                        "f1": result.get("f1"),
                        "base_cnt": result.get("base cnt"),
                        "comp_cnt": result.get("comp cnt"),
                    })
                rows.append(row)

    fieldnames = ["strain", "caller", "setting", "refdist", "chunksize", "pctseq", "pctsize", "pctovl",
                  "TP_base", "TP_comp", "FP", "FN", "precision", "recall", "f1",
                  "base_cnt", "comp_cnt", "error"]
    with open(args.csv_out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    print(f"\nDONE. {len(rows)} runs written to {args.csv_out}")


if __name__ == "__main__":
    main()