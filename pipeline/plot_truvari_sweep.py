#!/usr/bin/env python3
"""
plot_truvari_sweep.py -- line-plot visualization of truvari_param_sweep.py
results (precision/recall/F1 per SV caller across parameter settings).

Usage:
    python3 plot_truvari_sweep.py truvari_param_sweep_results.csv

Only strains in V2_STRAINS are included (the actual pf3k_seq_v2.txt LOO
panel) -- edit that list if your panel changes. Uses median (not mean)
across strains per caller/setting, since the strain-level distributions
are right-skewed (a few high-F1 callers pull the mean up relative to the
median) -- median is the more robust summary here.
"""

import sys
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ------------------------------------------------------------------
# CONFIG -- edit these to match your data
# ------------------------------------------------------------------
V2_STRAINS = ["Pf7G8", "PfCD01", "PfDd2", "PfGA01", "PfGB4", "PfGN01",
              "PfHB3", "PfIT", "PfKE01", "PfKH01", "PfKH02", "PfSN01"]

SETTING_ORDER = ["current", "truvari_default", "refdist_200", "refdist_500", "refdist_2000",
                  "pctseq_0.3", "pctseq_0.7", "pctsize_0.5", "pctsize_0.7",
                  "pctovl_0.3", "moderate_all"]

# the setting actually used/reported in the paper -- gets the highlighted
# column and the CHOSEN_LABEL below instead of its raw setting name
CHOSEN_SETTING = "current"
CHOSEN_LABEL = "CHOSEN"

PARAMS = {
    "current":         "refdist=1000, pctseq=0.0, pctsize=0.3, pctovl=0.0",
    "truvari_default": "refdist=500, pctseq=0.7, pctsize=0.7, pctovl=0.0",
    "refdist_200":     "refdist=200, pctseq=0.0, pctsize=0.3, pctovl=0.0",
    "refdist_500":     "refdist=500, pctseq=0.0, pctsize=0.3, pctovl=0.0",
    "refdist_2000":    "refdist=2000, pctseq=0.0, pctsize=0.3, pctovl=0.0",
    "pctseq_0.3":      "refdist=1000, pctseq=0.3, pctsize=0.3, pctovl=0.0",
    "pctseq_0.7":      "refdist=1000, pctseq=0.7, pctsize=0.3, pctovl=0.0",
    "pctsize_0.5":     "refdist=1000, pctseq=0.0, pctsize=0.5, pctovl=0.0",
    "pctsize_0.7":     "refdist=1000, pctseq=0.0, pctsize=0.7, pctovl=0.0",
    "pctovl_0.3":      "refdist=1000, pctseq=0.0, pctsize=0.3, pctovl=0.3",
    "moderate_all":    "refdist=500, pctseq=0.3, pctsize=0.5, pctovl=0.0",
}

CALLERS = ["pan_direct_sv", "vg_hap_sv_loo", "manta_loo", "manta_linear",
           "dysgu_loo", "dysgu_linear", "delly_loo", "delly_linear"]

PALETTE = ["#2BAE84", "#3366CC", "#8153A6", "#E87DBF", "#FF7033",
           "#F4A736", "#E63946", "#3FB1C2", "#8ACB4A", "#A5426D", "#6EC4E8"]

# font sizes (bumped up throughout vs matplotlib defaults)
FS_TITLE = 18
FS_LABEL = 15
FS_TICK = 12
FS_XTICK = 11
FS_LEGEND = 12
FS_SUPTITLE = 12

AGG_FUNC = "median"  # or "mean"


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: python3 plot_truvari_sweep.py <results.csv> [output.png]")
    csv_path = Path(sys.argv[1])
    out_path = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("truvari_sweep_plot.png")

    df = pd.read_csv(csv_path)
    df_valid = df[df["error"].isna()].copy()
    df_valid = df_valid[df_valid["strain"].isin(V2_STRAINS)]
    n_strains = df_valid["strain"].nunique()

    xtick_labels = [f"{s}\n({PARAMS[s]})" for s in SETTING_ORDER]
    chosen_idx = SETTING_ORDER.index(CHOSEN_SETTING)
    xtick_labels[chosen_idx] = f"{CHOSEN_LABEL}\n{PARAMS[CHOSEN_SETTING]}"

    agg = (df_valid.groupby(["caller", "setting"])[["precision", "recall", "f1"]]
           .agg(AGG_FUNC).reset_index())
    agg["setting"] = pd.Categorical(agg["setting"], categories=SETTING_ORDER, ordered=True)
    agg = agg.sort_values(["caller", "setting"])

    fig, axes = plt.subplots(3, 1, figsize=(15, 18), sharex=True)
    metric_titles = {"precision": "Precision", "recall": "Recall", "f1": "F1"}

    for ax, metric in zip(axes, ["precision", "recall", "f1"]):
        for i, caller in enumerate(CALLERS):
            sub = agg[agg["caller"] == caller].sort_values("setting")
            ax.plot(sub["setting"].astype(str), sub[metric], marker="o",
                     label=caller, color=PALETTE[i % len(PALETTE)],
                     linewidth=2.5, markersize=8)
        ax.axvspan(chosen_idx - 0.5, chosen_idx + 0.5, color="#2BAE84", alpha=0.08)
        ax.set_ylabel(f"{AGG_FUNC.capitalize()} {metric_titles[metric]}\n(across {n_strains} strains)",
                       fontsize=FS_LABEL)
        ax.set_title(f"{metric_titles[metric]} by truvari parameter setting, per SV caller",
                      fontsize=FS_TITLE)
        ax.tick_params(axis="y", labelsize=FS_TICK)
        ax.grid(axis="y", alpha=0.3)

    axes[-1].set_xticks(range(len(SETTING_ORDER)))
    axes[-1].set_xticklabels(xtick_labels, rotation=40, ha="right", fontsize=FS_XTICK)

    # single shared legend in the right margin, not overlapping any panel
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="center left", bbox_to_anchor=(1.0, 0.5),
               fontsize=FS_LEGEND, frameon=False)

    fig.suptitle(
        f"n={n_strains} strains (pf3k_seq_v2.txt LOO panel only)",
        fontsize=FS_SUPTITLE, y=0.995,
    )
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"saved {out_path}")


if __name__ == "__main__":
    main()
