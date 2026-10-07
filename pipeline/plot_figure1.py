#!/usr/bin/env python3
"""
Figure 1 (bp, panels A-E) and the supplementary node figure (panels A-C only)

  A  coverage histogram                      (panacus hist)
  B  cumulative pangenome growth             (panacus histgrowth)
  C  novel sequence added per genome         (derived from B)
  D  small variants (<50 bp), core genome    (vcfeval F-measure per strain)
  E  structural variants (>=50 bp), core genome (truvari F-measure per strain)

Usage (run in the folder holding the input files):
    python plot_figure1.py bp      # MAIN Figure 1: A-C from bp.hist + bp.growth, plus D and E
    python plot_figure1.py node    # SUPPLEMENTARY figure: A-C only, from node.hist + node.growth

Optional (bp only; not needed for node):
    --eval  vcfeval_raw_strain_data_v2.csv     per-strain F-measures (D, E)
    --stats paired_comparison_stats_v3.csv     BH-adjusted paired Wilcoxon p-values (brackets in D, E)

Inputs for A-C:
    panacus hist       --count {bp,node} -s paths.haplotypes.txt -S PfPan.gfa > {count}.hist
    panacus histgrowth --count {bp,node} -l 1,2,1,1,1 -q 0,0,1,0.5,0.1 -S -a \
                       -s paths.haplotypes.txt PfPan.gfa > {count}.growth

Outputs: Figure1_bp.tiff / FigureS_node.tiff (300 dpi, LZW), plus .png and .pdf

Fonts: the figure is drawn 15 in wide. At a 7 in (full-width) print size the
smallest text (12 pt here) comes out at ~5.6 pt and the panel titles at ~7 pt.
"""
import argparse
from io import StringIO

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ap = argparse.ArgumentParser()
ap.add_argument('count', choices=['bp', 'node'])
ap.add_argument('--eval', default='vcfeval_raw_strain_data_v2.csv')
ap.add_argument('--stats', default='paired_comparison_stats_v3.csv')
args = ap.parse_args()
COUNT = args.count
UNIT = '#bp' if COUNT == 'bp' else '#nodes'
WITH_EVAL = COUNT == 'bp'      # D and E (variant-calling benchmarks) only go in the bp figure

# ------------------------------------------------------------------ fonts ----
FONT = dict(title=15.5, axis_title=14.5, tick=12.5, xtick_eval=12, legend=12,
            legend_small=11.5, stars=22, ns=15)

plt.rcParams.update({
    'font.size': FONT['tick'],
    'axes.titlesize': FONT['title'],
    'axes.labelsize': FONT['axis_title'],
    'xtick.labelsize': FONT['tick'],
    'ytick.labelsize': FONT['tick'],
    'legend.fontsize': FONT['legend'],
    'axes.unicode_minus': False,
})

# ---------------------------------------------------------------- palette ----
GREEN, BLUE, PURPLE, PINK, ORANGE = "#2BAE84", "#3366CC", "#8153A6", "#E87DBF", "#FF7033"
GROWTH_COLORS = [GREEN, BLUE, PURPLE, PINK, ORANGE]

SNV_LEVELS = [
    "Linear Mapping\n(GATK)",
    "Pangenome Graph\n(LOOV: Surject + GATK)",
    "Pangenome Graph\n(LOOV: VG Call)",
    "Pangenome Graph\n(Cactus)",
]
SV_LEVELS = [
    "Linear Mapping\n(Delly)", "Pangenome Graph\n(LOOV: Delly Surject)",
    "Linear Mapping\n(Manta)", "Pangenome Graph\n(LOOV: Manta Surject)",
    "Linear Mapping\n(Dysgu)", "Pangenome Graph\n(LOOV: Dysgu Surject)",
    "Pangenome Graph\n(LOOV: VG Call)", "Pangenome Graph\n(Cactus)",
]
APPROACH_COLOURS = {
    "Linear Mapping\n(GATK)": "#3366CC",
    "Pangenome Graph\n(LOOV: Surject + GATK)": "#2BAE84",
    "Pangenome Graph\n(LOOV: VG Call)": "#8153A6",
    "Pangenome Graph\n(Cactus)": "#FF7033",
    "Linear Mapping\n(Delly)": "#3366CC",
    "Pangenome Graph\n(LOOV: Delly Surject)": "#2BAE84",
    "Linear Mapping\n(Manta)": "#F4A736",
    "Linear Mapping\n(Dysgu)": "#E63946",
    "Pangenome Graph\n(LOOV: Manta Surject)": "#8ACB4A",
    "Pangenome Graph\n(LOOV: Dysgu Surject)": "#A5426D",
}
# legend: two rows, filled column-wise -> each column is a linear / graph pair
LEGEND_ORDER = [
    "Linear Mapping\n(GATK)", "Pangenome Graph\n(LOOV: Surject + GATK)",
    "Linear Mapping\n(Delly)", "Pangenome Graph\n(LOOV: Delly Surject)",
    "Linear Mapping\n(Manta)", "Pangenome Graph\n(LOOV: Manta Surject)",
    "Linear Mapping\n(Dysgu)", "Pangenome Graph\n(LOOV: Dysgu Surject)",
    "Pangenome Graph\n(LOOV: VG Call)", "Pangenome Graph\n(Cactus)",
]
CORE_ACCENT = GREEN


# ------------------------------------------------------------- helpers -------
def read_panacus(path, **kw):
    with open(path) as fh:
        lines = "".join(l for l in fh if not l.startswith('#'))
    return pd.read_csv(StringIO(lines), **kw)


def humanize_number(i, precision=0):
    order, x = 0, i
    if abs(i) > 0:
        order = int(np.log10(abs(i))) // 3
        x = i / 10**(order * 3)
    return '{:,.{prec}f}{:}'.format(x, ['', 'K', 'M', 'B', 'D'][order], prec=precision)


def calibrate_yticks_text(yticks):
    prec = 0
    txt = [humanize_number(y, prec) for y in yticks]
    while len(set(txt)) < len(txt):
        prec += 1
        txt = [humanize_number(y, prec) for y in yticks]
    return txt


def style_y_human(ax):
    yt = ax.get_yticks()
    ax.set_yticks(yt)
    ax.set_yticklabels(calibrate_yticks_text(yt))


def int_xticks(ax, labels):
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels)


# -------------------------------------------------- data for panels A-C ------
hist = read_panacus(f'{COUNT}.hist', sep='\t', header=[0, 1], index_col=0)
hist.columns = [COUNT]
hist = hist[hist.index.notna()]
hist.index = hist.index.astype(int)
hist = hist.sort_index()
hist = hist[hist.index >= 1]                       # drop the 0 bin

growth = read_panacus(f'{COUNT}.growth', sep='\t', header=list(range(4)), index_col=0)
growth.columns = growth.columns.map(lambda x: (x[0], x[1], int(x[2]), float(x[3])))
growth = growth.reindex(sorted(growth.columns, key=lambda c: (c[3], c[2])), axis=1)
growth.index = growth.index.astype(int)
growth = growth[growth.index >= 1]                 # drop the 0 row

# -------------------------------------------------- data for panels D, E -----
if WITH_EVAL:
    raw = pd.read_csv(args.eval)
    raw = raw[raw['Region'] == 'Core Genome'].copy()
    raw['Approach'] = raw['Approach'].str.replace(' (', '\n(', n=1, regex=False)
    snv = raw[raw['Variant_type'] == 'Small Variants (<50bp)']
    sv = raw[raw['Variant_type'] == 'Structural Variants (>=50bp)']
    assert set(snv['Approach']) == set(SNV_LEVELS), set(snv['Approach']) ^ set(SNV_LEVELS)
    assert set(sv['Approach']) == set(SV_LEVELS), set(sv['Approach']) ^ set(SV_LEVELS)

    stats = pd.read_csv(args.stats)
    stats = stats[stats['Region'] == 'Core Genome']


def stars_for(linear, loov):
    """BH-adjusted paired Wilcoxon p for 'LOOV vs Linear' (same caller)."""
    key = f"{loov.replace(chr(10), ' ')} vs {linear.replace(chr(10), ' ')}"
    row = stats[stats['Comparison'] == key]
    if row.empty:
        return None
    p = float(row['wilcoxon_p_BH'].iloc[0])
    return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'


def bracket(ax, data, levels, linear, loov, gap, tick):
    label = stars_for(linear, loov)
    if label is None:
        return
    i, j = levels.index(linear), levels.index(loov)
    ymax = data[data['Approach'].isin([linear, loov])]['F_measure'].max()
    y = ymax + gap
    ax.plot([i, i, j, j], [y - tick, y, y, y - tick], color='#333333', lw=1.4,
            solid_capstyle='butt', zorder=5, clip_on=False)
    if label == 'ns':
        ax.text((i + j) / 2, y + tick * 0.5, 'ns', ha='center', va='bottom',
                fontsize=FONT['ns'], color='#111111')
    else:
        ax.text((i + j) / 2, y - tick * 0.2, label, ha='center', va='bottom',
                fontsize=FONT['stars'], color='#111111')


def eval_panel(ax, data, levels, title, ylim, yticks, brackets, bracket_gap, bracket_tick):
    n = len(levels)
    ax.axvspan(-0.5, n - 0.5, color=CORE_ACCENT, alpha=0.07, lw=0, zorder=0)
    groups = [data.loc[data['Approach'] == a, 'F_measure'].dropna().values for a in levels]
    bp = ax.boxplot(groups, positions=range(n), widths=0.6, patch_artist=True,
                    showfliers=False, whis=1.5, zorder=2,
                    medianprops=dict(color='#333333', lw=2.4),
                    boxprops=dict(edgecolor='#404040', lw=1.2),
                    whiskerprops=dict(color='#404040', lw=1.2),
                    capprops=dict(lw=0))
    for patch, a in zip(bp['boxes'], levels):
        patch.set_facecolor(APPROACH_COLOURS[a])
        patch.set_alpha(0.88)
    for k, g in enumerate(groups):          # one dot per strain, centred on the box
        ax.scatter(np.full(len(g), k), g, s=28, facecolor='white', edgecolor='#4D4D4D',
                   linewidth=0.8, alpha=0.9, zorder=3)
    ax.set_xlim(-0.6, n - 0.4)
    ax.set_ylim(*ylim)
    ax.set_yticks(yticks)
    ax.set_yticklabels([f'{t:.1f}' for t in yticks])
    ax.set_xticks(range(n))
    ax.set_xticklabels(levels, rotation=35, ha='right', rotation_mode='anchor',
                       fontsize=FONT['xtick_eval'], linespacing=0.95)
    ax.set_ylabel('F-measure', fontweight='bold')
    ax.set_title(title, loc='left', fontweight='bold', pad=10)
    ax.yaxis.grid(True, color='#EBEBEB', lw=0.9, zorder=1)
    ax.set_axisbelow(True)
    ax.xaxis.grid(False)
    for s in ax.spines.values():
        s.set_color('#CCCCCC')
        s.set_linewidth(1.0)
    ax.tick_params(colors='#333333', color='#CCCCCC', length=5)
    for (lin, loov) in brackets:
        bracket(ax, data, levels, lin, loov, bracket_gap, bracket_tick)


# ------------------------------------------------------------ the figure -----
LEFT, RIGHT = 0.060, 0.985
if WITH_EVAL:
    fig = plt.figure(figsize=(15, 12.2))
    gs_top = fig.add_gridspec(1, 3, left=LEFT, right=RIGHT, top=0.943, bottom=0.670, wspace=0.32)
    gs_bot = fig.add_gridspec(1, 2, left=LEFT, right=RIGHT, top=0.571, bottom=0.234,
                              width_ratios=[4, 8], wspace=0.14)
else:                                   # supplementary: A-C only, one short row
    fig = plt.figure(figsize=(15, 4.9))
    gs_top = fig.add_gridspec(1, 3, left=LEFT, right=RIGHT, top=0.855, bottom=0.215, wspace=0.32)

# ---- A: coverage histogram
ax = fig.add_subplot(gs_top[0])
ax.bar(range(len(hist)), hist[COUNT], color=GREEN, width=0.8)
int_xticks(ax, hist.index.astype(str))
style_y_human(ax)
ax.set_title('A: Coverage distribution', loc='left', fontweight='bold', pad=10)
ax.set_ylabel(UNIT)
ax.set_xlabel('coverage (# samples containing sequence)')

# ---- B: cumulative growth
ax = fig.add_subplot(gs_top[1])
for i, (t, ct, c, q) in enumerate(growth.columns):
    ax.bar(range(len(growth)), growth[(t, ct, c, q)], color=GROWTH_COLORS[i % 5], width=0.8,
           label=f'coverage ≥ {c}, quorum ≥ {q*100:.0f}%')
int_xticks(ax, growth.index.astype(str))
ax.set_ylim(0, growth.max().max() * 1.62)           # headroom: legend sits above the bars
style_y_human(ax)
ax.set_title('B: Cumulative growth', loc='left', fontweight='bold', pad=10)
ax.set_ylabel(UNIT)
ax.set_xlabel('genomes added')
ax.legend(loc='upper left', fontsize=FONT['legend_small'], frameon=True, framealpha=0.95,
          borderpad=0.5, labelspacing=0.3, handlelength=1.4)

# ---- C: novel sequence per added genome (coverage>=1, quorum>=0%)
ax = fig.add_subplot(gs_top[2])
cum = growth[growth.columns[0]]
marginal = cum.diff()
marginal.iloc[0] = cum.iloc[0]
ax.bar(range(len(marginal)), marginal, color=GREEN, width=0.8,
       label='coverage ≥ 1, quorum ≥ 0%')
int_xticks(ax, marginal.index.astype(str))
style_y_human(ax)
ax.set_title('C: Novel sequence per genome', loc='left', fontweight='bold', pad=10)
ax.set_ylabel(UNIT)
ax.set_xlabel('genomes added')
ax.legend(loc='upper right', fontsize=FONT['legend_small'], handlelength=1.4)

for a in fig.axes:                                   # light, consistent frame for A-C
    for s in a.spines.values():
        s.set_linewidth(1.0)

if WITH_EVAL:
    # ---- D: small variants, core genome
    axD = fig.add_subplot(gs_bot[0])
    eval_panel(axD, snv, SNV_LEVELS, 'D: Core genome, small variants (<50 bp)',
               ylim=(0.5, 1.08), yticks=np.arange(0.5, 1.01, 0.1),
               brackets=[("Linear Mapping\n(GATK)", "Pangenome Graph\n(LOOV: Surject + GATK)")],
               bracket_gap=0.028, bracket_tick=0.010)

    # ---- E: structural variants, core genome
    axE = fig.add_subplot(gs_bot[1])
    eval_panel(axE, sv, SV_LEVELS, 'E: Core genome, structural variants (≥50 bp)',
               ylim=(0.0, 1.0), yticks=np.arange(0.0, 1.01, 0.2),
               brackets=[("Linear Mapping\n(Delly)", "Pangenome Graph\n(LOOV: Delly Surject)"),
                         ("Linear Mapping\n(Manta)", "Pangenome Graph\n(LOOV: Manta Surject)"),
                         ("Linear Mapping\n(Dysgu)", "Pangenome Graph\n(LOOV: Dysgu Surject)")],
               bracket_gap=0.05, bracket_tick=0.016)

    # ---- shared legend for D and E
    handles = [Patch(facecolor=APPROACH_COLOURS[a], edgecolor='#404040', linewidth=0.8,
                     alpha=0.9, label=a) for a in LEGEND_ORDER]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, 0.012), ncol=5,
               fontsize=FONT['legend'], frameon=True, edgecolor='#DDDDDD',
               handlelength=1.1, handleheight=1.1, columnspacing=1.6, labelspacing=0.55,
               borderpad=0.7)

out = 'Figure1_bp' if WITH_EVAL else 'FigureS_node'
fig.savefig(f'{out}.tiff', dpi=300, facecolor='white', pil_kwargs={'compression': 'tiff_lzw'})
fig.savefig(f'{out}.png', dpi=200, facecolor='white')
fig.savefig(f'{out}.pdf', facecolor='white')
print(f'wrote {out}.tiff / .png / .pdf')
