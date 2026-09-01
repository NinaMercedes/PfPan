#!/usr/bin/env python3
from io import StringIO
from functools import partial
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

N_HEADERS = 4
CUSTOM_COLORS = ["#2BAE84", "#3366CC", "#8153A6", "#E87DBF", "#FF7033"]
N_SAMPLES = 13

def read_csv(path, comment='#', **kwargs):
    with open(path) as fh:
        lines = "".join(l for l in fh if not l.startswith(comment))
    return pd.read_csv(StringIO(lines), **kwargs)

def humanize_number(i, precision=0):
    order, x = 0, i
    if abs(i) > 0:
        order = int(np.log10(abs(i))) // 3
        x = i / 10**(order*3)
    human_r = ['', 'K', 'M', 'B', 'D']
    return '{:,.{prec}f}{:}'.format(x, human_r[order], prec=precision)

def calibrate_yticks_text(yticks):
    prec = 0
    txt = list(map(partial(humanize_number, precision=prec), yticks))
    while len(set(txt)) < len(txt):
        prec += 1
        txt = list(map(partial(humanize_number, precision=prec), yticks))
    return txt

# ---- load hist (final: -s paths.haplotypes.txt -S) ----
hist = read_csv('bp.hist', sep='\t', header=[0,1], index_col=0)
hist.columns = ['bp']
hist = hist[hist.index.notna()]
hist.index = hist.index.astype(int)
hist = hist.sort_index()
hist = hist[hist.index >= 1]   # drop the 0 column

# ---- load growth ----
growth = read_csv('bp.growth', sep='\t', header=list(range(N_HEADERS)), index_col=0)
growth.columns = growth.columns.map(lambda x: (x[0], x[1], int(x[2]), float(x[3])))
growth = growth.reindex(sorted(growth.columns, key=lambda c: (c[3], c[2])), axis=1)
growth.index = growth.index.astype(int)
growth = growth[growth.index >= 1]   # drop the 0 row

sns.set_theme(style='white')
plt.rcParams.update({'font.size': 13})

fig, axs = plt.subplots(1, 3, figsize=(21, 6))

# Panel A: coverage histogram
ax = axs[0]
ax.bar(hist.index.astype(str), hist['bp'], color=CUSTOM_COLORS[0])
ax.set_xticks(range(len(hist.index)))
ax.set_xticklabels(hist.index.astype(str), rotation=65, fontsize=10)
yticks = ax.get_yticks()
ax.set_yticks(yticks)
ax.set_yticklabels(calibrate_yticks_text(yticks))
ax.set_title('A: Node coverage distribution (bp)', fontsize=14, loc='left', fontweight='bold')
ax.set_ylabel('#bp')
ax.set_xlabel('coverage (# samples containing sequence)')

# Panel B: cumulative growth
ax = axs[1]
for i, (t, ct, c, q) in enumerate(growth.columns):
    label = f'coverage ≥ {c}, quorum ≥ {q*100:.0f}%'
    ax.bar(growth.index.astype(str), growth[(t, ct, c, q)], color=CUSTOM_COLORS[i % len(CUSTOM_COLORS)], label=label)
ax.set_xticks(range(len(growth.index)))
ax.set_xticklabels(growth.index.astype(str), rotation=65, fontsize=10)
yticks = ax.get_yticks()
ax.set_yticks(yticks)
ax.set_yticklabels(calibrate_yticks_text(yticks))
ax.set_title('B: Cumulative pangenome growth (bp)', fontsize=14, loc='left', fontweight='bold')
ax.set_ylabel('#bp')
ax.set_xlabel('genomes added')
ax.legend(loc='upper left', fontsize=9)

# Panel C: F_new (marginal growth), coverage>=1, quorum>=0%
ax = axs[2]
col = growth.columns[0]
cum = growth[col]
marginal = cum.diff()
marginal.iloc[0] = cum.iloc[0]
ax.bar(marginal.index.astype(str), marginal, color=CUSTOM_COLORS[0], label='coverage ≥ 1, quorum ≥ 0%')
ax.set_xticks(range(len(marginal.index)))
ax.set_xticklabels(marginal.index.astype(str), rotation=65, fontsize=10)
yticks = ax.get_yticks()
ax.set_yticks(yticks)
ax.set_yticklabels(calibrate_yticks_text(yticks))
ax.set_title('C: Novel sequence per additional genome (bp)', fontsize=14, loc='left', fontweight='bold')
ax.set_ylabel('#bp')
ax.set_xlabel('genomes added')
ax.legend(loc='upper right', fontsize=9)

plt.tight_layout()
plt.savefig('/mnt/user-data/outputs/Figure1_bp_FINAL.png', dpi=200)
plt.savefig('/mnt/user-data/outputs/Figure1_bp_FINAL.pdf')
print('done')
