#!/usr/bin/env python3
"""
Figure 2: P. falciparum variant summary and haplotype relationships.

  A  Variant type counts
  B  Allele frequency distribution (stacked by type)
  C  Structural variant length (DEL / INS), coloured to match A and B
  D  PCA of SV genotypes. Colour = geographic region, shape = continent
  E  Jaccard similarity network, same colour/shape encoding as D

Filtering reproduces the R process_variants() step exactly:
  missing_frac <= 0.9, AF != 0, |LEN| < 10,000, then classify by LEN.

Usage:
    python PfPan_pca_network.py PfPan_all_variants_info_GT.tsv  [out_prefix]
"""

import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
from matplotlib.ticker import FuncFormatter
from sklearn.decomposition import PCA
import networkx as nx

TSV = sys.argv[1] if len(sys.argv) > 1 else 'PfPan_all_variants_info_GT.tsv'
OUT = sys.argv[2] if len(sys.argv) > 2 else 'fig2_sv_summary'

# ── Typography ────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family':      'DejaVu Sans',
    'font.size':        18,
    'axes.labelsize':   20,
    'axes.labelweight': 'bold',
    'xtick.labelsize':  17,
    'ytick.labelsize':  17,
    'legend.fontsize':  17,
    'legend.title_fontsize': 18,
    'axes.linewidth':   1.4,
    'xtick.major.width': 1.2,
    'ytick.major.width': 1.2,
    'xtick.major.size': 5,
    'ytick.major.size': 5,
})
PANEL_FS = 30
POINT_FS = 17

# ── Variant-type colours (unchanged from the original R figure) ───────────────
TYPE_ORDER  = ['SNP', 'Small Insertion', 'Small Deletion',
               'Structural Insertion', 'Structural Deletion']
TYPE_COLORS = {
    'SNP':                  '#2BAE84',
    'Small Insertion':      '#3366CC',
    'Small Deletion':       '#8153A6',
    'Structural Insertion': '#FF7033',
    'Structural Deletion':  '#E87DBF',
}

# ── Region colours: sampled from the population PCA figure ───────────────────
REGION_COLORS = {
    'West Africa':     '#E63946',
    'Central Africa':  '#2BAE84',
    'East Africa':     '#3366CC',
    'Southeast Asia':  '#F4A736',
    'South America':   '#E87DBF',
    'Central America': '#3FB1C2',   # HB3 only; not in the population panel
}
REGION_ORDER = list(REGION_COLORS)

# Shape encodes continent, so region is readable without colour
CONTINENT = {
    'West Africa': 'Africa', 'Central Africa': 'Africa', 'East Africa': 'Africa',
    'Southeast Asia': 'Asia',
    'South America': 'Americas', 'Central America': 'Americas',
}
CONT_MARKER = {'Africa': 'o', 'Asia': '^', 'Americas': 's'}

SAMPLE_REGION = {
    'Pf7G8':  'South America',    # Brazil
    'PfCD01': 'Central Africa',   # DRC
    'PfDd2':  'Southeast Asia',   # Indochina
    'PfGA01': 'West Africa',      # Gabon
    'PfGB4':  'West Africa',      # Ghana
    'PfGN01': 'West Africa',      # Guinea
    'PfHB3':  'Central America',  # Honduras
    'PfIT':   'South America',    # Brazil (Itajuba)
    'PfKE01': 'East Africa',      # Kenya
    'PfKH01': 'Southeast Asia',   # Cambodia
    'PfKH02': 'Southeast Asia',   # Cambodia
    'PfSN01': 'West Africa',      # Senegal
}
SAMPLES = list(SAMPLE_REGION)
REF_COLOR = '#8C8C8C'

def scol(s):
    return REF_COLOR if s == 'Pf3D7' else REGION_COLORS[SAMPLE_REGION[s]]

def smark(s):
    return 'D' if s == 'Pf3D7' else CONT_MARKER[CONTINENT[SAMPLE_REGION[s]]]

def short(s):
    return s if s == 'Pf3D7' else s[2:]   # drop the "Pf" prefix for labels

# ══════════════════════════════════════════════════════════════════════════════
# Data: replicate R process_variants()
# ══════════════════════════════════════════════════════════════════════════════
cols = ['CHR', 'POS', 'REF', 'ALT', 'QUAL', 'AC', 'AN', 'AF', 'NS',
        'SVTYPE', 'SVLEN'] + SAMPLES
df = pd.read_csv(TSV, sep='\t', header=None, names=cols, low_memory=False)
df['AF']     = pd.to_numeric(df['AF'], errors='coerce')
df['LEN']    = df['ALT'].str.len() - df['REF'].str.len()
df['ABSLEN'] = df['LEN'].abs()
miss = (df[SAMPLES].isin(['.']) | df[SAMPLES].isna()).sum(axis=1) / len(SAMPLES)
df = df[(miss <= 0.9) & (df['AF'] != 0) & (df['ABSLEN'] < 10000)].copy()

df['TYPE'] = np.select(
    [df['LEN'] == 0, df['LEN'] >= 50, df['LEN'] > 0, df['LEN'] <= -50],
    ['SNP', 'Structural Insertion', 'Small Insertion', 'Structural Deletion'],
    default='Small Deletion')

# SV genotype matrix for PCA / network (as in the original Python code)
filt = df[(df['ABSLEN'] > 50) & (df['TYPE'] != 'SNP')].copy()
geno = filt[SAMPLES].replace('.', np.nan).astype(float)
filt = filt[(geno == 1).sum(axis=1) >= 2]
geno = filt[SAMPLES].replace('.', np.nan).astype(float)
geno['Pf3D7'] = 0.0
ALL = SAMPLES + ['Pf3D7']

mat = np.where(np.isnan(geno[ALL].values), 0.5, geno[ALL].values)
alt = (mat >= 0.75).astype(float)
shared = alt.T @ alt
tot    = alt.sum(axis=0)
jac = np.divide(shared, tot[:, None] + tot[None, :] - shared,
                out=np.zeros_like(shared), where=(tot[:, None] + tot[None, :] - shared) > 0)

pca     = PCA(n_components=5)
X       = pca.fit_transform(mat.T)
var_exp = pca.explained_variance_ratio_ * 100

print(f'Variants after filtering: {len(df):,}')
print(f'SVs in PCA matrix:        {len(filt):,}')
print('Variance explained:', ', '.join(f'PC{i+1} {v:.1f}%' for i, v in enumerate(var_exp)))

# ══════════════════════════════════════════════════════════════════════════════
# Layout
# ══════════════════════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(21, 19), facecolor='white')
outer = gridspec.GridSpec(2, 1, figure=fig, height_ratios=[1, 1.3],
                          hspace=0.30, left=0.07, right=0.98, top=0.955, bottom=0.115)
top = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer[0],
                                       width_ratios=[1, 1.15, 1.15], wspace=0.34)
bot = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[1],
                                       width_ratios=[1.1, 1], wspace=0.12)
gs_c = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=top[2], hspace=0.12)

ax_a = fig.add_subplot(top[0])
ax_b = fig.add_subplot(top[1])
ax_cd = fig.add_subplot(gs_c[0])
ax_ci = fig.add_subplot(gs_c[1], sharex=ax_cd)
ax_d = fig.add_subplot(bot[0])
ax_e = fig.add_subplot(bot[1])

def clean(ax):
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)

def letter(ax, l, x=-0.16, y=1.04):
    ax.text(x, y, l, transform=ax.transAxes, fontsize=PANEL_FS,
            fontweight='bold', va='bottom', ha='left')

thousands = FuncFormatter(lambda v, _: f'{int(v):,}')

# ── A: variant type counts ────────────────────────────────────────────────────
counts = df['TYPE'].value_counts().reindex(TYPE_ORDER)
bars = ax_a.bar(range(5), counts.values, color=[TYPE_COLORS[t] for t in TYPE_ORDER],
                width=0.7, edgecolor='none')
for b, n in zip(bars, counts.values):
    ax_a.text(b.get_x() + b.get_width() / 2, n + counts.max() * 0.015, f'{n:,}',
              ha='center', va='bottom', fontsize=14, rotation=0)
ax_a.set_xticks(range(5))
ax_a.set_xticklabels(['SNP', 'Small\ninsertion', 'Small\ndeletion',
                      'Structural\ninsertion', 'Structural\ndeletion'],
                     rotation=45, ha='right', rotation_mode='anchor', fontsize=16)
ax_a.set_ylabel('Count')
ax_a.set_ylim(0, counts.max() * 1.12)
ax_a.yaxis.set_major_formatter(thousands)
clean(ax_a)
letter(ax_a, 'A', x=-0.30)

# ── B: AF distribution, stacked as in ggplot ──────────────────────────────────
bins_af = np.arange(0, 1.0 + 0.05, 0.05)
ax_b.hist([df.loc[df['TYPE'] == t, 'AF'] for t in TYPE_ORDER], bins=bins_af,
          stacked=True, color=[TYPE_COLORS[t] for t in TYPE_ORDER],
          edgecolor='white', linewidth=0.6, label=TYPE_ORDER)
ax_b.set_xlabel('Allele frequency')
ax_b.set_ylabel('Count')
ax_b.set_xlim(0, 1)
ax_b.yaxis.set_major_formatter(thousands)
ax_b.legend(frameon=False, loc='upper right', fontsize=16, handlelength=1.2,
            handleheight=1.0, title='Variant type', title_fontsize=17)
clean(ax_b)
letter(ax_b, 'B', x=-0.26)

# ── C: SV length, colours tied to panel A ─────────────────────────────────────
sv = df[df['TYPE'].isin(['Structural Insertion', 'Structural Deletion'])]  # >=50 bp, matches panel A
bins_log = np.logspace(np.log10(50), np.log10(10000), 41)
for ax, t, lab in [(ax_cd, 'Structural Deletion', 'Deletions'),
                   (ax_ci, 'Structural Insertion', 'Insertions')]:
    vals = sv.loc[sv['TYPE'] == t, 'ABSLEN']
    ax.hist(vals, bins=bins_log, color=TYPE_COLORS[t], edgecolor='white', linewidth=0.5)
    ax.set_xscale('log')
    ax.set_ylabel('Count')
    ax.text(0.97, 0.90, f'{lab} (n = {len(vals):,})', transform=ax.transAxes,
            ha='right', va='top', fontsize=17, fontweight='bold', color=TYPE_COLORS[t])
    clean(ax)
ax_ci.set_xlabel('SV length (bp)')
ax_ci.set_xticks([50, 100, 200, 500, 1000, 2000, 5000, 10000])
ax_ci.set_xticklabels(['50', '100', '200', '500', '1k', '2k', '5k', '10k'])
plt.setp(ax_cd.get_xticklabels(), visible=False)
letter(ax_cd, 'C', x=-0.26, y=1.08)

# ── D: PCA ────────────────────────────────────────────────────────────────────
# Label offsets in points (dx, dy); hand-set so the tight African cluster stays readable
LABEL_OFF = {
    'Pf7G8':  (14, 0),   'PfHB3':  (14, 0),
    'PfCD01': (14, 0),   'PfGB4':  (14, 0),
    'Pf3D7':  (-14, 4),  'PfGA01': (-14, -2),
    'PfKE01': (14, 4),   'PfGN01': (-14, -2),
    'PfSN01': (14, -4),
    'PfIT':   (-14, 8),  'PfKH02': (4, -20),
    'PfDd2':  (0, 16),   'PfKH01': (14, 0),
}
ax_d.axhline(0, color='#E3E3E3', lw=1, zorder=0)
ax_d.axvline(0, color='#E3E3E3', lw=1, zorder=0)
for i, s in enumerate(ALL):
    ax_d.scatter(X[i, 0], X[i, 1], s=380, marker=smark(s), color=scol(s),
                 edgecolors='white' if s != 'Pf3D7' else '#4D4D4D',
                 linewidths=1.6, zorder=4)
    dx, dy = LABEL_OFF[s]
    ax_d.annotate(short(s), (X[i, 0], X[i, 1]), xytext=(dx, dy),
                  textcoords='offset points', fontsize=POINT_FS, fontweight='bold',
                  color='#555555' if s == 'Pf3D7' else '#222222',
                  fontstyle='italic' if s == 'Pf3D7' else 'normal',
                  ha='left' if dx > 0 else ('right' if dx < 0 else 'center'),
                  va='center' if dy == 0 or abs(dx) > 0 else ('bottom' if dy > 0 else 'top'),
                  zorder=5)
ax_d.set_xlabel(f'PC1 ({var_exp[0]:.1f}%)')
ax_d.set_ylabel(f'PC2 ({var_exp[1]:.1f}%)')
pad_x = (X[:, 0].max() - X[:, 0].min()) * 0.14
pad_y = (X[:, 1].max() - X[:, 1].min()) * 0.10
ax_d.set_xlim(X[:, 0].min() - pad_x * 1.4, X[:, 0].max() + pad_x)
ax_d.set_ylim(X[:, 1].min() - pad_y * 1.3, X[:, 1].max() + pad_y)
clean(ax_d)
letter(ax_d, 'D', x=-0.12, y=1.02)

# scree inset in the empty upper-right area
ins = ax_d.inset_axes([0.60, 0.66, 0.36, 0.28])
ins.bar(range(1, 6), var_exp[:5], color='#9A9A9A', width=0.65)
ins.set_xticks(range(1, 6))
ins.set_xticklabels([f'PC{i}' for i in range(1, 6)], fontsize=14)
ins.tick_params(axis='y', labelsize=14)
ins.set_ylabel('Variance (%)', fontsize=15, fontweight='normal')
for sp in ('top', 'right'):
    ins.spines[sp].set_visible(False)

# ── E: similarity network ─────────────────────────────────────────────────────
idx = [ALL.index(s) for s in SAMPLES]
J = jac[np.ix_(idx, idx)]
G = nx.Graph()
G.add_nodes_from(SAMPLES)
for i, a in enumerate(SAMPLES):
    for j in range(i + 1, len(SAMPLES)):
        G.add_edge(a, SAMPLES[j], weight=J[i, j])

pos = nx.spring_layout(G, weight='weight', seed=42, k=3.2, iterations=400)
w = np.array([G[u][v]['weight'] for u, v in G.edges()])
norm = Normalize(vmin=w.min(), vmax=w.max())
for (u, v), wt in sorted(zip(G.edges(), w), key=lambda t: t[1]):
    f = norm(wt)
    ax_e.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
              color=plt.cm.Purples(0.30 + f * 0.65), lw=0.6 + f * 7.5,
              alpha=0.20 + f * 0.75, solid_capstyle='round', zorder=1)

NODE_SCALE = 1.6
sv_count = {s: tot[ALL.index(s)] for s in SAMPLES}
for s in SAMPLES:
    ax_e.scatter(*pos[s], s=sv_count[s] * NODE_SCALE, marker=smark(s), color=scol(s),
                 edgecolors='white', linewidths=2, zorder=3)

# labels pushed radially outward from the layout centre so they never sit on a node
centre = np.mean([pos[s] for s in SAMPLES], axis=0)
for s in SAMPLES:
    v = np.array(pos[s]) - centre
    v = v / (np.linalg.norm(v) + 1e-9)
    r = np.sqrt(sv_count[s] * NODE_SCALE) / 2 + 10
    ax_e.annotate(short(s), pos[s], xytext=(v[0] * r, v[1] * r), textcoords='offset points',
                  ha='left' if v[0] > 0.3 else ('right' if v[0] < -0.3 else 'center'),
                  va='bottom' if v[1] > 0.3 else ('top' if v[1] < -0.3 else 'center'),
                  fontsize=POINT_FS, fontweight='bold', color='#222222', zorder=5)

ax_e.set_aspect('equal', adjustable='datalim')
ax_e.margins(0.16)
ax_e.axis('off')
letter(ax_e, 'E', x=-0.02, y=1.02)

lo, hi = int(tot[idx].min() // 50 * 50), int(np.ceil(tot[idx].max() / 50) * 50)
size_h = [Line2D([0], [0], marker='o', ls='', color='#BDBDBD', markeredgecolor='white',
                 markersize=np.sqrt(n * NODE_SCALE), label=f'{n:,}')
          for n in (lo, (lo + hi) // 2, hi)]
leg_s = ax_e.legend(handles=size_h, title='Alt SVs', loc='lower left',
                    bbox_to_anchor=(-0.04, -0.10), frameon=False, labelspacing=1.3,
                    borderpad=0.2, handletextpad=1.0, fontsize=16)
edge_h = [Line2D([0], [0], color=plt.cm.Purples(0.30 + f * 0.65), lw=0.6 + f * 7.5,
                 label=f'{w.min() + f * (w.max() - w.min()):.2f}') for f in (1, 0.5, 0.05)]
ax_e.legend(handles=edge_h, title='Jaccard similarity', loc='lower right',
            bbox_to_anchor=(1.07, -0.10), frameon=False, fontsize=16)
ax_e.add_artist(leg_s)

# ── Shared geographic legend under D and E ────────────────────────────────────
geo_h = [Line2D([0], [0], marker=CONT_MARKER[CONTINENT[r]], ls='', markersize=17,
                 color=REGION_COLORS[r], markeredgecolor='white', label=r)
         for r in REGION_ORDER]
geo_h.append(Line2D([0], [0], marker='D', ls='', markersize=15, color=REF_COLOR,
                    markeredgecolor='#4D4D4D', label='Pf3D7 (reference)'))
fig.legend(handles=geo_h, loc='lower center', bbox_to_anchor=(0.53, 0.005),
           ncol=7, frameon=False, fontsize=18, title='Geographic origin',
           title_fontsize=19, columnspacing=1.2, handletextpad=0.4)

fig.savefig(f'{OUT}.png',  dpi=200, facecolor='white')
fig.savefig(f'{OUT}.tiff', dpi=300, facecolor='white', pil_kwargs={'compression': 'tiff_lzw'})
fig.savefig(f'{OUT}.pdf',  facecolor='white')
print(f'Saved {OUT}.png / .tiff / .pdf')
