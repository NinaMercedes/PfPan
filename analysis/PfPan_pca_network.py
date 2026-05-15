import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
from matplotlib.patches import Ellipse
from sklearn.decomposition import PCA
import networkx as nx
from adjustText import adjust_text
import warnings
warnings.filterwarnings('ignore')

# ── Palettes ──────────────────────────────────────────────────────────────────
pals12 = ["#2BAE84","#3366CC","#8153A6","#E87DBF","#FF7033",
          "#F4A736","#D6A419","#3FB1C2","#7B5EA7","#A5426D","#6EC4E8","#8B5E3C"]

pals5 = ["#2BAE84","#3366CC","#8153A6","#E87DBF","#FF7033"]

type_order  = ['SNP','Small Insertion','Small Deletion',
               'Structural Insertion','Structural Deletion']
type_colors = {
    'SNP':                  pals5[0],
    'Small Insertion':      pals5[1],
    'Small Deletion':       pals5[2],
    'Structural Insertion': pals5[3],
    'Structural Deletion':  pals5[4],
}
sv_colors = {'Insertion': pals5[0], 'Deletion': pals5[4]}

sample_cols = ['Pf7G8','PfCD01','PfDd2','PfGA01','PfGB4','PfGN01','PfHB3',
               'PfIT','PfKE01','PfKH01','PfKH02','PfSN01']

# ── Load data ─────────────────────────────────────────────────────────────────
# variants = your dataframe
df = variants.copy()

# ── Top panels: variant summary data ─────────────────────────────────────────
vc      = df.copy()
sv_data = vc[vc['ABSLEN'] > 50].copy()
sv_data['SVTYPE2'] = np.where(sv_data['LEN'] > 0, 'Insertion', 'Deletion')

# ── SV matrix for PCA / network ───────────────────────────────────────────────
filt = df[(df['ABSLEN'] > 50) & (df['TYPE'] != 'SNP')].copy()
geno = filt[sample_cols].replace('.', np.nan).astype(float)
filt['n_alt'] = (geno == 1).sum(axis=1)
filt = filt[filt['n_alt'] >= 2].reset_index(drop=True)
geno = filt[sample_cols].replace('.', np.nan).astype(float)
geno['Pf3D7'] = 0.0
all_samples = sample_cols + ['Pf3D7']
n_all = len(all_samples)

mat_imp  = np.where(np.isnan(geno[all_samples].values), 0.5, geno[all_samples].values)
alt_bin  = (mat_imp >= 0.75).astype(float)
shared_b = alt_bin.T @ alt_bin
totals_b = alt_bin.sum(axis=0)

jaccard = np.zeros((n_all, n_all))
for i in range(n_all):
    for j in range(n_all):
        d = totals_b[i] + totals_b[j] - shared_b[i,j]
        jaccard[i,j] = shared_b[i,j] / d if d > 0 else 0

pca     = PCA(n_components=5)
X_pca   = pca.fit_transform(mat_imp.T)
var_exp = pca.explained_variance_ratio_ * 100

net_idx    = [all_samples.index(s) for s in sample_cols]
jac_net    = jaccard[np.ix_(net_idx, net_idx)]
totals_net = {s: int(totals_b[all_samples.index(s)]) for s in sample_cols}

sample_colors = {s: pals12[i] for i, s in enumerate(sample_cols)}
sample_colors['Pf3D7'] = '#999999'
sea_cluster = ['PfDd2','PfIT','PfKH01','PfKH02']

# ══════════════════════════════════════════════════════════════════════════════
# Figure: 2 rows — top = A B C, bottom = D E
# ══════════════════════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(22, 16), facecolor='white')
gs_outer = gridspec.GridSpec(
    2, 1, figure=fig,
    height_ratios=[1, 1.15],
    hspace=0.42,
    left=0.06, right=0.97, top=0.94, bottom=0.07
)
gs_top = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs_outer[0], wspace=0.38)
gs_bot = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_outer[1], wspace=0.30)

ax_a = fig.add_subplot(gs_top[0])
ax_b = fig.add_subplot(gs_top[1])
ax_c = fig.add_subplot(gs_top[2])
ax_d = fig.add_subplot(gs_bot[0])
ax_e = fig.add_subplot(gs_bot[1])

TITLE_FS = 13
LABEL_FS = 11
TICK_FS  = 10
BAR_FS   = 11

# ════════ Panel A — Variant types bar ════════════════════════════════════════
type_counts    = vc['TYPE'].value_counts()
present_types  = [t for t in type_order if t in type_counts.index]
counts         = [type_counts[t] for t in present_types]
cols_a         = [type_colors[t] for t in present_types]

bars = ax_a.bar(present_types, counts, color=cols_a,
                edgecolor='white', linewidth=0.5, width=0.6)
for bar, cnt in zip(bars, counts):
    ax_a.text(bar.get_x() + bar.get_width()/2,
              bar.get_height() + 30,
              f'{cnt:,}', ha='center', va='bottom',
              fontsize=BAR_FS, color='#333333')

ax_a.set_xlabel('Variant type', fontsize=LABEL_FS)
ax_a.set_ylabel('Count', fontsize=LABEL_FS)
ax_a.set_xticklabels(present_types, rotation=35, ha='right', fontsize=TICK_FS)
ax_a.tick_params(axis='y', labelsize=TICK_FS)
ax_a.set_ylim(0, max(counts) * 1.15)
for sp in ['top','right']:   ax_a.spines[sp].set_visible(False)
for sp in ['left','bottom']: ax_a.spines[sp].set_color('#CCCCCC')
ax_a.set_title('A   Variant types', fontsize=TITLE_FS,
               fontweight='normal', loc='left', pad=6)

# ════════ Panel B — SV length histogram ══════════════════════════════════════
bins     = np.linspace(50, sv_data['ABSLEN'].max(), 41)
ins_data = sv_data[sv_data['SVTYPE2'] == 'Insertion']['ABSLEN']
del_data = sv_data[sv_data['SVTYPE2'] == 'Deletion']['ABSLEN']

ax_b.hist(ins_data, bins=bins, color=sv_colors['Insertion'],
          edgecolor='white', linewidth=0.2, label='Insertion', alpha=0.9)
ax_b.hist(del_data, bins=bins, color=sv_colors['Deletion'],
          edgecolor='white', linewidth=0.2, label='Deletion',  alpha=0.9)

ax_b.set_xlabel('SV length (bp)', fontsize=LABEL_FS)
ax_b.set_ylabel('Count',          fontsize=LABEL_FS)
ax_b.tick_params(labelsize=TICK_FS)
ax_b.legend(fontsize=TICK_FS, frameon=False, loc='upper right')
for sp in ['top','right']:   ax_b.spines[sp].set_visible(False)
for sp in ['left','bottom']: ax_b.spines[sp].set_color('#CCCCCC')
ax_b.set_title('B   SV length distribution  (ABSLEN > 50 bp)',
               fontsize=TITLE_FS, fontweight='normal', loc='left', pad=6)

# ════════ Panel C — AF distribution ══════════════════════════════════════════
bins_af = np.arange(0, 1.05, 0.05)
for t in present_types:
    sub = vc[vc['TYPE'] == t]['AF']
    ax_c.hist(sub, bins=bins_af, color=type_colors[t],
              edgecolor='white', linewidth=0.2,
              label=t, alpha=0.88)

ax_c.set_xlabel('Allele frequency', fontsize=LABEL_FS)
ax_c.set_ylabel('Count',            fontsize=LABEL_FS)
ax_c.tick_params(labelsize=TICK_FS)
ax_c.legend(fontsize=8.5, frameon=False, loc='upper right',
            title='Type', title_fontsize=9)
for sp in ['top','right']:   ax_c.spines[sp].set_visible(False)
for sp in ['left','bottom']: ax_c.spines[sp].set_color('#CCCCCC')
ax_c.set_title('C   Allele frequency distribution',
               fontsize=TITLE_FS, fontweight='normal', loc='left', pad=6)

# ════════ Panel D — PCA ═══════════════════════════════════════════════════════
xs = X_pca[:, 0]
ys = X_pca[:, 1]

sea_idx_pca = [all_samples.index(s) for s in sea_cluster]
sea_pc = X_pca[sea_idx_pca, :2]
cx, cy = sea_pc[:,0].mean(), sea_pc[:,1].mean()
w = (sea_pc[:,0].max() - sea_pc[:,0].min()) * 2.0
h = (sea_pc[:,1].max() - sea_pc[:,1].min()) * 3.5
ax_d.add_patch(Ellipse((cx, cy), w, h, angle=10,
                        facecolor=pals12[4]+'22', edgecolor=pals12[4],
                        linewidth=1.8, linestyle='--', zorder=1))
ax_d.text(cx, cy - h/2 - 0.55, 'SE Asian cluster',
          ha='center', fontsize=10.5, color=pals12[4], fontstyle='italic')

for i, s in enumerate(all_samples):
    ax_d.scatter(xs[i], ys[i],
                 s=300 if s != 'Pf3D7' else 200,
                 color=sample_colors[s],
                 marker='D' if s == 'Pf3D7' else 'o',
                 edgecolors='#555555' if s == 'Pf3D7' else 'white',
                 linewidths=1.5, zorder=4)

texts = []
for i, s in enumerate(all_samples):
    t = ax_d.text(xs[i], ys[i], s,
                  fontsize=10.5,
                  color='#555555' if s == 'Pf3D7' else '#111111',
                  fontweight='bold' if s in sea_cluster else 'normal',
                  fontstyle='italic' if s == 'Pf3D7' else 'normal',
                  zorder=5)
    texts.append(t)

adjust_text(texts, x=xs, y=ys, ax=ax_d,
            expand_points=(2.5, 2.5), expand_text=(1.6, 1.6),
            arrowprops=dict(arrowstyle='-', color='#BBBBBB', lw=0.9),
            force_points=0.8, force_text=0.6, lim=500)

ax_d.axhline(0, color='#EBEBEB', lw=0.8, zorder=0)
ax_d.axvline(0, color='#EBEBEB', lw=0.8, zorder=0)
ax_d.set_xlabel(f'PC1  ({var_exp[0]:.1f}% variance explained)', fontsize=LABEL_FS)
ax_d.set_ylabel(f'PC2  ({var_exp[1]:.1f}% variance explained)', fontsize=LABEL_FS)
ax_d.tick_params(labelsize=TICK_FS)
for sp in ['top','right']:   ax_d.spines[sp].set_visible(False)
for sp in ['left','bottom']: ax_d.spines[sp].set_color('#CCCCCC')

ax_sc = ax_d.inset_axes([0.73, 0.73, 0.25, 0.24])
ax_sc.bar(range(1,6), var_exp[:5], color=pals12[2], alpha=0.85, width=0.6)
ax_sc.set_xticks(range(1,6))
ax_sc.set_xticklabels([f'PC{i}' for i in range(1,6)], fontsize=7.5)
ax_sc.set_ylabel('%', fontsize=7.5)
ax_sc.set_title('Variance\nexplained', fontsize=8, pad=2)
ax_sc.tick_params(labelsize=7.5)
for sp in ['top','right']:   ax_sc.spines[sp].set_visible(False)
for sp in ['left','bottom']: ax_sc.spines[sp].set_color('#CCCCCC')

ax_d.set_title('D   PCA — SV genotype matrix  (NA imputed as 0.5)',
               fontsize=TITLE_FS, fontweight='normal', loc='left', pad=8)

# ════════ Panel E — Network ═══════════════════════════════════════════════════
ax_e.set_facecolor('white')

G = nx.Graph()
for s in sample_cols:
    G.add_node(s)
for i, s1 in enumerate(sample_cols):
    for j, s2 in enumerate(sample_cols):
        if j <= i: continue
        G.add_edge(s1, s2, weight=jac_net[i,j])

pos   = nx.spring_layout(G, weight='weight', seed=42, k=3.2, iterations=400)
edges = list(G.edges())
wts   = [G[u][v]['weight'] for u,v in edges]
norm  = Normalize(vmin=min(wts), vmax=max(wts))

for (u,v), w in zip(edges, wts):
    nx.draw_networkx_edges(G, pos, edgelist=[(u,v)], ax=ax_e,
                           width=0.5 + norm(w)*6.5,
                           alpha=0.15 + norm(w)*0.72,
                           edge_color=[plt.cm.Purples(0.25 + norm(w)*0.65)])

node_list  = list(G.nodes())
node_sizes = [totals_net[s] for s in node_list]
node_cols  = [sample_colors[s] for s in node_list]
nx.draw_networkx_nodes(G, pos, nodelist=node_list, ax=ax_e,
                       node_size=node_sizes, node_color=node_cols,
                       edgecolors='white', linewidths=1.8)

px = np.array([pos[s][0] for s in node_list])
py = np.array([pos[s][1] for s in node_list])
ntexts = []
for s, x, y in zip(node_list, px, py):
    t = ax_e.text(x, y, s, fontsize=10.5, color='#111111',
                  fontweight='bold' if s in sea_cluster else 'normal',
                  ha='center', va='center', zorder=6)
    ntexts.append(t)

adjust_text(ntexts, x=px, y=py, ax=ax_e,
            expand_points=(2.8, 2.8), expand_text=(1.8, 1.8),
            arrowprops=dict(arrowstyle='-', color='#CCCCCC', lw=0.9),
            force_points=1.0, force_text=0.7, lim=500)

for sv_n, label in [(400,'400'),(500,'500'),(600,'600')]:
    ax_e.scatter([], [], s=sv_n, color='#CCCCCC',
                 edgecolors='white', linewidths=1.2, label=f'{label} SVs')
leg1 = ax_e.legend(title='Alt SV count', title_fontsize=10,
                    fontsize=9.5, frameon=True, loc='lower left',
                    framealpha=0.92, edgecolor='#DDDDDD')

edge_legend = [
    Line2D([0],[0], color=plt.cm.Purples(0.88), lw=5,   label='High similarity'),
    Line2D([0],[0], color=plt.cm.Purples(0.58), lw=2.8, label='Mid similarity'),
    Line2D([0],[0], color=plt.cm.Purples(0.32), lw=1.0, label='Low similarity'),
]
ax_e.legend(handles=edge_legend, title='Jaccard similarity',
            title_fontsize=10, fontsize=9.5, frameon=True,
            loc='lower right', framealpha=0.92, edgecolor='#DDDDDD')
ax_e.add_artist(leg1)
ax_e.axis('off')
ax_e.set_title('E   Haplotype similarity network',
               fontsize=TITLE_FS, fontweight='normal', loc='left', pad=8)

# ── Save ──────────────────────────────────────────────────────────────────────
fig.suptitle(
    'P. falciparum structural variant summary and haplotype relationships',
    fontsize=13, fontweight='normal', y=0.97, color='#222222'
)

plt.savefig('sv_summary_figure.png', dpi=200, bbox_inches='tight', facecolor='white')
plt.savefig('sv_summary_figure.tiff', dpi=300, bbox_inches='tight', facecolor='white')