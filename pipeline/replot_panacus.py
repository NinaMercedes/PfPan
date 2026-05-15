#!/usr/bin/env python3

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from sys import stdout, stderr, exit
from functools import partial
from io import StringIO
from os import fdopen
import re

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from sklearn.linear_model import LinearRegression
from scipy.optimize import curve_fit
import seaborn as sns
from matplotlib.transforms import Bbox

PAT_PANACUS = re.compile(r'^#.+panacus (\S+) (.+)')
N_HEADERS = 4
SUPPORTED_FILE_FORMATS = plt.gcf().canvas.get_supported_filetypes().keys()

ids = pd.IndexSlice

class colors:
    WARNING = '\033[93m'
    END_COLOR = '\033[0m'

# custom palette
CUSTOM_COLORS = [
    "#2BAE84", "#3366CC", "#8153A6", "#E87DBF", "#FF7033",
    "#F4A736", "#D6A419", "#3FB1C2", "#8ACB4A", "#A5426D", "#6EC4E8"
]


def read_csv(filename, comment='#', **kwargs):
    lines = "".join([line for line in filename if not line.startswith(comment)])
    return pd.read_csv(StringIO(lines), **kwargs)

def clean_multicolumn_labels(df):
    column_header = []
    for c in df.columns:
        if isinstance(c, tuple):
            c = tuple((not x.startswith('Unnamed:') and x or '' for x in c))
        elif isinstance(c, str) and c.startswith('Unnamed:'):
            c = ''
        column_header.append(c)
    df.columns = pd.Index(column_header, tupleize_cols=True)
    return df

def humanize_number(i, precision=0):
    order = 0
    x = i
    if abs(i) > 0:
        order = int(np.log10(abs(i)))//3
        x = i/10**(order*3)
    human_r= ['', 'K', 'M', 'B', 'D']
    return '{:,.{prec}f}{:}'.format(x, human_r[order], prec=precision)

def calibrate_yticks_text(yticks):
    prec = 0
    yticks_text = list(map(partial(humanize_number, precision=prec), yticks))
    while len(set(yticks_text)) < len(yticks_text):
        prec += 1
        yticks_text = list(map(partial(humanize_number, precision=prec), yticks))
    return yticks_text

def fit(X, Y, func):
    popt, pcov = curve_fit(func, X, Y, p0=[1, 1], maxfev=1000*len(Y))
    return popt, func(X, *popt)

def fit_gamma(Y):
    X = np.arange(len(Y))+1
    return fit(X, Y, lambda x, *y: y[0]*x**y[1])

def fit_alpha(Y):
    X = np.arange(len(Y))+2
    return fit(X, Y, lambda x, *y:  y[0]*x**(-y[1]))

def plot_hist(df, ax, loc='lower left'):
    df.plot.bar(ax=ax, label=df.columns[0][1], color=CUSTOM_COLORS[0])
    ax.set_xticks(ax.get_xticks())
    ax.set_xticklabels(ax.get_xticks(), rotation=65, fontsize=20)
    yticks = ax.get_yticks()
    ax.set_yticks(yticks)
    ax.set_yticklabels(calibrate_yticks_text(yticks), fontsize=20)
    ax.set_title(f'Coverage histogram for #{df.columns[0][1]}s', fontsize=20)
    ax.set_ylabel(f'#{df.columns[0][1]}s', fontsize=20)
    ax.get_legend().remove()

def plot_growth(df, axs, loc='lower left', estimate_growth=False):
    popts = []
    df = df.reindex(sorted(df.columns, key=lambda x: (x[3], x[2])), axis=1)
    for i, (t, ct, c, q) in enumerate(df.columns):
        label = f'coverage ≥ {c}, quorum ≥ {q*100:.0f}%'
        df[(t, ct, c, q)].plot.bar(
            color=CUSTOM_COLORS[i % len(CUSTOM_COLORS)],
            label=label,
            ax=axs[0]
        )
        if c <= 1 and q <= 1/df.shape[0]:
            if estimate_growth:
                popt, curve = fit_gamma(df.loc[1:, (t, ct, c,q)].array)
                popts.append((c, q, popt, i))
                axs[0].plot(df.loc[1:].index, curve, '--', color='black',
                            label=f'coverage ≥ {c}, quorum ≥ {q*100:.0f}%, k₁={humanize_number(popt[0],1)}, γ={popt[1]:.3f}')
            else:
                popts.append((c, q, None, i))
    axs[0].set_xticklabels(axs[0].get_xticklabels(), rotation=65, fontsize=20)
    yticks = axs[0].get_yticks()
    axs[0].set_yticks(yticks)
    axs[0].set_yticklabels(calibrate_yticks_text(yticks), fontsize=20)
    axs[0].set_title(f'{df.columns[0][0]} plot for #{df.columns[0][1]}s', fontsize=20)
    axs[0].set_ylabel(f'#{df.columns[0][1]}s', fontsize=20)
    axs[0].set_xlabel('samples', fontsize=20)
    axs[0].legend(loc=loc, fontsize=16)

    if popts:
        for c, q, _, i in popts:
            x = np.zeros(df.shape[0]-1)
            x[1:] = df.loc[df.index[1]:df.index[-2], (t, ct, c, q)]
            label = f'coverage ≥ {c}, quorum ≥ {q*100:.0f}%'
            (df.loc[df.index[1]:, (t, ct, c, q)] - x).plot.bar(
                color=CUSTOM_COLORS[i % len(CUSTOM_COLORS)],
                label=label,
                ax=axs[1]
            )
            if estimate_growth:
                popt, _ = fit_alpha((df.loc[df.index[2]:, (t, ct, c, q)] - x[1:]).array)
                k2 = popt[0]
                alpha = popt[1]
                Y = k2*np.arange(1, df.shape[0]+1)**(-alpha)
                axs[1].plot(Y, '--', color='black',
                            label=f'coverage ≥ {c}, quorum ≥ {q*100:.0f}%, k₂={humanize_number(k2,1)}, α={alpha:.3f}')
        axs[1].set_xticklabels(axs[1].get_xticklabels(), rotation=65, fontsize=20)
        yticks = axs[1].get_yticks()
        axs[1].set_yticks(yticks)
        axs[1].set_yticklabels(calibrate_yticks_text(yticks), fontsize=20)
        axs[1].set_title(f'F_new (#{df.columns[0][1]}s)', fontsize=20)
        axs[1].set_ylabel(f'#{df.columns[0][1]}s', fontsize=20)
        axs[1].set_xlabel('samples', fontsize=20)
        axs[1].legend(loc=loc, fontsize=16)

def get_subplot_dim(df):
    growths = [x for x in df.columns.levels[0] if x.endswith('growth')]
    non_cum = 0
    if growths:
        non_cum = df.loc[:, ids[growths, :, :, :]].columns.map(
            lambda c: c[2] <= 1 and c[3] <= 1/df.shape[0]
        ).any() and len(growths) or 0
    return len(df.columns.levels[1]), len(df.columns.levels[0]) + non_cum, non_cum

def full_extent(ax, pad=0.0):
    ax.figure.canvas.draw()
    items = ax.get_xticklabels() + ax.get_yticklabels()
    items += [ax, ax.title, ax.xaxis.label, ax.yaxis.label]
    items += [ax, ax.title]
    bbox = Bbox.union([item.get_window_extent() for item in items])
    return bbox.expanded(1.0 + pad, 1.0 + pad)

def save_split_figures(axs, f, format, prefix):
    for i, ax_row in enumerate(axs):
        for j, ax in enumerate(ax_row):
            extent = full_extent(ax).transformed(f.dpi_scale_trans.inverted())
            with open(f'{prefix}{i}_{j}.{format}', 'wb+') as out:
                plt.savefig(out, bbox_inches=extent, format=format)

if __name__ == '__main__':
    print(colors.WARNING + "WARNING: panacus-visualize is DEPRECATED and will soon be REMOVED. Please use the report command of panacus to generate plots." + colors.END_COLOR, file=stderr)

    parser = ArgumentParser(formatter_class=ADHF, description='Visualize growth stats.')
    parser.add_argument('stats', type=open, help='Growth/Histogram table computed by panacus')
    parser.add_argument('-e', '--estimate_growth_params', action='store_true', help='Estimate growth parameters')
    parser.add_argument('-l', '--legend_location', choices=['lower left','lower right','upper left','upper right'], default='upper left', help='Legend location')
    parser.add_argument('-s', '--figsize', nargs=2, type=int, default=[10, 6], help='Figure size')
    parser.add_argument('-f', '--format', default='pdf' in SUPPORTED_FILE_FORMATS and 'pdf' or list(SUPPORTED_FILE_FORMATS)[0], choices=SUPPORTED_FILE_FORMATS, help='Output format')
    parser.add_argument('--split_subfigures', action='store_true', help='Split into multiple files')
    parser.add_argument('--split_prefix', default='out_', help='Prefix for split files')

    args = parser.parse_args()

    df = clean_multicolumn_labels(read_csv(args.stats, sep='\t', header=list(range(N_HEADERS)), index_col=[0], comment='#'))
    if df.columns[0][0] not in ['hist', 'growth', 'ordered-growth']:
        print('Unsupported table type, exiting.', file=stderr)
        exit(1)
    df.columns = df.columns.map(lambda x: (x[0], x[1], x[2] and int(x[2].replace("A", "")), x[3] and float(x[3].replace("R", ""))))

    n, m, non_cum_plots = get_subplot_dim(df)

    # clean white style
    sns.set_theme(style='white')
    plt.rcParams.update({'font.size': 20})

    f, axs = plt.subplots(n, m, figsize=(args.figsize[0] * m * 1.25 , args.figsize[1] * n * 2 ))

    if m == 1 and n == 1:
        axs = np.array([[axs]])
    elif m == 1:
        axs = axs.reshape(axs.size, 1)
    elif n == 1:
        axs = axs.reshape(1, axs.size)

    for t in df.columns.levels[0]:
        for j, c in enumerate(df.columns.levels[1]):
            df_tc = df.loc[:, ids[t, c, :, :]]
            if t == 'hist':
                plot_hist(df_tc, axs[j, 0], loc=args.legend_location)
            elif t == 'growth':
                offset = 'hist' in df.columns.levels[0] and 1 or 0
                axs_tc = axs[j, offset:offset+1]
                if non_cum_plots:
                    axs_tc = axs[j, offset:offset+2]
                plot_growth(df_tc, axs_tc, loc=args.legend_location, estimate_growth=args.estimate_growth_params)
            elif t == 'ordered-growth':
                if args.estimate_growth_params:
                    print('Cannot estimate growth with ordered-growth', file=stderr)
                    exit(1)
                axs_tc = axs[j, -1:]
                if non_cum_plots:
                    axs_tc = axs[j, -2:]
                if df_tc.index[0] == '0' and df_tc.loc['0'].isna().all():
                    df_tc.drop(['0'], inplace=True)
                plot_growth(df_tc, axs_tc, loc=args.legend_location, estimate_growth=False)

    plt.tight_layout()
    if not args.split_subfigures:
        with fdopen(stdout.fileno(), 'wb', closefd=False) as out:
            plt.savefig(out, format=args.format)
    else:
        save_split_figures(axs, f, args.format, args.split_prefix)

    plt.close()
