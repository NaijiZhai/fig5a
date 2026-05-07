#!/usr/bin/env python3
from __future__ import annotations

import csv
import os
from pathlib import Path

import matplotlib.pyplot as plt

# By default, read WORK_DIR/profiles_simple. For a new run, just export WORK_DIR.
WORK_DIR = Path(os.environ.get('WORK_DIR', 'clean_pnas2017_fig5a'))
IN_DIR = WORK_DIR / 'profiles_simple'
# Three figure rows; file names match the outputs from make_profiles.py.
ROWS = [
    ('BPDE-dG tXR-seq', IN_DIR / 'BPDE_dG_GM12878_1h.profile.tsv'),
    ('CPD XR-seq', IN_DIR / 'CPD_NHF1_1h.profile.tsv'),
    ('(6-4)PP XR-seq', IN_DIR / 'PP64_NHF1_1h.profile.tsv'),
]
# Colors are chosen to roughly match the paper.
COLORS = {'NTS': '#c9474d', 'TS': '#6f9dc6'}

# Keep text editable in vector outputs and make exported figures cleaner.
plt.rcParams.update({
    'font.family': 'Arial',
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})


def read_profile(path, scale):
    # profile.tsv already has bin-level means; this only reads them and applies scaling.
    rows = []
    with open(path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            rows.append({
                'region': row['region'],
                'strand_class': row['strand_class'],
                'rel_mid': float(row['rel_mid']) / 1000.0,
                'rpkm': float(row['mean_RPKM']) * scale,
            })
    return rows


def draw(scale, suffix):
    # Share axes across panels to keep the layout close to the paper.
    fig, axes = plt.subplots(3, 2, figsize=(7.05, 6.9), sharex='col', sharey=True)
    all_data = [(label, read_profile(path, scale)) for label, path in ROWS]
    # Fix the y range so raw and paper_like figures are directly comparable.
    ymax = 3.5
    yticks = [0, 1, 2, 3]
    for i, (label, data) in enumerate(all_data):
        for j, region in enumerate(('TSS', 'TES')):
            ax = axes[i, j]
            # NTS is red and TS is blue. Sort by rel_mid before drawing the lines.
            for strand in ('NTS', 'TS'):
                sub = sorted([x for x in data if x['region'] == region and x['strand_class'] == strand],
                             key=lambda x: x['rel_mid'])
                ax.plot([x['rel_mid'] for x in sub], [x['rpkm'] for x in sub], color=COLORS[strand], linewidth=1.25,
                        label=strand)
            # Dashed line marks the TSS or TES anchor.
            ax.axvline(0, color='black', linestyle=(0, (5.5, 4.5)), linewidth=1.45)
            ax.grid(True, color='#e2e2e2', linewidth=1.0)
            ax.set_axisbelow(True)
            ax.set_ylim(0, ymax)
            ax.set_yticks(yticks)
            for spine in ax.spines.values():
                spine.set_color('#555555')
                spine.set_linewidth(0.95)
            ax.tick_params(axis='both', labelsize=10, width=1.0, length=3.5, colors='#444444')
        axes[i, 0].set_ylabel('RPKM', fontsize=13)
        axes[i, 1].text(1.06, 0.5, label, transform=axes[i, 1].transAxes, rotation=-90, va='center', ha='left',
                        fontsize=16)
    # Put the legend in the upper-right panel, roughly matching the paper layout.
    axes[0, 1].legend(frameon=False, loc='upper left', bbox_to_anchor=(0.17, 0.97), fontsize=14, handlelength=1.4)
    for ax in axes[:, 0]:
        ax.set_xlim(-5, 15)
        ax.set_xticks([-5, 0, 5, 10, 15])
    for ax in axes[:, 1]:
        ax.set_xlim(-15, 5)
        ax.set_xticks([-15, -10, -5, 0, 5])
    axes[2, 0].set_xticklabels(['-5', 'TSS', '5', '10', '15'], fontsize=13)
    axes[2, 1].set_xticklabels(['-15', '-10', '-5', 'TES', '5'], fontsize=13)
    fig.text(0.49, 0.045, 'Relative position (kb)', ha='center', va='center', fontsize=17)
    fig.text(0.005, 0.995, 'A', ha='left', va='top', fontsize=38, color='black')
    fig.subplots_adjust(left=0.125, right=0.845, bottom=0.13, top=0.985, wspace=0.12, hspace=0.15)
    fig.savefig(IN_DIR / f'fig5a_{suffix}.png', dpi=300)
    fig.savefig(IN_DIR / f'fig5a_{suffix}.pdf')
    plt.close(fig)


# raw is the BED-based RPKM. paper_like is a global 1.6x scale for visual comparison.
draw(1.0, 'raw')
# draw(1.6, 'paper_like_x1.6')
print(IN_DIR / 'fig5a_raw.png')
# print(IN_DIR / 'fig5a_paper_like_x1.6.png')
