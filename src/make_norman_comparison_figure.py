#!/usr/bin/env python3
"""Render complete common-inference counts and paired Cycling score directions."""
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def main():
    summary = json.loads((ROOT / 'data_norman_continuation/shared_inference/summary.json').read_text())
    arms = ['all_cell_sae_seed0', 'all_cell_sae_seed1', 'all_cell_sae_seed2',
            'control_only_sae_seed0', 'direct_pca_seed0', 'capped_direct_nmf_seed0']
    panels = list(summary['programs'])
    counts = np.array([[summary['programs'][p][a]['n_associations'] if
                        summary['programs'][p][a]['matched'] else np.nan for a in arms] for p in panels])
    with (ROOT / 'data_norman_continuation/shared_inference/comparisons.csv').open() as stream:
        rows = [r for r in csv.DictReader(stream) if r['arm'] == 'direct_pca_seed0' and r['program'] == 'Cycling']
    assert len(rows) == 105
    y = np.array([float(r['auc']) for r in rows])
    s = np.array([float(r['S_auc']) for r in rows])
    g = np.array([float(r['G2M_auc']) for r in rows])
    both = (y >= .55) & (s <= .45) & (g <= .45)
    plt.rcParams.update({'font.size': 9, 'svg.hashsalt': 'norman-six-model-comparison-v1',
                         'font.family': 'DejaVu Sans'})
    fig = plt.figure(figsize=(13, 5.5), layout='constrained')
    grid = fig.add_gridspec(1, 3, width_ratios=[1.35, 1, 1])
    ax = fig.add_subplot(grid[0])
    cmap = plt.get_cmap('Blues').copy()
    cmap.set_bad('#e1e1e1')
    ax.imshow(np.ma.masked_invalid(counts), cmap=cmap, vmin=0, vmax=54, aspect='auto')
    ax.set_xticks(range(6), ['SAE 0', 'SAE 1', 'SAE 2', 'Control\nSAE 0', 'PCA', 'NMF*'], rotation=40, ha='right')
    ax.set_yticks(range(9), [p.replace('_', ' ') for p in panels])
    for i in range(9):
        for j in range(6):
            value = counts[i, j]
            ax.text(j, i, 'U' if np.isnan(value) else str(int(value)), ha='center', va='center',
                    color='white' if value > 32 else '#1a1a1a', fontsize=9)
    ax.set_title('A  Joint association counts', loc='left', fontweight='bold')
    ax.set_xlabel('U: unscored   *NMF reached cap without convergence', fontsize=8)
    for index, (values, name) in enumerate([(s, 'S'), (g, 'G2M')], start=1):
        axis = fig.add_subplot(grid[index])
        axis.axhspan(.55, .95, xmin=0, xmax=.45, color='#f7eee8', alpha=.7)
        axis.scatter(values[~both], y[~both], s=23, c='#457b9d', alpha=.8, linewidths=0)
        axis.scatter(values[both], y[both], s=30, c='#b44530', alpha=.9, linewidths=0,
                     label=f'Both phase AUCs ≤ .45\nwith PCA AUC ≥ .55 (n={both.sum()})')
        axis.axhline(.55, color='#555555', linestyle='--', linewidth=.8)
        axis.axvline(.45, color='#555555', linestyle='--', linewidth=.8)
        axis.axhline(.5, color='#999999', linewidth=.6)
        axis.axvline(.5, color='#999999', linewidth=.6)
        axis.set(xlim=(0, 1), ylim=(.25, .9), xlabel=f'External {name} score AUC', ylabel='Cycling PCA feature AUC')
        axis.set_title(f'{"B" if index == 1 else "C"}  Paired directions: {name}', loc='left', fontweight='bold')
        axis.legend(loc='lower left', fontsize=7, frameon=False)
        axis.spines[['top', 'right']].set_visible(False)
    fig.suptitle('Common inference retains different measurements under the same program label', fontsize=12)
    dest = ROOT / 'docs/figures'
    for suffix in ['png', 'svg']:
        metadata = {'Date': None} if suffix == 'svg' else {'Software': 'SAE-scRNA-seq-Research'}
        fig.savefig(dest / f'fig_norman_common_inference.{suffix}', dpi=200, metadata=metadata)
    plt.close(fig)


if __name__ == '__main__':
    main()
