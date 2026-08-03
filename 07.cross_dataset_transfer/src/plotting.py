"""
plotting.py -- Nature-style before/after distribution figures for harmonization.

Publication defaults: Arial/Helvetica, ~7pt, thin spines, no top/right axes,
vector-friendly. Two figures:
  * plot_before_after_grid: KDE of raw vs harmonized, datasets x {RNA, Protein}.
  * save_publication_figure: write PDF + PNG at 300 dpi.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


def set_nature_style():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 8,
        "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.dpi": 150, "savefig.dpi": 300, "pdf.fonttype": 42, "ps.fonttype": 42,
        "legend.frameon": False,
    })


def _kde(ax, vals, color, label, clip=(0.5, 99.5)):
    vals = np.asarray(vals, float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return
    if vals.size > 20000:
        vals = np.random.default_rng(0).choice(vals, 20000, replace=False)
    lo, hi = np.percentile(vals, clip)
    grid = np.linspace(lo, hi, 256)
    # simple Gaussian KDE without scipy dependency here
    from scipy.stats import gaussian_kde
    try:
        dens = gaussian_kde(vals)(grid)
        ax.fill_between(grid, dens, color=color, alpha=0.35, lw=0)
        ax.plot(grid, dens, color=color, lw=1.0, label=label)
    except Exception:
        ax.hist(vals, bins=60, density=True, color=color, alpha=0.4, label=label)


def plot_before_after_grid(
    raw: Dict[str, Dict[str, np.ndarray]],
    harmonized: Dict[str, Dict[str, np.ndarray]],
    ref: Dict[str, np.ndarray],
    datasets, title="Harmonization to CCLE reference",
):
    """
    raw/harmonized: {dataset: {"RNA": array, "Protein": array}}
    ref: {"RNA": ccle_rna_values, "Protein": ccle_protein_values}
    Renders a (len(datasets) x 4) grid: RNA raw | RNA harmonized | Prot raw | Prot harmonized.
    """
    set_nature_style()
    n = len(datasets)
    fig, axes = plt.subplots(n, 4, figsize=(7.2, 1.6 * n), constrained_layout=True)
    if n == 1:
        axes = axes[None, :]
    col_titles = ["RNA (raw)", "RNA (harmonized)", "Protein (raw)", "Protein (harmonized)"]
    blue, orange = "#3477eb", "#e8820c"
    for r, d in enumerate(datasets):
        panels = [
            (raw[d]["RNA"], ref["RNA"], "RNA"),
            (harmonized[d]["RNA"], ref["RNA"], "RNA"),
            (raw[d]["Protein"], ref["Protein"], "Protein"),
            (harmonized[d]["Protein"], ref["Protein"], "Protein"),
        ]
        for c, (vals, refvals, _t) in enumerate(panels):
            ax = axes[r, c]
            _kde(ax, refvals, orange, "CCLE ref")
            _kde(ax, vals, blue, d)
            if r == 0:
                ax.set_title(col_titles[c])
            if c == 0:
                ax.set_ylabel(d, fontweight="bold")
            ax.set_yticks([])
            if r == 0 and c == 3:
                ax.legend(loc="upper right")
    fig.suptitle(title, fontweight="bold")
    return fig


def save_publication_figure(fig, out_stem):
    out_stem = Path(out_stem)
    out_stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out_stem.with_suffix(".png"), bbox_inches="tight", dpi=300)
    return out_stem
