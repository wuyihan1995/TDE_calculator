#!/usr/bin/env python3
"""
Plot Ed convergence as a function of number of simulated directions.

Overall + per-element cumulative mean ± std in a 1×2 subplot.
Saves plot_Ed.png at dpi=300.

Data: results/config_<N>/Ed_direction_<M>.txt  (<Ed> <PKA_type>)
"""

import os
import glob
import numpy as np
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── Academic plot defaults ──────────────────────────────────────────────────
rc = {
    'xtick.direction': 'inout',
    'ytick.direction': 'inout',
    'ytick.minor.visible': False,
    'xtick.minor.visible': False,
    'axes.linewidth': 1.5,
    'axes.labelweight': 'light',
    'xtick.major.width': 1.5,
    'ytick.major.width': 1.5,
    'xtick.minor.width': 1.5,
    'ytick.minor.width': 1.5,
    'xtick.major.size': 10.0,
    'ytick.major.size': 10.0,
    'xtick.minor.size': 6.0,
    'ytick.minor.size': 6.0,
    'axes.spines.left': True,
    'axes.spines.bottom': True,
    'axes.spines.right': False,
    'axes.spines.top': False,
    'xtick.bottom': True,
    'xtick.top': False,
    'ytick.left': True,
    'ytick.right': False,
    'font.family': 'sans-serif',
    'font.size': 18.0,
    'font.weight': 'light',
    'font.sans-serif': 'Arial',
    'axes.labelsize': 18.0,
    'axes.titlesize': 18.0,
    'xtick.labelsize': 18,
    'ytick.labelsize': 18,
    'axes.labelpad': 2.0,
    'legend.fontsize': 18,
    'legend.title_fontsize': 18.0,
    'legend.shadow': False,
    'legend.loc': 'best',
    'legend.labelspacing': 0.5,
    'legend.handleheight': 0.7,
    'legend.handlelength': 1.6,
    'legend.handletextpad': 0.5,
    'legend.framealpha': 0.6,
    'legend.frameon': False,
    'legend.fancybox': True,
    'legend.borderaxespad': 0.3,
    'legend.borderpad': 0.3,
    'legend.columnspacing': 1.0,
    'legend.edgecolor': '0.0',
    'legend.numpoints': 1,
    'legend.scatterpoints': 1,
    'legend.markerscale': 1.0,
    'axes.grid': False,
    'axes.grid.axis': 'both',
    'grid.linewidth': 1.2,
    'grid.color': '.6',
    'grid.linestyle': ':',
    'grid.alpha': 0.9,
    'text.color': '0.15',
    'xtick.color': '0.15',
    'ytick.color': '0.15',
    'figure.dpi': 200.0,
    'figure.subplot.wspace': 0.2,
    'figure.subplot.hspace': 0.2,
    'figure.autolayout': True,
    'figure.constrained_layout.use': False,
    'figure.constrained_layout.h_pad': 0.04167,
    'figure.constrained_layout.w_pad': 0.04167,
    'figure.constrained_layout.hspace': 0.02,
    'figure.constrained_layout.wspace': 0.02,
    'figure.titleweight': 'light',
    'figure.subplot.left': 0.125,
    'figure.subplot.right': 0.9,
    'figure.subplot.top': 0.88,
}

MyColors = ["#576fa0", "#9f9f9f", "#e3b87f", "#b57979", "#a7b9d7",
            "#91AD5A", "#DDB5E1", "#B7E2DE", "#C1BFE4"]
MyPalette = sns.color_palette(palette=MyColors, n_colors=len(MyColors), desat=None)
sns.set_theme(context='talk', style='ticks', palette=MyPalette,
              font='sans-serif', font_scale=1, color_codes=True, rc=rc)


# ── Load data ────────────────────────────────────────────────────────────────
data = []

results_dir = os.path.dirname(os.path.abspath(__file__))

for cdir in sorted(glob.glob(os.path.join(results_dir, "config_*"))):
    config_id = int(os.path.basename(cdir).split("_")[1])
    for fpath in sorted(glob.glob(os.path.join(cdir, "Ed_direction_*.txt"))):
        fname = os.path.basename(fpath)
        direction = int(fname.replace("Ed_direction_", "").replace(".txt", ""))
        with open(fpath) as fh:
            content = fh.read().strip()
        parts = content.split()
        if len(parts) < 2:
            continue
        try:
            ed = float(parts[0])
            pka_type = int(parts[1])
        except ValueError:
            continue
        data.append((ed, pka_type, config_id, direction))

if not data:
    raise SystemExit("No data found.")

# Sort by config then direction
data.sort(key=lambda x: (x[2], x[3]))
n_total = len(data)

eds = np.array([d[0] for d in data])
types = np.array([d[1] for d in data])

# ── Cumulative statistics ────────────────────────────────────────────────────
# Overall
overall_mean = np.zeros(n_total)
overall_std  = np.zeros(n_total)
for i in range(n_total):
    overall_mean[i] = np.mean(eds[:i+1])
    overall_std[i]  = np.std(eds[:i+1])

# Per element
unique_types = sorted(set(types))
elem_colors  = [MyPalette[i % len(MyPalette)] for i in range(len(unique_types))]
# Read element names from the parent directory's alloy.conf (if available)
elem_labels = {1: "Hf", 2: "Nb", 3: "Zr", 4: "Ti", 5: "Ta"}
conf_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "alloy.conf")
if os.path.exists(conf_path):
    with open(conf_path) as fh:
        for line in fh:
            if line.startswith("ELEMENT_NAMES="):
                raw = line.strip().split("=", 1)[1].strip()
                val = raw.strip('"')
                names = val.split()
                elem_labels = {i+1: name for i, name in enumerate(names)}
                break
elem_mean    = {t: np.full(n_total, np.nan) for t in unique_types}
elem_std     = {t: np.full(n_total, np.nan) for t in unique_types}
elem_seen    = {t: [] for t in unique_types}

for i in range(n_total):
    t = types[i]
    elem_seen[t].append(eds[i])
    for et in unique_types:
        if elem_seen[et]:
            arr = np.array(elem_seen[et])
            elem_mean[et][i] = np.mean(arr)
            elem_std[et][i]  = np.std(arr)

# ── Plot ─────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(13, 5.0))
x = np.arange(1, n_total + 1)

#  Left: overall
ax1 = plt.subplot(121)
ax1.plot(x, overall_mean, color=MyColors[0], linewidth=2)
ax1.fill_between(x, overall_mean - overall_std, overall_mean + overall_std,
                 color=MyColors[0], linewidth=0, alpha=0.2)
ax1.set_xlabel("Number of simulated directions")
ax1.set_ylabel("Threshold displacement energy (eV)", labelpad=8)
ax1.set_title("Overall", fontweight='light', pad=8)
mean_last10 = np.mean(overall_mean[-10:])
ax1.text(0.6, 0.90, f"mean (last 10): {mean_last10:.1f} eV",
         transform=ax1.transAxes, ha='center', fontsize=14,)
         #bbox=dict(boxstyle='round,pad=0.3', facecolor=None, edgecolor=None, alpha=0.8))

#  Right: per element
ax2 = plt.subplot(122)
lines = []
for i, t in enumerate(unique_types):
    l = ax2.plot(x, elem_mean[t], color=elem_colors[i], linewidth=1.5,
                 label=elem_labels.get(t, f"Type {t}"))
    lines.append(l[0])
    # ax2.fill_between(x, elem_mean[t] - elem_std[t], elem_mean[t] + elem_std[t],
    #                  color=elem_colors[i], linewidth=0, alpha=0.2)
ax2.set_xlabel("Number of simulated directions")
ax2.set_ylabel("Threshold displacement energy (eV)", labelpad=8)
ax2.set_title("Per element", fontweight='light', pad=8)
ax2.legend(lines, [f"{elem_labels.get(t, f'Type {t}')} ({elem_mean[t][-1]:.1f} eV)"
                   for t in unique_types],
           frameon=False, loc='upper left', bbox_to_anchor=(1, 0.95))

#  Common
for ax in [ax1, ax2]:
    ax.set_xlim(1, n_total)

plt.tight_layout()
out = os.path.join(results_dir, "plot_Ed.png")
plt.savefig(out, dpi=300, bbox_inches='tight')
print(f"Saved: {out}")
