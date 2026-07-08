#!/usr/bin/env python3
"""Collect Ed results from all configs and compute element-resolved statistics.

Usage: python3 collect_results.py [results_dir] ["elem1 elem2 ..."]
Output: Overall + per-element statistics printed to stdout.
"""
import sys
import os
import numpy as np

results_dir = sys.argv[1] if len(sys.argv) > 1 else "results"
elem_names = sys.argv[2].split() if len(sys.argv) > 2 else []

all_ed = []          # (ed_value, pka_type, config_id)
configs_found = set()

for entry in sorted(os.listdir(results_dir)):
    config_subdir = os.path.join(results_dir, entry)
    if not os.path.isdir(config_subdir) or not entry.startswith("config_"):
        continue
    config_id = entry.replace("config_", "")
    configs_found.add(config_id)

    for fname in sorted(os.listdir(config_subdir)):
        if not fname.startswith("Ed_direction_"):
            continue
        fpath = os.path.join(config_subdir, fname)
        with open(fpath) as fh:
            content = fh.read().strip()
        parts = content.split()
        try:
            ed = float(parts[0])
            pka_type = int(parts[1]) if len(parts) > 1 else 0
            all_ed.append((ed, pka_type, config_id))
        except (ValueError, IndexError):
            pass

if not all_ed:
    print("No valid Ed values found.")
    sys.exit(1)

ed_arr = np.array([v[0] for v in all_ed])
types_arr = np.array([v[1] for v in all_ed])

## Element name lookup
def type_name(t):
    if elem_names and 1 <= t <= len(elem_names):
        return elem_names[t - 1]
    return f"type_{t}"

## Overall statistics
n_valid = len(all_ed)
errors = 0  # no separate error count needed since we skip bad entries
mean_ed = np.mean(ed_arr)
std_ed = np.std(ed_arr)

print("=" * 55)
print("  Ed Summary Statistics")
print("=" * 55)
print(f"  Configs processed:    {len(configs_found)} ({sorted(configs_found)})")
print(f"  Valid directions:     {n_valid}")
print(f"  Mean Ed:   {mean_ed:.2f} eV")
print(f"  Std Ed:    {std_ed:.2f} eV")
print(f"  Min Ed:    {np.min(ed_arr):.2f} eV")
print(f"  Max Ed:    {np.max(ed_arr):.2f} eV")
print(f"  Median Ed: {np.median(ed_arr):.2f} eV")
print(f"  Q1 (25%):  {np.percentile(ed_arr, 25):.2f} eV")
print(f"  Q3 (75%):  {np.percentile(ed_arr, 75):.2f} eV")
print()

## Per-element statistics
unique_types = sorted(set(types_arr))
print(f"  Element-specific Ed (eV):")
print(f"  {'Element':<10s} {'mean':>8s} {'std':>8s} {'min':>8s} {'max':>8s} {'median':>8s} {'N':>5s}")
print(f"  {'-'*10} {'-'*8} {'-'*8} {'-'*8} {'-'*8} {'-'*8} {'-'*5}")
for t in unique_types:
    mask = types_arr == t
    vals = ed_arr[mask]
    tn = type_name(int(t))
    print(f"  {tn:<10s} {np.mean(vals):8.2f} {np.std(vals):8.2f} "
          f"{np.min(vals):8.2f} {np.max(vals):8.2f} "
          f"{np.median(vals):8.2f} {len(vals):5d}")
print(f"  {'overall':<10s} {mean_ed:8.2f} {std_ed:8.2f} "
      f"{np.min(ed_arr):8.2f} {np.max(ed_arr):8.2f} "
      f"{np.median(ed_arr):8.2f} {n_valid:5d}")
print()

## Per-config statistics
print(f"  Per-config Ed (eV):")
print(f"  {'Config':<10s} {'mean':>8s} {'std':>8s} {'min':>8s} {'max':>8s} {'N':>5s}")
print(f"  {'-'*10} {'-'*8} {'-'*8} {'-'*8} {'-'*8} {'-'*5}")
for cid in sorted(configs_found, key=int):
    mask = np.array([v[2] == cid for v in all_ed])
    vals = ed_arr[mask]
    if len(vals) > 0:
        print(f"  config_{cid:<4s} {np.mean(vals):8.2f} {np.std(vals):8.2f} "
              f"{np.min(vals):8.2f} {np.max(vals):8.2f} {len(vals):5d}")
print()

## Overall histogram
print(f"  Ed distribution (10 bins):")
hist, edges = np.histogram(ed_arr, bins=10)
max_h = max(hist) if max(hist) > 0 else 1
for i in range(len(hist)):
    bar = "#" * int(hist[i] * 40 / max_h)
    print(f"    {edges[i]:6.1f} - {edges[i+1]:6.1f} eV: {hist[i]:3d} {bar}")

print()
print(f"  All Ed values (eV):")
print(f"  [{', '.join(f'{v:.1f}' for v in ed_arr.tolist())}]")
print("=" * 55)
