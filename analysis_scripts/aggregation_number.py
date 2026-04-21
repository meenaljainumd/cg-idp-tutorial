"""
Weighted aggregation number (Nw) timeseries.
Computes the weighted aggregation number of peptides over the trajectory
using BB-BB contact-based graph clustering.
Usage:
    python aggregation_number.py <TPR> <XTC>
    python aggregation_number.py run.tpr run.xtc 
Requires general_metrics_functions.py in the same directory.
"""

import sys
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import MDAnalysis as mda
from tqdm import tqdm
from general_metrics_functions import agg_num

# ── Input ──────────────────────────────────────────────────────
tpr     = sys.argv[1]
xtc     = sys.argv[2]

# ── Load trajectory ────────────────────────────────────────────
u = mda.Universe(tpr, xtc)
pep_num = len(u.select_atoms("name BB").fragments)
print(f"Number of peptides: {pep_num}")

# ── Compute Nw per frame ───────────────────────────────────────
time_ns = []
nw_list = []

for ts in tqdm(u.trajectory, desc="Computing Nw"):
    time_ns.append(ts.time / 1000.0)

    agg_sizes = agg_num(u, ts, selection_string="name BB S1 S2 S3", cutoff=7)

    # Weighted aggregation number: Nw = sum(ni * ci) / sum(ci)
    # where ni = aggregate size, ci = number of peptides in that aggregate
    sizes = np.array(agg_sizes)
    nw = np.sum(sizes * sizes) / np.sum(sizes) if np.sum(sizes) > 0 else 1.0
    nw_list.append(nw)

time_ns = np.array(time_ns)
nw_arr  = np.array(nw_list)

# ── Plot ───────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(time_ns, nw_arr, linewidth=0.5, color="tab:purple")
ax.set_xlabel("Time (ns)", fontsize=12)
ax.set_ylabel("Weighted aggregation number (Nw)", fontsize=12)
ax.set_title("Peptide Aggregation Number", fontsize=12)
ax.set_ylim(bottom=0)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()

outpng = os.path.join( "aggregation_number.png")
plt.savefig(outpng, dpi=300, bbox_inches="tight")
print(f"Saved: {outpng}")

# ── Summary ────────────────────────────────────────────────────
print(f"\nMean Nw = {np.mean(nw_arr):.2f}")
print(f"Max  Nw = {np.max(nw_arr):.2f}")
