"""
Beta-sheet fraction timeseries.

Computes the fraction of peptide pairs forming contiguous beta-sheet
contacts (≥4 consecutive BB-BB contacts)
over the trajectory.

Usage:
    python beta_sheet_fraction.py <TPR> <XTC> <CSV_FOLDER> <IMAGE_FOLDER>

Example:
    python beta_sheet_fraction.py run.tpr run.xtc results/ figures/

Requires general_metrics_functions.py in the same directory.
"""

import sys
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
from tqdm import tqdm
import networkx as nx
from general_metrics_functions import contiguous_beta_sheet

# ── Input ──────────────────────────────────────────────────────
if len(sys.argv) != 5:
    print("Usage: python beta_sheet_fraction.py <TPR> <XTC> <CSV_FOLDER> <IMAGE_FOLDER>")
    sys.exit(1)

tpr     = sys.argv[1]
xtc     = sys.argv[2]
csv_dir = sys.argv[3]
img_dir = sys.argv[4]

os.makedirs(csv_dir, exist_ok=True)
os.makedirs(img_dir, exist_ok=True)

# ── Parameters ─────────────────────────────────────────────────
BB_BB_CUTOFF = 7.0        # Angstrom cutoff for BB-BB contacts
MIN_CONTIG   = 4          # minimum contiguous contacts for beta-sheet

# ── Load trajectory ────────────────────────────────────────────
u = mda.Universe(tpr, xtc)
bb = u.select_atoms("name BB")
pep_num = len(bb.fragments)
beads_per_pep = len(bb.fragments[0].select_atoms("name BB"))

print(f"Number of peptides: {pep_num}")
print(f"BB beads per peptide: {beads_per_pep}")

# ── Compute beta-sheet fraction per frame ──────────────────────
time_ns = []
beta_fraction = []       # fraction of peptides in a beta-sheet
n_beta_pairs = []        # number of peptide pairs with beta contacts
largest_sheet = []       # size of the largest beta-sheet cluster

for ts in tqdm(u.trajectory, desc="Computing beta-sheet"):
    time_ns.append(ts.time / 1000.0)

    # Distance array of all BB beads
    bb_bb = distance_array(bb.positions, bb.positions, ts.dimensions)

    # Remove intrapeptide contacts and lower triangle
    for pep_id in range(pep_num):
        bb_bb[
            pep_id * beads_per_pep : pep_id * beads_per_pep + beads_per_pep,
            : pep_id * beads_per_pep + beads_per_pep
        ] += BB_BB_CUTOFF + 10

    # Boolean contact matrix
    bb_bb_bool = (bb_bb <= BB_BB_CUTOFF).astype(int)

    # Check each peptide pair for contiguous beta-sheet contacts
    beta_sheet_array = np.zeros((pep_num, pep_num), dtype=int)

    for pep_i in range(pep_num):
        for pep_j in range(pep_i + 1, pep_num):
            ij_contacts = bb_bb_bool[
                pep_i * beads_per_pep : (pep_i + 1) * beads_per_pep,
                pep_j * beads_per_pep : (pep_j + 1) * beads_per_pep,
            ]
            has_beta, _, _ = contiguous_beta_sheet(
                ij_contacts, min_contiguous_contacts=MIN_CONTIG
            )
            beta_sheet_array[pep_i, pep_j] = int(has_beta)

    # Count pairs with beta contacts
    n_pairs = int(np.sum(beta_sheet_array))
    n_beta_pairs.append(n_pairs)

    # Build graph to find connected beta-sheet clusters
    G = nx.Graph(beta_sheet_array)
    beta_components = [
        list(c) for c in sorted(nx.connected_components(G), key=len, reverse=True)
    ]

    # Peptides in beta-sheets = those in clusters of size >= 2
    peps_in_beta = sum(len(c) for c in beta_components if len(c) >= 2)
    beta_fraction.append(peps_in_beta / pep_num)

    # Largest beta-sheet cluster size
    largest = max(len(c) for c in beta_components) if beta_components else 0
    largest_sheet.append(largest if largest >= 2 else 0)

time_ns       = np.array(time_ns)
beta_fraction = np.array(beta_fraction)
n_beta_pairs  = np.array(n_beta_pairs)
largest_sheet = np.array(largest_sheet)

# ── Save CSV ───────────────────────────────────────────────────
outcsv = os.path.join(csv_dir, "beta_sheet_fraction.csv")
np.savetxt(outcsv,
           np.column_stack([time_ns, beta_fraction, n_beta_pairs, largest_sheet]),
           header="time_ns,beta_fraction,n_beta_pairs,largest_sheet_size",
           delimiter=",",
           comments="")
print(f"Saved: {outcsv}")

# ── Plot ───────────────────────────────────────────────────────
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 7), sharex=True)

# Panel 1: fraction of peptides in beta-sheets
ax1.plot(time_ns, beta_fraction * 100, linewidth=0.5, color="tab:red")
ax1.set_ylabel("Peptides in β-sheet (%)", fontsize=12)
ax1.set_title("Beta-sheet Formation", fontsize=12)
ax1.set_ylim(-5, 105)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)

# Panel 2: largest beta-sheet cluster size
ax2.plot(time_ns, largest_sheet, linewidth=0.5, color="tab:orange")
ax2.set_xlabel("Time (ns)", fontsize=12)
ax2.set_ylabel("Largest β-sheet cluster (peptides)", fontsize=12)
ax2.set_ylim(bottom=0)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)

plt.tight_layout()

outpng = os.path.join(img_dir, "beta_sheet_fraction.png")
plt.savefig(outpng, dpi=300, bbox_inches="tight")
print(f"Saved: {outpng}")

# ── Summary ────────────────────────────────────────────────────
print(f"\nMean beta-sheet fraction = {np.mean(beta_fraction)*100:.1f}%")
print(f"Max  largest sheet size  = {np.max(largest_sheet)}")
