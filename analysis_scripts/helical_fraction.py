"""
Helical fraction timeseries.

Computes the fraction of residues with helical (i, i+4) dipole contacts
using BBp and BBm dummy bead distances from the ProMPT force field.
A helical contact at residue i exists if BBm(i)–BBp(i+4) or BBp(i)–BBm(i+4)
distance < cutoff.

Usage:
    python helical_fraction.py <TPR> <XTC> <CSV_FOLDER> <IMAGE_FOLDER>

Example:
    python helical_fraction.py run.tpr run.xtc results/ figures/

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
from general_metrics_functions import is_helical_dipole_contacts

# ── Input ──────────────────────────────────────────────────────
if len(sys.argv) != 5:
    print("Usage: python helical_fraction.py <TPR> <XTC> <CSV_FOLDER> <IMAGE_FOLDER>")
    sys.exit(1)

tpr     = sys.argv[1]
xtc     = sys.argv[2]
csv_dir = sys.argv[3]
img_dir = sys.argv[4]

os.makedirs(csv_dir, exist_ok=True)
os.makedirs(img_dir, exist_ok=True)

# ── Parameters ─────────────────────────────────────────────────
HELICAL_CUTOFF = 3.0  # Angstrom cutoff for BBp-BBm dipole contact

# ── Load trajectory ────────────────────────────────────────────
u = mda.Universe(tpr, xtc)

bbp = u.select_atoms("name BBp")
if bbp.n_atoms == 0:
    print("ERROR: No BBp beads found. This script requires ProMPT dummy beads.")
    sys.exit(1)

pep_num = len(bbp.fragments)
n_res = len(bbp.fragments[0].select_atoms("name BBp").resids)
n_helical_positions = n_res - 4  # (i, i+4) contacts

print(f"Number of peptides: {pep_num}")
print(f"Residues per peptide: {n_res}")
print(f"Helical contact positions per peptide: {n_helical_positions}")

# ── Compute helical fraction per frame ─────────────────────────
time_ns = []
mean_helical_fraction = []   # mean over all peptides per frame
per_peptide_helical = []     # per-peptide helical fraction per frame

for ts in tqdm(u.trajectory, desc="Computing helical contacts"):
    time_ns.append(ts.time / 1000.0)

    # Returns a list of length (n_res - 4), each element is an array
    # of shape (pep_num,) with 1 if helical contact exists at that position
    helical_result = is_helical_dipole_contacts(
        u, ts, helical_contact_cutoff=HELICAL_CUTOFF
    )

    # Stack into array: shape (n_helical_positions, pep_num)
    helical_array = np.array(helical_result)

    # Per-peptide helical fraction = fraction of (i,i+4) positions with contacts
    per_pep = np.mean(helical_array, axis=0)  # shape (pep_num,)
    per_peptide_helical.append(per_pep)

    # Mean helical fraction across all peptides
    mean_helical_fraction.append(np.mean(per_pep))

time_ns = np.array(time_ns)
mean_helical_fraction = np.array(mean_helical_fraction)
per_peptide_helical = np.array(per_peptide_helical)  # shape (n_frames, pep_num)

# ── Save CSV ───────────────────────────────────────────────────
# Save mean helical fraction timeseries
outcsv = os.path.join(csv_dir, "helical_fraction.csv")
np.savetxt(outcsv,
           np.column_stack([time_ns, mean_helical_fraction]),
           header="time_ns,mean_helical_fraction",
           delimiter=",",
           comments="")
print(f"Saved: {outcsv}")

# Save per-peptide helical fraction
outcsv2 = os.path.join(csv_dir, "helical_fraction_per_peptide.csv")
header_cols = "time_ns," + ",".join([f"pep_{i}" for i in range(pep_num)])
np.savetxt(outcsv2,
           np.column_stack([time_ns, per_peptide_helical]),
           header=header_cols,
           delimiter=",",
           comments="")
print(f"Saved: {outcsv2}")

# ── Plot ───────────────────────────────────────────────────────
fig, ax = plt.subplots()

# Panel 1: mean helical fraction over time
ax.plot(time_ns, mean_helical_fraction * 100, linewidth=0.5, color="tab:green")
ax.set_ylabel("Mean helical fraction (%)", fontsize=12)
ax.set_title("Helical Content (i, i+4 dipole contacts)", fontsize=12)
ax.set_ylim(-5, 105)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()

outpng = os.path.join(img_dir, "helical_fraction.png")
plt.savefig(outpng, dpi=300, bbox_inches="tight")
print(f"Saved: {outpng}")

# ── Summary ────────────────────────────────────────────────────
print(f"\nMean helical fraction = {np.mean(mean_helical_fraction)*100:.1f}%")
print(f"Max  helical fraction = {np.max(mean_helical_fraction)*100:.1f}%")
