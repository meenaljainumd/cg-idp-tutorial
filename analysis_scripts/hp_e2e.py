"""
Heparin dp18 end-to-end distance timeseries.

Computes the distance between the first and last backbone beads
of the heparin chain over the trajectory.

Usage:
    python heparin_end_to_end.py <TPR> <XTC> <CSV_FOLDER> <IMAGE_FOLDER>

Example:
    python heparin_end_to_end.py run.tpr run.xtc results/ figures/
"""

import sys
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import MDAnalysis as mda

# ── Input ──────────────────────────────────────────────────────
if len(sys.argv) != 5:
    print("Usage: python heparin_end_to_end.py <TPR> <XTC> <CSV_FOLDER> <IMAGE_FOLDER>")
    sys.exit(1)

tpr     = sys.argv[1]
xtc     = sys.argv[2]
csv_dir = sys.argv[3]
img_dir = sys.argv[4]

os.makedirs(csv_dir, exist_ok=True)
os.makedirs(img_dir, exist_ok=True)

# ── Load trajectory ────────────────────────────────────────────
u = mda.Universe(tpr, xtc)

heparin = u.select_atoms("resname GDS IDO")
if heparin.n_atoms == 0:
    print("ERROR: No heparin (GDS/IDO) found. Is this a with_hp system?")
    sys.exit(1)

# First and last backbone beads for end-to-end distance
hp_bb = heparin.select_atoms("name BB")
if hp_bb.n_atoms < 2:
    hp_bb = heparin

first_bead = hp_bb[0]
last_bead  = hp_bb[-1]

print(f"Heparin: {heparin.n_atoms} atoms, {heparin.n_residues} residues")
print(f"End-to-end: {first_bead.name} (resid {first_bead.resid}) "
      f"<-> {last_bead.name} (resid {last_bead.resid})")

# ── Compute ────────────────────────────────────────────────────
time_ns, ree = [], []

for ts in u.trajectory:
    time_ns.append(ts.time / 1000.0)
    ree.append(np.linalg.norm(last_bead.position - first_bead.position))

time_ns = np.array(time_ns)
ree     = np.array(ree)

# ── Save CSV ───────────────────────────────────────────────────
outcsv = os.path.join(csv_dir, "heparin_ree.csv")
np.savetxt(outcsv,
           np.column_stack([time_ns, ree]),
           header="time_ns,Ree_A",
           delimiter=",",
           comments="")
print(f"Saved: {outcsv}")

# ── Plot ───────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(time_ns, ree / 10.0, linewidth=0.5, color="tab:blue")
ax.set_xlabel("Time (ns)", fontsize=12)
ax.set_ylabel("End-to-end distance (nm)", fontsize=12)
ax.set_title("End-to-end distance of heparin 18", fontsize=12)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()

outpng = os.path.join(img_dir, "heparin_e2e.png")
plt.savefig(outpng, dpi=300, bbox_inches="tight")
print(f"Saved: {outpng}")

# ── Summary ────────────────────────────────────────────────────
print(f"\nMean E2E of heparin = {np.mean(ree)/10:.2f} nm")
print(f"Std  E2E of heparin = {np.std(ree)/10:.2f} nm")

