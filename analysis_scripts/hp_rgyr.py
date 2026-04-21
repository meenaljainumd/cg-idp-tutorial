"""
Heparin dp18 radius of gyration (Rg) timeseries.
Computes Rg of the heparin chain over the trajectory.

Usage:
    python heparin_rg.py <TPR> <XTC> <CSV_FOLDER> <IMAGE_FOLDER>
Example:
    python heparin_rg.py run.tpr run.xtc results/ figures/
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
    print("Usage: python heparin_rg.py <TPR> <XTC> <CSV_FOLDER> <IMAGE_FOLDER>")
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

print(f"Heparin: {heparin.n_atoms} atoms, {heparin.n_residues} residues")

# ── Compute ────────────────────────────────────────────────────
time_ns, rg = [], []

for ts in u.trajectory:
    time_ns.append(ts.time / 1000.0)
    rg.append(heparin.radius_of_gyration())

time_ns = np.array(time_ns)
rg      = np.array(rg)

# ── Save CSV ───────────────────────────────────────────────────
outcsv = os.path.join(csv_dir, "heparin_rg.csv")
np.savetxt(outcsv,
           np.column_stack([time_ns, rg]),
           header="time_ns,Rg_A",
           delimiter=",",
           comments="")
print(f"Saved: {outcsv}")

# ── Plot ───────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(time_ns, rg / 10.0, linewidth=0.5, color="tab:green")
ax.set_xlabel("Time (ns)", fontsize=12)
ax.set_ylabel(r"R$_g$ (nm)", fontsize=12)
ax.set_title("Radius of Gyration of heparin 18", fontsize=12)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()

outpng = os.path.join(img_dir, "heparin_rg.png")
plt.savefig(outpng, dpi=300, bbox_inches="tight")
print(f"Saved: {outpng}")

# ── Summary ────────────────────────────────────────────────────
print(f"\nMean Rg = {np.mean(rg)/10:.2f} nm")
print(f"Std  Rg = {np.std(rg)/10:.2f} nm")
