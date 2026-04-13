"""
aa2cg.py — Convert an all-atom peptide PDB to a coarse-grained ProMPT-ready .gro

Usage:
    python aa2cg.py -f input.pdb [-o complete.gro] [-ff elnedyn22] [-b 12]

Pipeline:
    1. martinize2             : all-atom PDB  ->  CG PDB (Martini 2)
    2. gmx editconf / genconf : CG PDB        ->  CG .gro in a cubic box
    3. ProMPT bead expansion  : add BBp/BBm dummies, side-chain dummies,
                                and the GLU S2 bead
    4. Write final .gro       : positions, names, resnames, resnums

Requirements:
    - martinize2 installed (`pip install vermouth`)
    - GROMACS (`gmx`) in $PATH
    - The input PDB file in the working directory
"""

import os
import argparse
import numpy as np
from MDAnalysis import Universe


# ------------------------------------------------------------------ #
#  Step 1+2: All-atom PDB -> CG .gro via martinize2 + gmx            #
# ------------------------------------------------------------------ #
def run_martinize(pdb_in, cg_pdb_out, ff, box):
    """Run martinize2 on the input PDB, then editconf to a cubic box.

    Flags used:
      -ss C     : treat all residues as coil (these are IDPs)
      -noscfix  : disable side-chain corrections (not defined for elnedyn22)
    """
    os.system(
                f"martinize2 -f {pdb_in} -x {cg_pdb_out} -o topol_cg.top "
                f"-ff {ff} -ss C -noscfix"
            )

    half = box / 2.0
    os.system(
                f"gmx_mpi editconf -f {cg_pdb_out} -o CG.gro "
                f"-center {half} {half} {half} -box {box} {box} {box}"
            )

    os.system("gmx_mpi genconf -f CG.gro -o CG.gro")


# ------------------------------------------------------------------ #
#  Step 3: Expand Martini bead names (SC* -> S*)                     #
# ------------------------------------------------------------------ #
def expand_cg_beads(cg_gro_in, mart_s1_out):
    """Rename SC1..SC4 -> S1..S4 and write an intermediate .gro."""
    u = Universe(cg_gro_in)
    for old, new in [("SC1", "S1"), ("SC2", "S2"),
                     ("SC3", "S3"), ("SC4", "S4")]:
        sel = u.select_atoms(f"name {old}")
        sel.names = [new] * sel.n_atoms
    u.atoms.write(mart_s1_out)


# ------------------------------------------------------------------ #
#  Step 3 (cont): Build the expanded ProMPT bead list                #
# ------------------------------------------------------------------ #
def build_new_positions(mart_s1_file):
    """Read intermediate .gro and build the expanded ProMPT bead list."""
    with open(mart_s1_file, "r") as f:
        lines = f.readlines()[2:-1]
    A = [line.split() for line in lines]

    CurrPos = np.array(A)[:, 3:].astype(float)
    CurrNames = np.array(A)[:, 1].astype(str)
    CurrResNames = [i[-3:] for i in np.array(A)[:, 0]]
    CurrResNums = [int(i[:-3]) for i in np.array(A)[:, 0]]

    NewPos, NewNames, NewResNames, NewResNums = [], [], [], []

    def append_dummy_pair(idx, prefix):
        """Add the standard +z and +y dummy beads for a given parent bead."""
        bppos = np.array([CurrPos[idx][0], CurrPos[idx][1],
                          CurrPos[idx][2] + 0.141])
        bmpos = np.array([CurrPos[idx][0], CurrPos[idx][1] + 0.142,
                          CurrPos[idx][2]])
        NewPos.extend([bppos, bmpos])
        NewNames.extend([f"{prefix}p", f"{prefix}m"])
        NewResNums.extend([CurrResNums[idx], CurrResNums[idx]])
        NewResNames.extend([CurrResNames[idx], CurrResNames[idx]])

    for i in range(len(CurrPos)):
        NewPos.append(CurrPos[i])
        NewNames.append(CurrNames[i])
        NewResNames.append(CurrResNames[i])
        NewResNums.append(CurrResNums[i])

        # Backbone dummies for every residue
        if CurrNames[i] == "BB":
            append_dummy_pair(i, "BB")

        # Side-chain dummies for polar / aromatic residues
        if CurrResNames[i] in ["SER", "THR", "ASN", "GLN"] and CurrNames[i] == "S1":
            append_dummy_pair(i, "S1")
        if CurrResNames[i] in ["TRP", "HIS"] and CurrNames[i] == "S2":
            append_dummy_pair(i, "S2")
        if CurrResNames[i] in ["TYR", "HIS"] and CurrNames[i] == "S3":
            append_dummy_pair(i, "S3")

        # Extra GLU S2 bead, placed 0.27 nm along the BB->S1 vector
        if CurrResNames[i] == "GLU" and CurrNames[i] == "S1":
            Vec = CurrPos[i] - CurrPos[i - 1]
            UnitVec = Vec / np.linalg.norm(Vec)
            NewVec = CurrPos[i] + 0.27 * UnitVec
            NewPos.append(NewVec)
            NewNames.append("S2")
            NewResNums.append(CurrResNums[i])
            NewResNames.append(CurrResNames[i])

    return NewPos, NewNames, NewResNames, NewResNums


# ------------------------------------------------------------------ #
#  Step 4: Write the final .gro file directly                        #
# ------------------------------------------------------------------ #
def write_gro(positions, names, resnames, resnums, out_name,
              box=12.0, title="ProMPT CG peptide"):
    """Write a GROMACS .gro file directly from arrays.

    The .gro format is a plain-text fixed-width format:
        line 1     : title
        line 2     : number of atoms
        lines 3..N : %5d resnum, %-5s resname, %5s name, %5d atomnum,
                     %8.3f x, %8.3f y, %8.3f z   (positions in nm)
        last line  : box vectors (nm)
    """
    n = len(positions)
    with open(out_name, "w") as f:
        f.write(f"{title}\n")
        f.write(f"{n:5d}\n")
        for i in range(n):
            x, y, z = positions[i]
            f.write(
                     f"{resnums[i]:5d}{resnames[i]:<5s}{names[i]:>5s}"
                     f"{i + 1:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n"
                    )
        f.write(f"{box:10.5f}{box:10.5f}{box:10.5f}\n")


# ------------------------------------------------------------------ #
#  Entry point                                                       #
# ------------------------------------------------------------------ #
def main():
    parser = argparse.ArgumentParser(
        description="Convert an all-atom peptide PDB to a ProMPT CG .gro file."
    )
    parser.add_argument("-f", "--pdb", required=True,
                        help="Input all-atom PDB file (e.g. 2MI1.pdb)")
    parser.add_argument("-o", "--out", default="complete.gro",
                        help="Output CG .gro file (default: complete.gro)")
    parser.add_argument("-ff", "--forcefield", default="elnedyn22",
                        help="Force field passed to martinize2 (default: elnedyn22)")
    parser.add_argument("-b", "--box", type=float, default=12.0,
                        help="Cubic box edge length in nm (default: 12)")
    args = parser.parse_args()

    if not os.path.isfile(args.pdb):
        raise FileNotFoundError(f"Input PDB not found: {args.pdb}")

    print(f"[1/3] Running martinize2 on {args.pdb}")
    run_martinize(args.pdb, "CG.pdb", args.forcefield, args.box)

    print(f"[2/3] Expanding CG bead names (SC* -> S*)")
    expand_cg_beads("CG.gro", "mart_s1.gro")

    print(f"[3/3] Building ProMPT bead list and writing {args.out}")
    NewPos, NewNames, NewResNames, NewResNums = build_new_positions("mart_s1.gro")
    write_gro(NewPos, NewNames, NewResNames, NewResNums,
              args.out, box=args.box)

    print(f"Done. Final CG structure written to: {args.out}")


if __name__ == "__main__":
    main()
