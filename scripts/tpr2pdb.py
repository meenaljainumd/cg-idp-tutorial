"""
tpr2pdb.py — Extract a PDB file from a GROMACS .tpr + .xtc pair
               for visualization in VMD.

VMD needs atom connectivity information to render bonds correctly, which
.xtc files do not contain. By first writing a single-frame PDB from the
.tpr topology, VMD can load the trajectory (.xtc) on top of it with the
correct connectivity.

Usage:
    python tpr2pdb.py -s eq.tpr -x eq.xtc -o eq.pdb

Then visualize in VMD with:
    vmd -f eq.pdb eq.xtc
"""

import argparse
import MDAnalysis


def main():
    parser = argparse.ArgumentParser(
        description="Extract a PDB file from a GROMACS .tpr + .xtc pair."
    )
    parser.add_argument("-s", "--tpr", default="eq.tpr",
                        help="Input .tpr file (default: eq.tpr)")
    parser.add_argument("-x", "--xtc", default="eq.xtc",
                        help="Input .xtc file (default: eq.xtc)")
    parser.add_argument("-o", "--out", default="eq.pdb",
                        help="Output .pdb file (default: eq.pdb)")
    args = parser.parse_args()

    u = MDAnalysis.Universe(args.tpr, args.xtc)
    u.select_atoms("all").write(args.out)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
