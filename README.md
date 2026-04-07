This tutorial walks students through setting up and running coarse-grained molecular dynamics (CG-MD) simulations of intrinsically disordered peptides using the ProMPT force field (https://pubs.acs.org/doi/10.1021/acs.jctc.2c00269) in GROMACS. It was developed as course material for **BIOE464**, taught by **Dr. Silvina Matysiak** at the University of Maryland.

Students will learn to set up and run CG-MD simulations of three
peptide systems, somatostatin (SS14), cyclic somatostatin (CSS14), and pituitary adenylate cyclase-activating polypeptide (PACAP27), with and without heparin. For details on how heparin was parameterized, please refer to https://doi.org/10.1039/D4CP02331E. 

## Systems

| System  | Description              | Length |
|---------|--------------------------|--------|
| SS14    | Somatostatin-14          | 14 aa  |
| CSS14   | Cyclic somatostatin-14   | 14 aa  |
| PACAP27 | PACAP-27                 | 27 aa  |

Each peptide folder contains two setups: `with_hp/` (with heparin dp18)
and `no_hp/` (peptide in water).

## Repository layout
'''
cg-idp-tutorial
│   README.md
│
└───forcefield/
│   │   ff_1hp_8PACAP27.itp   # system-level ff file for 8 PACAP27 + 1 heparin
│   │   ff_1hp_8SS14.itp      # system-level ff file for 8 SS14 + 1 heparin
│   │   heparin.itp           # topology for dp18 heparin
│   │   martini_v2.0_ions.itp # topology for MARTINI monovalent ions
│   │   water.em.itp          # polarizable water topology (energy min)
│   │   water.md.itp          # polarizable water topology (MD)
│   │
│   └───PEP_SS14/             # topology files (pep1.itp ... pep8.itp)
│   │                         # for the 8 SS14 peptides
│   │
│   └───PEP_CSS14/            # topology files (pep1.itp ... pep8.itp)
│   │                         # for the 8 cyclic SS14 peptides
│   │
│   └───PEP_PACAP27/          # topology files (pep1.itp ... pep8.itp)
│   │                         # for the 8 PACAP27 peptides
│   │
│   └───posre_itps/           # position restraint files for each peptide
│   │                         # and for heparin
│   │
│   └───Tables/               # tabulated potentials for ProMPT:
│                             # amino-acid-specific angles (table_a1...a40),
│                             # backbone dihedrals
│
└───gro_files/
│   │   ss14.gro              # CG starting configuration, SS14
│   │   css14.gro             # CG starting configuration, cyclic SS14
│   │   pacap27.gro           # CG starting configuration, PACAP27
│   │   heparin_dp18.gro      # CG starting configuration, dp18 heparin
│   │   water_8.gro           # polarizable water box for solvation
│
└───mdp_files/
│   │   em.mdp                # energy minimization parameters
│   │   eq.mdp                # NPT equilibration parameters
│
└───scripts/
│   │   aa2cg.py              # convert an all-atom peptide PDB to a
│   │                         # ProMPT-ready CG .gro using martinize2
│   │   build_ff_itp.ipynb    # notebook to build the system-level ff.itp
│   │   genitp_batch_of_peptides.py   # generate pep1.itp ... pep8.itp
│   │                                 # for a batch of identical peptides
│   │   tpr2pdb.py            # extract a PDB from a .tpr + .xtc pair
│   │                         # for visualization in VMD
│   │   submit.sh             # sample SLURM submission script for
│   │                         # running production MD on UMD Zaratan
│   │   martini_v2.P_supp-material.itp # supplementary file for build_ff_itp.ipynb
│
└───SS14/                     # per-system setup for SS14
│   │   8pep.gro              # 8 peptides inserted in a 12x12x12 box
│   │   8pep_1hp.gro          # 8 peptides + 1 dp18 heparin
│   │
│   └───no_hp/                # 8 peptides, no heparin
│   │   │   eq.gro            # equilibrated CG structure
│   │   │   topol.top         # system topology
│   │
│   └───with_hp/              # 8 peptides + 1 dp18 heparin
│       │   eq.gro
│       │   topol.top
│
└───CSS14/                    # per-system setup for cyclic SS14
│   │   8pep.gro
│   │   8pep_1hp.gro
│   │
│   └───no_hp/
│   │   │   eq.gro
│   │   │   topol.top
│   │
│   └───with_hp/
│       │   eq.gro
│       │   topol.top
│
└───PACAP27/                  # per-system setup for PACAP27
    │   8pep.gro
    │   8pep_1hp.gro
    │
    └───no_hp/
    │   │   eq.gro
    │   │   topol.top
    │
    └───with_hp/
        │   eq.gro
        │   topol.top
'''

## Requirements

- GROMACS 2019.4 
- Python 3 with numpy and MDAnalysis
- martinize2 (`pip install vermouth`) — for `scripts/aa2cg.py`

# Tutorial: Setting up a CG-MD simulation of 8 peptides (± heparin)

This tutorial walks you through every step required to set up a coarse-grained
molecular dynamics simulation of 8 intrinsically disordered peptides, with or
without one molecule of dp18 heparin, using the ProMPT force field. We'll use
**PACAP27** as the running example — the same steps apply to SS14 and CSS14.

The tutorial assumes you are working from the root of the `cg-idp-tutorial`
repository and have GROMACS 2019.4+ available in your `$PATH`.

## Prerequisites

Install Python dependencies:

```bash
pip install numpy MDAnalysis vermouth
```

`vermouth` provides the `martinize2` command used in Step 1.

Make sure GROMACS is loaded:

```bash
gmx --version
```

## Step 1: Generate a CG representation from an all-atom PDB

If you already have a CG starting structure (for example, the ones in
`gro_files/`), you can skip to Step 2. Otherwise, use `scripts/aa2cg.py`
to convert an all-atom PDB file into a ProMPT-ready CG `.gro` file.

```bash
python ./scripts/aa2cg.py -f 2MI1.pdb -o complete.gro -ff elnedyn22 -b 12
```

The script runs `martinize2` to build a Martini 2 coarse-grained structure,
then expands it into a ProMPT representation by adding backbone dummy beads
(`BBp`, `BBm`), side-chain dummies for polar and aromatic residues, and the
extra GLU S2 bead.

**Arguments:**

| Flag | Description | Default |
|------|-------------|---------|
| `-f`, `--pdb` | Input all-atom PDB file | required |
| `-o`, `--out` | Output CG `.gro` file | `complete.gro` |
| `-ff`, `--forcefield` | Force field passed to martinize2 | `elnedyn22` |
| `-b`, `--box` | Cubic box edge length in nm | `12` |

The three starting CG structures used in this tutorial (`ss14.gro`,
`css14.gro`, `pacap27.gro`) were generated this way and are provided in
`gro_files/` so you can reproduce the full workflow without running Step 1.

## Step 2: Build the force field `.itp` files

The ProMPT force field needs two types of `.itp` files for each system:

1. A **system-level force field file** (`ff_1hp_*.itp`) that defines nonbonded
   interactions, including cation-π terms.
2. A set of **per-peptide topology files** (`pep1.itp` ... `pep8.itp`) that
   define each peptide's bonded terms.

### 2a. Build the system-level `ff.itp`

Open `scripts/build_ff_itp.ipynb` in Jupyter and change two variables:

```python
n_peptides = 8                       # number of peptides in the system
pep_seq    = "AGCKNFFWKTFTSC"        # sequence of the peptide (one-letter code)
```

Run all cells. The notebook generates the system-level force field file
containing cation-π interaction parameters. Move the resulting `.itp` into
`forcefield/` with a descriptive name, for example:

```bash
mv ff.itp forcefield/ff_1hp_8SS14.itp
```

Pre-built versions for PACAP27 and SS14 are already included in `forcefield/`
as `ff_1hp_8PACAP27.itp` and `ff_1hp_8SS14.itp`.

### 2b. Build the per-peptide `.itp` files

First, create a one-line text file containing your peptide sequence, for
example `seq.txt`:

```
HSDGIFTDSYSRYRKQMAVKKYLAAVL
```

Then run `scripts/genitp_batch_of_peptides.py`. The available arguments are:

| Flag | Description | Default |
|------|-------------|---------|
| `-o`, `--output_dir` | Output directory | required |
| `-f`, `--input_seq` | Path to text file containing the peptide sequence | required |
| `-count`, `--peptide_count` | Number of peptides in the system | required |
| `-bbdih`, `--backbone_dih_type` | Type of backbone dihedral (6 = double-well potential) | `6` |
| `-bbdih_fc`, `--backbone_dih_force_constant` | Force constant for backbone dihedrals (kJ/mol) | `2` |
| `-cys`, `--make_cys_bridge` | Turn on a disulfide bridge between the two cysteines in the peptide | `0` |

**Important:** for **CSS14** (cyclic somatostatin-14), set `-cys 1` to enable
the Cys3–Cys14 disulfide bridge that makes the peptide cyclic:

```bash
python ./scripts/genitp_batch_of_peptides.py \
    -o forcefield/PEP_CSS14/ \
    -f seq.txt \
    -count 8 \
    -bbdih 6 \
    -bbdih_fc 0 \
    -cys 1
```

For SS14 and PACAP27, use `-cys 0` (and point `-o` to the corresponding
`PEP_SS14/` or `PEP_PACAP27/` directory).

This produces 8 files, `pep1.itp` through `pep8.itp`, one per peptide in the
system. Pre-built versions for all three peptides are already in
`forcefield/PEP_SS14/`, `forcefield/PEP_CSS14/`, and `forcefield/PEP_PACAP27/`.

## Step 3: Heparin topology and starting structure

The heparin topology (`heparin.itp`) and starting structure
(`heparin_dp18.gro`) were adapted from
[suhasgotla/heparin_amyloid_self-assembly](https://github.com/suhasgotla/heparin_amyloid_self-assembly)
and are provided in `forcefield/` and `gro_files/` respectively.

## Step 4: Build the simulation box

We'll build two box configurations per peptide system: one with only peptides
(`no_hp`), and one with peptides + 1 heparin (`with_hp`). The example below is
for PACAP27 — swap in `ss14.gro` or `css14.gro` for the other systems.

### 4a. Insert 8 peptides into a 12 × 12 × 12 nm box

```bash
cd PACAP27
gmx insert-molecules \
    -ci ../gro_files/pacap27.gro \
    -nmol 8 \
    -box 12 12 12 \
    -radius 0.5 \
    -o 8pep.gro
```

This `8pep.gro` is the starting configuration for the `no_hp` directory.

### 4b. Add 1 heparin molecule on top of the 8 peptides

```bash
gmx insert-molecules \
    -ci ../gro_files/heparin_dp18.gro \
    -nmol 1 \
    -f 8pep.gro \
    -radius 0.5 \
    -o 8pep_1hp.gro
```

This `8pep_1hp.gro` is the starting configuration for the `with_hp` directory.

## Step 5: Solvate, neutralize, and run energy minimization + equilibration

The next steps are identical for `no_hp` and `with_hp` — only the input
`.gro` file differs:

- **`no_hp/`** uses `8pep.gro`
- **`with_hp/`** uses `8pep_1hp.gro`

Enter the appropriate sub-directory (e.g. `PACAP27/with_hp/`) and make sure
the corresponding starting `.gro` file and `topol.top` are present, then:

### 5a. Solvate with polarizable water

The energy minimization step uses a slightly different water topology
(`water.em.itp`) than the production MD (`water.md.itp`). Before solvating,
make sure `topol.top` points to the energy-minimization version:

```bash
sed -i "s:water.md.itp:water.em.itp:g" topol.top
```

Then solvate:

```bash
gmx solvate \
    -cp 8pep_1hp.gro \
    -cs ../../gro_files/water_8.gro \
    -p topol.top \
    -o solv.gro
```

(Replace `8pep_1hp.gro` with `8pep.gro` for the `no_hp` case.)

### 5b. Neutralize with counterions

```bash
gmx grompp \
    -f ../../mdp_files/em.mdp \
    -c solv.gro -r solv.gro \
    -p topol.top \
    -o ions.tpr \
    -maxwarn 1

gmx genion \
    -s ions.tpr \
    -o solv_neutral.gro \
    -p topol.top \
    -neutral
```

### 5c. Energy minimization

```bash
gmx grompp \
    -f ../../mdp_files/em.mdp \
    -c solv_neutral.gro -r solv_neutral.gro \
    -p topol.top \
    -o em.tpr \
    -maxwarn 1

gmx mdrun \
    -v \
    -s em.tpr \
    -deffnm em \
    -tableb ../../forcefield/Tables/table*
```

The `-tableb` flag is essential — ProMPT uses tabulated potentials for
amino-acid-specific bonded interactions, and the tables in `forcefield/Tables/`
must be passed to every `mdrun` call.

### 5d. NPT equilibration

Switch the water topology back to the MD version and generate an index file for simulation with heparin (create an index group HP):

```bash
sed -i "s:water.em.itp:water.md.itp:g" topol.top
gmx make_ndx -f em.gro -o index.ndx
# r GDS | r IDO
# name GDS_IDP HP 
# q
```

Then run equilibration (with restraints on peptide and heparin backbones):

```bash
gmx grompp \
    -f ../../mdp_files/eq.mdp \
    -c em.gro \
    -n index.ndx \
    -p topol.top \
    -o eq.tpr \
    -maxwarn 1

gmx mdrun \
    -v \
    -s eq.tpr \
    -deffnm eq \
    -tableb ../../forcefield/Tables/table*
```

The resulting `eq.gro` files for each system are already provided in
`SS14/`, `CSS14/`, and `PACAP27/` (under both `no_hp/` and `with_hp/`) so you
can compare your output against the reference.

## Visualizing your trajectory in VMD

After equilibration, you can visualize the trajectory in VMD. Because `.xtc`
files do not store bond connectivity, you need to first extract a PDB file
from the `.tpr` topology using `scripts/tpr2pdb.py`:

```bash
python ../../scripts/tpr2pdb.py -s eq.tpr -x eq.xtc -o eq.pdb
```

Then open VMD with the PDB as the topology reference and the XTC as the
trajectory:

```bash
vmd -f eq.pdb eq.xtc
```

## Running on Zaratan

To run production MD on the University of Maryland's HPC cluster (zaratan) a sample SLURM submission script is provided in the `scripts/` folder. For more details, please refer to https://hpcc.umd.edu/kb/submitting-jobs/.

```bash
sbatch scripts/submit.sh
```

Before submitting, open `scripts/submit.sh` and edit the job name,
wall time, email, allocation account, and working directory to match
your setup.

## Script attribution

- `scripts/build_ff_itp.ipynb` and `scripts/genitp_batch_of_peptides.py`
  were written by **Dr. Suhas Gotla**.
- `scripts/aa2cg.py` was originally written by **Dr. Abhilash Sahoo** and updated by **Meenal Jain**.

## Further reading

For the original tutorial on setting up ProMPT simulations of amyloid peptides with heparin, see the companion repository:
[suhasgotla/heparin_amyloid_self-assembly](https://github.com/suhasgotla/heparin_amyloid_self-assembly). 