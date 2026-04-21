This tutorial walks students through setting up and running coarse-grained molecular dynamics (CG-MD) simulations of intrinsically disordered peptides using the ProMPT force field (https://pubs.acs.org/doi/10.1021/acs.jctc.2c00269) in GROMACS. It is developed as course material for **BIOE464**, taught by **Dr. Silvina Matysiak** at the University of Maryland.

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
```
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
│   │   md.mdp                # production MD parameters (500 ns)
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
│   │   │   run.tpr           # ready-to-run production MD input
│   │
│   └───with_hp/              # 8 peptides + 1 dp18 heparin
│       │   eq.gro
│       │   topol.top
│       │   run.tpr
│
└───CSS14/                    # per-system setup for cyclic SS14
│   │   8pep.gro
│   │   8pep_1hp.gro
│   │
│   └───no_hp/
│   │   │   eq.gro
│   │   │   topol.top
│   │   │   run.tpr
│   │
│   └───with_hp/
│       │   eq.gro
│       │   topol.top
│       │   run.tpr
│
└───PACAP27/                  # per-system setup for PACAP27
    │   8pep.gro
    │   8pep_1hp.gro
    │
    └───no_hp/
    │   │   eq.gro
    │   │   topol.top
    │   │   run.tpr
    │
    └───with_hp/
        │   eq.gro
        │   topol.top
        │   run.tpr
```

## Requirements
- Access to UMD's Zaratan HPC cluster (GROMACS 2019.4 module)
- Python 3 with numpy, MDAnalysis, and vermouth (installed in Step 0)

# Tutorial: Setting up a CG-MD simulation of 8 peptides (± heparin)

This tutorial walks you through every step required to set up a coarse-grained
molecular dynamics simulation of 8 intrinsically disordered peptides, with or
without one molecule of dp18 heparin, using the ProMPT force field. We'll use
**SS14** as the running example, the same steps apply to CSS14 and PACAP27.

The entire workflow runs on UMD's [Zaratan](https://hpcc.umd.edu/) HPC cluster.

## Step 0: Transfer files to Zaratan

Download the repository as a ZIP from the GitHub page. From your laptop, transfer the ZIP to your Zaratan scratch directory:

```bash
scp -rp cg-idp-tutorial-main.zip <username>@login.zaratan.umd.edu:/home/<username>/scratch/
```

Enter your UMD password and approve the Duo push notification when prompted. Replace `<username>` with your Zaratan username.

Then log in to Zaratan and unzip the archive:

```bash
ssh <username>@login.zaratan.umd.edu
cd /home/<username>/scratch/
unzip cg-idp-tutorial-main.zip
cd cg-idp-tutorial-main
```

## Prerequisites (on Zaratan)

Load Python and install the Python dependencies (one-time setup):

```bash
module load python
pip install --user numpy MDAnalysis vermouth
```

The `--user` flag installs packages into your home directory since you don't have sudo access on Zaratan. `vermouth` provides the `martinize2` command used in Step 1.

Load GROMACS:

```bash
module load gromacs
gmx_mpi --version
```

On Zaratan, the GROMACS command is `gmx_mpi` (not `gmx`), since it's the MPI-enabled build.

## Step 1: Generate a CG representation from an all-atom PDB

If you already have a CG starting structure (for example, the ones in
`gro_files/`), you can skip to Step 2. Otherwise, use `scripts/aa2cg.py`
to convert an all-atom PDB file into a ProMPT-ready CG `.gro` file.

```bash
python ./scripts/aa2cg.py -f ./PDBs/2MI1.pdb -o cg_prompt.gro -ff elnedyn22 -b 12
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

## Step 2: Force field `.itp` files (pre-built)

The ProMPT force field needs two types of `.itp` files for each system:

1. A **system-level force field file** (`ff_1hp_*.itp`) that defines nonbonded
   interactions, including cation-π terms.
2. A set of **per-peptide topology files** (`pep1.itp` ... `pep8.itp`) that
   define each peptide's bonded terms.

**Pre-built versions for all three peptides are already included in the repository** (`forcefield/ff_1hp_8*.itp` and `forcefield/PEP_*/`), so you can proceed directly to Step 3. If you want to build these files from scratch see the [Optional: Building itp files yourself] section at the end of this tutorial. Note that building the system-level `ff.itp` requires a Jupyter notebook, so this is easier to do on your laptop than on Zaratan.

## Step 3: Heparin topology and starting structure

The heparin topology (`heparin.itp`) and starting structure
(`heparin_dp18.gro`) were adapted from
[suhasgotla/heparin_amyloid_self-assembly](https://github.com/suhasgotla/heparin_amyloid_self-assembly)
and are provided in `forcefield/` and `gro_files/` respectively.

## Step 4: Build the simulation box

We'll build two box configurations per peptide system: one with only peptides
(`no_hp`), and one with peptides + 1 heparin (`with_hp`). 

### 4a. Insert 8 peptides into a 12 × 12 × 12 nm box

```bash
cd SS14/
```
```bash
gmx_mpi insert-molecules -ci ../gro_files/ss14.gro -nmol 8 -box 12 12 12 -radius 0.5 -o 8pep.gro
```

This `8pep.gro` is the starting configuration for the `no_hp` directory.

### 4b. Add a heparin molecule 

```bash
gmx_mpi insert-molecules -ci ../gro_files/heparin_dp18.gro -nmol 1 -f 8pep.gro -radius 0.5 -o 8pep_1hp.gro
```

This `8pep_1hp.gro` is the starting configuration for the `with_hp` directory.

## Step 5: Solvate, neutralize, energy minimize and equilibrate

The next steps are identical for `no_hp` and `with_hp`, only the input `.gro` file differs:

- **`no_hp/`** uses `8pep.gro`
- **`with_hp/`** uses `8pep_1hp.gro`

Enter the appropriate sub-directory (e.g. `SS14/with_hp/`) and make sure
the corresponding starting `.gro` file and `topol.top` are present, then:

```bash
cd with_hp
```

### 5a. Solvate with polarizable water

```bash
gmx_mpi solvate -cp ../8pep_1hp.gro -cs ../../gro_files/water_8.gro -p topol.top -o solv.gro
```

(Replace `../8pep_1hp.gro` with `../8pep.gro` for the `no_hp` case.)

### 5b. Neutralize with counterions
Build the `.tpr` for ion placement:

```bash
gmx_mpi grompp -f ../../mdp_files/em.mdp -c solv.gro -r solv.gro -p topol.top -o ions.tpr -maxwarn 1
```

Then run `genion`. When prompted for the group to replace with ions, select **PW** (polarizable water):

```bash
echo "PW" | gmx_mpi genion -s ions.tpr -o solv_neutral.gro -p topol.top -neutral
```

> **About `-maxwarn`:** `-maxwarn 1` tells `grompp` to proceed despite one warning. If you see more warnings, increase the number (e.g. `-maxwarn 2`), but **read each warning first**. Most are harmless (e.g., net system charge before neutralization, which is exactly why we are about to run `genion`), but some indicate real problems with the topology or parameters.

### 5c. Energy minimization

The energy minimization step uses a slightly different water topology
(`water.em.itp`) than the production MD (`water.md.itp`). Make sure `topol.top` points to the 
energy-minimization version:

```bash
sed -i "s:water.md.itp:water.em.itp:g" topol.top
```
Build the `em.tpr`:

```bash
gmx_mpi grompp -f ../../mdp_files/em.mdp -c solv_neutral.gro -r solv_neutral.gro -p topol.top -o em.tpr -maxwarn 1
```

Then run energy minimization:

```bash
gmx_mpi mdrun -v -s em.tpr -deffnm em -tableb ../../forcefield/Tables/table*
```

The `-tableb` flag is essential, ProMPT uses tabulated potentials for bonded interactions, and the tables in `forcefield/Tables/` must be passed to every `mdrun` call.


### 5d. NPT equilibration
Switch the water topology back to the MD version:

```bash
sed -i "s:water.em.itp:water.md.itp:g" topol.top
```

Generate an index file with a combined heparin group `HP` (only needed for `with_hp/`; skip this step for `no_hp/`):

```bash
gmx_mpi make_ndx -f em.gro -o index.ndx
# r GDS | r IDO
# name GDS_IDO HP
# q
```

Build the `eq.tpr` (with restraints on peptide and heparin backbones), make sure to load gromacs/2019.4 version for this

```bash
module purge  #removes all the loaded modules
module load gromacs/2019.4
gmx_mpi --version 
```

```bash
gmx_mpi grompp -f ../../mdp_files/eq.mdp -c em.gro -r em.gro -n index.ndx -p topol.top -o eq.tpr -maxwarn 1
```

Equilibration MD is long enough that it must be submitted as a batch job rather than run interactively. Submit using the provided SLURM script:

```bash
sbatch ../../scripts/eq_submit.sh
```

**Before submitting**, open `scripts/eq_submit.sh` and edit the job name, wall time, email, allocation account, and the `cd` path to match your setup. Details on the fields are in the comments inside the script.

Check job status with:

```bash
squeue -u $USER
```

When the job finishes, you'll have `eq.gro` and `eq.cpt` in the working directory. Reference `eq.gro` files for each system are also provided in `SS14/`, `CSS14/`, and `PACAP27/` (under both `no_hp/` and `with_hp/`) so you can compare your output against the expected result.

## Step 6: Production MD

A 500 ns production MD parameter file (`mdp_files/md.mdp`) and pre-built `run.tpr` files for each system are included, so students can submit the production run directly. From inside the peptide sub-directory (e.g. `SS14/with_hp/`):

```bash
sbatch ../../scripts/submit.sh
```

**Before submitting**, edit `scripts/submit.sh` the same way as `eq_submit.sh` (job name, wall time, email, account, working directory).

If you want to modify production parameters (simulation length, output frequency, etc.), regenerate `run.tpr` from the equilibrated structure first:

```bash
gmx_mpi grompp -f ../../mdp_files/md.mdp -c eq.gro -p topol.top -n index.ndx -o run.tpr -maxwarn 1
```
## Visualizing your trajectory in VMD

After the simulation finishes, you can visualize the trajectory in VMD. Because `.xtc` files do not store bond connectivity, you need to first extract a PDB file from the `.tpr` topology using `scripts/tpr2pdb.py`:

```bash
module purge 
module load python 
python ../../scripts/tpr2pdb.py -s run.tpr -x run.xtc -o run.pdb
```

Then open VMD with the PDB as the topology reference and the XTC as the trajectory:

```bash
vmd -f run.pdb run.xtc
```

VMD is typically not installed on Zaratan, so copy `run.pdb` and `run.xtc` back to your laptop first with `scp`.

## Analysis

After the production run finishes, you can analyze the trajectory using
the scripts in `analysis_scripts/`. Each script takes TPR and XTC files and saves a PNG plot in the current directory. 

```bash
module purge
module load python
```

### Running the analysis scripts

Navigate to the simulation directory you want to analyze:
```bash
cd SS14/with_hp
```
Run each analysis script (one by one). Here we show how to analyze heparin's
end-to-end distance over time:

```bash
python ../../analysis_scripts/hp_e2e.py run.tpr run.xtc
```
The PNG files will appear in the current directory (e.g. `SS14/with_hp/`).

> **Note:** The heparin scripts (`hp_e2e.py`, `hp_rg.py`) only work for
> `with_hp/` systems.

### Available scripts

| Script | Property |
|---|---|
| `hp_rg.py` | Heparin radius of gyration vs time |
| `hp_e2e.py` | Heparin end-to-end distance vs time |
| `aggregation_number.py` | Weighted aggregation number (Nw) vs time |
| `beta_sheet_fraction.py` | Beta-sheet fraction and largest sheet size vs time |
| `helical_fraction.py` | Helical contact fraction vs time |

### Viewing plots on Zaratan

1. Connect to UMD VPN (GlobalProtect) or eduroam
2. Go to https://portal.zaratan.umd.edu/
3. Sign in with your Zaratan username and password
4. Click **Files → Home Directory**
5. Navigate to `cg-idp-tutorial` → your system folder (e.g. `SS14/with_hp/`)
6. Click on the PNG file to preview it

This is how the browser window looks:

![portal_zaratan](/images/portal_zaratan_umd.png)

## Optional: Building itp files yourself

If you want to generate the force field `.itp` files from scratch (for a new peptide sequence), follow the two sub-steps below. **This section is easier to run on your laptop than on Zaratan** because Step A requires a Jupyter notebook. Once generated, copy the resulting `.itp` files back to Zaratan into the appropriate `forcefield/` sub-directories.

### A. Build the system-level `ff.itp`

Open `scripts/build_ff_itp.ipynb` in Jupyter and change two variables:

```python
n_peptides = 8                       # number of peptides in the system
pep_seq    = "AGCKNFFWKTFTSC"        # sequence of the peptide (one-letter code)
```

Run all cells. The notebook generates the system-level force field file
containing cation-π interaction parameters. Move the resulting `.itp` into
`forcefield/` with a descriptive name:

```bash
mv ff.itp forcefield/ff_1hp_8SS14.itp
```

### B. Build the per-peptide `.itp` files

First, create a one-line text file containing your peptide sequence, for example `seq.txt`:

```
AGCKNFFWKTFTSC
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

```bash
python ./scripts/genitp_batch_of_peptides.py \
    -o forcefield/PEP_SS14/ \
    -f seq.txt \
    -count 8 \
    -bbdih 6 \
    -bbdih_fc 0 \
    -cys 0
```

For PACAP27, use `-cys 0` and point `-o` to `forcefield/PEP_PACAP27/`.

**Important:** for **CSS14** (cyclic somatostatin-14), set `-cys 1` to enable the Cys3–Cys14 disulfide bridge that makes the peptide cyclic, and point `-o` to `forcefield/PEP_CSS14/`.

This produces 8 files, `pep1.itp` through `pep8.itp`, one per peptide in the system.

## Script attribution

- `scripts/build_ff_itp.ipynb` and `scripts/genitp_batch_of_peptides.py`
  were written by **Dr. Suhas Gotla**.
- `scripts/aa2cg.py` was originally written by **Dr. Abhilash Sahoo** and updated by **Meenal Jain**.

## Further reading

For the original tutorial on setting up ProMPT simulations of amyloid peptides with heparin, see the companion repository:
[suhasgotla/heparin_amyloid_self-assembly](https://github.com/suhasgotla/heparin_amyloid_self-assembly).