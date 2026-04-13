#!/bin/csh
#SBATCH -J SS14_eq                      # job name (change per system)
#SBATCH -o eq.o%j                       # stdout/stderr log file
#SBATCH -n 128                          # number of cores
#SBATCH -t 00-02:00:00                  # wall time (DD-HH:MM:SS) - 2h for eq
#SBATCH --mail-user=your_email@umd.edu  # email for job notifications
#SBATCH --mail-type=all                 # notify on start, end, and failure
#SBATCH -A matysiak-prj-eng             # allocation account

module purge
module load gromacs/2019.4/gcc

# Change this to the absolute path of your simulation directory
cd /home/your_username/scratch/cg-idp-tutorial/SS14/with_hp

# Run NPT equilibration with ProMPT tabulated potentials.
# The -cpi flag resumes from the checkpoint file if one exists,
# making the script safe to re-submit after a wall-time interruption.
mpirun -np 128 gmx_mpi mdrun -deffnm eq -v -cpi eq -tableb ../../forcefield/Tables/table*
