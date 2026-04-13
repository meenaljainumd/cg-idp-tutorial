#!/bin/csh
#SBATCH -J SS14                      # job name (change per system)
#SBATCH -o run.o%j                      # stdout/stderr log file
#SBATCH -n 128                          # number of cores
#SBATCH -t 00-12:15:00                  # wall time (DD-HH:MM:SS)
#SBATCH --mail-user=your_email@umd.edu  # email for job notifications
#SBATCH --mail-type=all                 # notify on start, end, and failure
#SBATCH -A matysiak-prj-eng             # allocation account

module load gromacs/2019.4/gcc

# Change this to the absolute path of your simulation directory
cd /home/your_username/scratch/cg-idp-tutorial/SS14/with_hp

# Run GROMACS with the ProMPT tabulated potentials.
# The -cpi flag tells mdrun to resume from the checkpoint file if one exists,
# which makes the script safe to re-submit after a wall-time interruption.
mpirun -np 128 gmx_mpi mdrun \
    -deffnm run \
    -v \
    -cpi run \
    -tableb ../../forcefield/Tables/table*
