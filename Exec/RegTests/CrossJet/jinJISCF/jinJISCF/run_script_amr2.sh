#!/bin/bash
#FILENAME: SR

#SBATCH -A MCH240024
#SBATCH --nodes=4
#SBATCH --ntasks=512
#SBATCH -J SR     # Job name
#SBATCH --time=24:00:00
#SBATCH -p wholenode

module --force purge
module load gcc
module load openmpi
module load cmake


# make TPLrealclean && make realclean && make TPL > tpl.log
#make -j 8 > make.log
mpirun -np $SLURM_NTASKS PeleC3d.gnu.MPI.ex rectflame.inp > PeleC.log

