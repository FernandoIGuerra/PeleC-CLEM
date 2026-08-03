#!/bin/bash
#FILENAME: SR

#SBATCH -A MCH240024
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH -J SR     # Job name
#SBATCH --time=12:00:00
#SBATCH -p wholenode

module --force purge
module load gcc
module load openmpi
module load cmake
module load eigen/3.3.9

# make TPLrealclean && make realclean && make TPL > tpl.log
#make -j 8 > make.log
mpirun -np $SLURM_NTASKS PeleC3d.gnu.MPI.ex rectflame.inp > PeleC.log

