#!/bin/bash

#SBATCH --job-name=cell_human
#SBATCH --time=1:00:00
#SBATCH --account=nn2849k
#SBATCH --mem-per-cpu=1G
#SBATCH --output=Output_File.out
#SBATCH --cpus-per-task=16

export OMP_NUM_THREADS=16

time ./a.out


