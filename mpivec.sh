#!/bin/bash
#
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=fast
#SBATCH --mem=10G

srun ./MPI_Vec --source 4500 --inputFile Wiki