#!/bin/bash
#
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=fast
#SBATCH --mem=10G

srun ./SSSP_MPI --source 1 --inputFile input_graphs/testGraphConverted --y_or_n yes

