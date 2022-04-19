#!/bin/bash
#
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=fast
#SBATCH --mem=10G

srun ./SSSP_Dijk --inputFile input_graphs/roadNet-CA