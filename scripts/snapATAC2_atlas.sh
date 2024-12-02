#!/bin/bash

#SBATCH --cpus-per-task=64
#SBATCH --mem=1500g
#SBATCH --time=72:00:00
#SBATCH --output=logs/snapATAC2_atlas-%j.out 
#SBATCH --partition=largemem

source /data/$USER/conda/etc/profile.d/conda.sh && source /data/$USER/conda/etc/profile.d/mamba.sh

conda activate single_cell_cpu

srun python scripts/snapATAC2_atlas.py
