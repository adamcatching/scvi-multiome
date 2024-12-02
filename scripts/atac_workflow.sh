#!/bin/sh

#SBATCH --cpus-per-task=64
#SBATCH --mem=400g
#SBATCH --time=48:00:00
#SBATCH --output=logs/atac_workflow-%j.out 
#SBATCH --partition=largemem

source /data/$USER/conda/etc/profile.d/conda.sh && source /data/$USER/conda/etc/profile.d/mamba.sh

conda activate single_cell_cpu

srun python scripts/atac_workflow.py
