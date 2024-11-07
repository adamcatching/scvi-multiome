# SCVI-Multiome

This is where the current pipeline for variational inference of multiomic single-cell analysis goes.

This pipeline is currently optimized to run on a Slurm-based HPC with the command, from within the scvi-multiome directory, sbatch snakemake.sh. The first three rules, move_data, cellbender, and pileup, are memory and GPU intensive and should only be run once.

Beginning to build one consolidated snakemake based workflow to take from CellRanger output to Differential Gene Expression and Differential Accessibility Regions. 
