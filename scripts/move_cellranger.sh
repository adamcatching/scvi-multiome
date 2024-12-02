#!/bin/bash

#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=32G
#SBATCH --time 24:00:00
#SBATCH --output=/data/CARD_singlecell/SN_atlas/cellbender/logs/move_files

# Benjamin Adam Catching
# NIA-CARD/Data Tecnica
# 2024-10-21

"""
Move raw matrix files to directory where they 
"""

# Find the sample locations where CellRanger has been run
input_file_base=/data/CARD_singlecell/Brain_atlas/SN_Multiome/

# Define where the output files will be
output_file_base=/data/CARD_singlecell/SN_atlas/cellbender/data

# Read in the sample sheet for each batch
cat /data/CARD_singlecell/SN_atlas/data/SNsamples.csv
while IFS="," read -r Sample_ID Homogenizing_batch Library_batch Sequencing_batch Repeated Use_batch Age PMI Ethnicity Race Brain_bank Short Diagnosis
do
    # Create directories of the samples (ONLY RUN ONCE)
    mkdir $output_file_base"/sample/"$Sample_ID

    # Copy cellranger output
    cp $input_file_base"batch"$Use_batch"/Multiome/"$Sample_ID"-ARC/outs/raw_feature_bc_matrix.h5" $output_file_base"/samples/"$Sample_ID"/raw_feature_bc_matrix.h5"
    #cp $input_file_base"batch"$Use_batch"/Multiome/"$Sample_ID"-ARC/outs/filtered_feature_bc_matrix.h5" $output_file_base"/samples/"$Sample_ID"/filtered_feature_bc_matrix.h5"
    #cp $input_file_base"batch"$Use_batch"/Multiome/"$Sample_ID"-ARC/outs/atac_peak_annotation.tsv" $output_file_base"/samples/"$Sample_ID"/atac_peak_annotation.tsv"
    #cp $input_file_base"batch"$Use_batch"/Multiome/"$Sample_ID"-ARC/outs/atac_fragments.tsv.gz" $output_file_base"/samples/"$Sample_ID"/atac_fragments.tsv.gz"
    #cp $input_file_base"batch"$Use_batch"/Multiome/"$Sample_ID"-ARC/outs/atac_fragments.tsv.gz.tbi" $output_file_base"/samples/"$Sample_ID"/atac_fragments.tsv.gz.tbi"
    #cp -r $input_file_base"batch"$Use_batch"/Multiome/"$Sample_ID"-ARC/outs/filtered_feature.tsv.gz" $output_file_base"/samples/"$Sample_ID"/filtered_feature.tsv.gz"
    #rm $output_file_base"/samples/"$Sample_ID"/atac_possorted_bam.bam"
done < <(tail -n +2 /data/CARD_singlecell/SN_atlas/data/SNsamples.csv)