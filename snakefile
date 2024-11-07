import pandas

num_workers = 8
input_table = 'input/SNsamples.csv'

"""
This should be updated to be more flexible 
"""
datasets = pandas.read_csv(input_table)['Sample'].tolist()

envs = {'singlecell': 'envs/sc_2.yml', 'single_cell_gpu': 'envs/single_cell_gpu.yml', 'cellbender': 'envs/cellbender.yml', 'muon': 'envs/muon.yml'}

# Define RNA thresholds
mito_percent_thresh = 15
ribo_percent_thresh = 10
doublet_thresh = 0.25
min_genes_per_cell = 500

# Define ATAC thresholds
min_peak_counts = 500
min_num_cell_by_counts = 10

localrules: all


rule all:
    input:
        expand('data/samples/{dataset}/01_{dataset}_anndata_object.h5ad', dataset=datasets)

"""rule move_data:
    script:
        'scripts/move_cellranger.sh'

rule cellbender:
    script:
        'scripts/cellbender_array.sh'

rule pileup:
    script:
        'scripts/pileup.sh'"""

rule preprocess:
    input:
        input_table=input_table,
        anndata='data/samples/{dataset}/cellbender_gex_counts_filtered.h5'
    output:
        anndata='data/samples/{dataset}/01_{dataset}_anndata_object_rna.h5ad'
    conda:
        envs['singlecell']
    params:
        sample='{dataset}', 
        data_root='data/CARD_singlecell/SN_atlas/cellbender'
    resources:
        runtime=120, mem_mb=64000, disk_mb=10000, slurm_partition='quick' 
    script:
        'scripts/preprocess.py'

rule merge_unfiltered:
    input:
        anndata=expand('data/samples/{dataset}/01_{dataset}_anndata_object_rna.h5ad', dataset=datasets)
    output:
        anndata='data/atlas/02_merged_anndata_rna.h5ad'
    conda:
        envs['singlecell']
    resources:
        runtime=120, mem_mb=1000000, disk_mb=10000, slurm_partition='largemem' 
    script:
        'scripts/merge_anndata.py'

rule plot_qc:
    input:
        anndata='data/atlas/02_merged_anndata_rna.h5ad'
    conda:
        envs['muon']
    resources:
        runtime=120, mem_mb=100000, disk_mb=10000, slurm_partition='largemem' 
    script:
        'scripts/plot_qc_metrics.py'

"""rule filter:
    input:        
        anndata='data/atlas/02_merged_anndata_rna.h5ad'
    output:
        seurat='data/atlas/03_merged_anndata_rna.h5ad'
    conda:
        envs['singlecell']
    resources:
        runtime=120, mem_mb=64000, disk_mb=10000, slurm_partition='quick' 
    script: 
        'scripts/rna_filter.py'

rule rna_workflow:
    input:
        anndata='data/atlas/03_merged_anndata_rna.h5mu',
    output:
        anndata='data/atlas/04_merged_anndata_rna.h5mu',
        model='data/models/rna/model_history.csv'
    conda:
        envs['single_cell_gpu']
    threads:
        num_workers * 16
    resources:
        runtime=2880, disk_mb=500000, mem_mb=300000, gpu=1, gpu_model='v100x'
    shell:
        'scripts/rna_workflow.py'
"""