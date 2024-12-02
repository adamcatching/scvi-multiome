import pandas

# Define the working directory, explictly
work_dir = '/data/CARD_singlecell/SN_atlas/'

# Number of threads to use when running the rules
num_workers = 8

# Define where the metadata data exists for each sample to be processed
input_table = 'input/SNsamples.csv'

# Read in the list of 
datasets = pandas.read_csv(input_table)['Sample'].tolist()

envs = {
    'singlecell': 'envs/single_cell_cpu.yml', 
    'single_cell_gpu': 'envs/single_cell_gpu.yml', 
    'cellbender': 'envs/cellbender.yml', 
    'muon': 'envs/muon.yml'
    }

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
        expand(
            work_dir+'data/samples/{dataset}/02_{dataset}_anndata_filtered_atac.h5ad', 
            dataset=datasets
            )

rule move_data:
    script:
        work_dir+'scripts/move_cellranger.sh'

"""
rule cellbender:
    script:
        work_dir+'scripts/cellbender_array.sh'

rule pileup:
    script:
        'scripts/pileup.sh'
"""   
# ADDING GVCF, QTL, work here

rule preprocess:
    input:
        input_table=input_table,
        rna_anndata = work_dir+'data/samples/{dataset}/cellbender_gex_counts_filtered.h5'
    output:
        rna_anndata = work_dir+'data/samples/{dataset}/01_{dataset}_anndata_object_rna.h5ad'
    conda:
        envs['singlecell']
    params:
        sample='{dataset}', 
        data_root='data/CARD_singlecell/SN_atlas/cellbender'
    resources:
        runtime=120, mem_mb=64000, disk_mb=10000, slurm_partition='quick' 
    script:
        work_dir+'scripts/preprocess.py'

rule merge_unfiltered:
    input:
        rna_anndata=expand(
            work_dir+'data/samples/{dataset}/01_{dataset}_anndata_object_rna.h5ad', 
            dataset=datasets
            )
    output:
        merged_rna_anndata = work_dir+'data/atlas/01_merged_anndata_rna.h5ad'
    conda:
        envs['singlecell']
    resources:
        runtime=240, mem_mb=1500000, disk_mb=10000, slurm_partition='largemem' 
    script:
        work_dir+'scripts/merge_anndata.py'

rule plot_qc_rna:
    input:
        merged_rna_anndata = work_dir+'data/atlas/01_merged_anndata_rna.h5ad'
    conda:
        envs['muon']
    resources:
        runtime=960, mem_mb=500000, disk_mb=10000, slurm_partition='largemem' 
    script:
        work_dir+'scripts/plot_qc_metrics.py'

rule filter_rna:
    input:        
        rna_anndata = work_dir+'data/samples/{dataset}/01_{dataset}_anndata_object_rna.h5ad'
    output:
        rna_anndata = work_dir+'data/samples/{dataset}/02_{dataset}_anndata_filtered_rna.h5ad'
    conda:
        envs['muon']
    params:
        mito_percent_thresh = mito_percent_thresh,
        doublet_thresh = doublet_thresh,
        min_genes_per_cell = min_genes_per_cell,
        ribo_percent_thresh = ribo_percent_thresh
    resources:
        runtime=120, mem_mb=100000, disk_mb=10000, slurm_partition='quick' 
    script: 
        work_dir+'scripts/rna_filter.py'

rule merge_filtered_rna:
    input:
        rna_anndata=expand(
            work_dir+'data/samples/{dataset}/02_{dataset}_anndata_filtered_rna.h5ad', 
            dataset=datasets
            )
    output:
        merged_rna_anndata = work_dir+'data/atlas/02_filtered_anndata_rna.h5ad'
    conda:
        envs['singlecell']
    resources:
        runtime=120, mem_mb=1000000, disk_mb=10000, slurm_partition='largemem' 
    script:
        work_dir+'scripts/merge_anndata.py'

rule atac_preprocess:
    input:
        input_table=input_table,
        atac_anndata=work_dir+'data/samples/{dataset}/raw_feature_bc_matrix.h5'
    output:
        atac_anndata=work_dir+'data/samples/{dataset}/01_{dataset}_anndata_object_atac.h5ad'
    conda:
        envs['singlecell']
    params:
        sample='{dataset}'
    resources:
        runtime=120, mem_mb=50000, disk_mb=10000, slurm_partition='quick' 
    script:
        work_dir+'scripts/atac_preprocess.py'

rule merge_unfiltered_atac:
    input:
        atac_anndata=expand(
            work_dir+'data/samples/{dataset}/01_{dataset}_anndata_object_atac.h5ad', 
            dataset=datasets
            )
    output:
        merged_atac_anndata=work_dir+'data/atlas/01_merged_anndata_atac.h5ad'
    conda:
        envs['singlecell']
    resources:
        runtime=480, mem_mb=1500000, disk_mb=10000, slurm_partition='largemem' 
    script:
        work_dir+'scripts/merge_atac.py'

rule plot_qc_atac:
    input:
        merged_atac_anndata=work_dir+'data/atlas/01_merged_anndata_atac.h5ad'
    conda:
        envs['singlecell']
    resources:
        runtime=240, mem_mb=1500000, disk_mb=10000, slurm_partition='largemem'
    script:
        work_dir+'scripts/atac_plot_qc.py'

rule filter_atac:
    input:
        atac_anndata = work_dir+'data/samples/{dataset}/01_{dataset}_anndata_object_atac.h5ad'
    output:
        atac_anndata = work_dir+'data/samples/{dataset}/02_{dataset}_anndata_filtered_atac.h5ad'
    conda:
        envs['singlecell']
    resources:
        runtime=30, mem_mb=50000, slurm_partition='quick'
    script:
        work_dir+'scripts/atac_filter.py'

rule merge_filtered_atac:
    input:
        atac_anndata = expand(
            work_dir+'data/samples/{dataset}/02_{dataset}_anndata_filtered_atac.h5ad', 
            dataset=datasets
            )
    output:
        merged_atac_anndata = work_dir+'data/atlas/02_filtered_anndata_atac.h5ad'
    conda:
        envs['singlecell']
    resources:
        runtime=240, mem_mb=1500000, disk_mb=10000, slurm_partition='largemem'
    script:
        work_dir+'scripts/merge_atac.py' 

rule rna_atac_filter:
    input:
        merged_rna_anndata = work_dir+'data/atlas/02_filtered_anndata_rna.h5ad',
        merged_atac_anndata = work_dir+'data/atlas/02_filtered_anndata_atac.h5ad'
    output:
        merged_rna_anndata = work_dir+'data/atlas/03_filtered_anndata_rna.h5ad',
        merged_atac_anndata = work_dir+'data/atlas/03_filtered_anndata_atac.h5ad'
    conda:
        envs['singlecell']
    resources:
        runtime=240, mem_mb=2000000, slurm_partition='largemem'
    script:
        'scripts/filter_rna_atac.py'





rule rna_model:
    input:
        merged_rna_anndata = work_dir+'data/atlas/03_filtered_anndata_rna.h5ad',
        marker_table = work_dir+'input/PanglaoDB_markers_27_Mar_2020.tsv'
    output:
        merged_rna_anndata = work_dir+'data/atlas/04_modeled_anndata_rna.h5ad',
        model_history = work_dir+'model_elbo/rna_model_history.csv'
    params:
        model = work_dir+'data/models/rna'
    conda:
        envs['single_cell_gpu']
    threads:
        64
    resources:  
        runtime=2880, disk_mb=500000, mem_mb=300000, gpu=4, gpu_model='v100x'
    script:
        'scripts/rna_model.py'

rule atac_model:
    input:
        merged_atac_anndata = work_dir+'data/atlas/02_filtered_anndata_atac.h5ad'
    output:
        merged_rna_anndata = work_dir+'data/atlas/04_modeled_anndata_atac.h5ad',
        model_history = work_dir+'model_elbo/atac_model_history.csv'
    params:
        model = work_dir+'data/models/atac'
    conda:
        envs['singlecell']
    threads:
        64
    resources:
        runtime=2880, mem_mb=1500000, slurm_partition='largemem' 
    script:
        'scripts/atac_model.py'

rule annotate:
    input:
        merged_rna_anndata = work_dir+'data/atlas/04_modeled_anndata_rna.h5ad',
    output:
        merged_rna_anndata = work_dir+'data/atlas/05_annotated_anndata_rna.h5ad'
    conda:
        envs['singlecell']
    resources:
        runtime=240, mem_mb=500000, slurm_partition='largemem'
    script:
        'scripts/annotate.py'

rule SCANVI_annot:
    input:
        merged_rna_anndata = work_dir+'data/atlas/05_annotated_anndata_rna.h5ad',
        model = work_dir+'data/models/rna'
    output:
        merged_rna_anndata = work_dir+'data/atlas/06_SCVI_anndata_rna.h5ad',
        model_history = work_dir+'model_elbo/SCANVI_model_history.csv'
    conda:
        envs['single_cell_gpu']
    threads:
        64
    resources:
        runtime=2880, disk_mb=500000, mem_mb=300000, gpu=1, gpu_model='v100x'
    script:
        'scripts/SCANVI_annot.py'




"""
rule apply_atac_annotation:
rule DGE:
rule DA:
# """