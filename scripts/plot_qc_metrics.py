import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc

# Save the AnnData object
adata = sc.read_h5ad(snakemake.input.anndata) # type: ignore
adata = adata['rna']

for sample in adata.obs['sample'].to_list():

    # Make plot directory
    try:
        os.mkdir(f'plots/{sample}')
    except FileExistsError:
        print

    sample_name = sample.split('-')[0]

    # Plot percent mitochondria 
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    fig.suptitle(f' Sample {sample_name} ', fontsize=18)
    sc.pl.violin(adata[adata.obs['sample'] == sample], ['pct_counts_mt'], jitter=0.5, ax=ax[0], show=False)
    ax[0].set_ylabel('percent')
    ax[0].set_xticks('')
    ax[0].set_xlim(-.75, .75)
    ax[0].set_title('Percent mitochondria per cell')

    sns.histplot(adata[adata.obs['sample'] == sample].obs['pct_counts_mt'], ax=ax[1])
    ax[1].set_yscale("log")
    ax[1].set_xlabel('percent')
    ax[1].set_ylabel('number of cells')
    ax[1].set_title('Percent mitochondria per cell')

    plt.savefig(f'plots/{sample}/mito_pct.png', dpi=300)

    # Plot percent ribosome
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    fig.suptitle(f' Sample {sample_name} ', fontsize=18)
    sc.pl.violin(adata[adata.obs['sample'] == sample], ['pct_counts_rb'], jitter=0.5, ax=ax[0], show=False)
    ax[0].set_ylabel('percent')
    ax[0].set_xlim(-.75, .75)
    ax[0].set_title('Percent ribosome genes per cell')

    sns.histplot(adata[adata.obs['sample'] == sample].obs['pct_counts_rb'], ax=ax[1])
    ax[1].set_yscale("log")
    ax[1].set_xlabel('percent')
    ax[1].set_ylabel('number of cells')
    ax[1].set_title('Percent ribosome genes per cell')
    plt.savefig(f'plots/{sample}/ribo_pct.png', dpi=300)

    # Plot number of genes per cell
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    fig.suptitle(f' Sample {sample_name} ', fontsize=18)
    sc.pl.violin(adata[adata.obs['sample'] == sample], ['n_genes_by_counts'], jitter=0.5, ax=ax[0], show=False)
    ax[0].set_ylabel('total counts')
    ax[0].set_xlim(-.75, .75)
    ax[0].set_title('Number of genes per cell')

    sns.histplot(adata[adata.obs['sample'] == sample].obs['n_genes_by_counts'], ax=ax[1])
    ax[1].set_yscale("log")
    ax[1].set_xlabel('total counts')
    ax[1].set_ylabel('number of cells')
    ax[1].set_title('Number of genes per cell')
    plt.savefig(f'plots/{sample}/num_genes_per_cell.png', dpi=300)

    sc.pl.scatter(
        adata[adata.obs['sample'] == sample], 
        "total_counts", 
        "n_genes_by_counts", 
        color="pct_counts_mt",
        show=False
        )
    plt.savefig(f'plots/{sample}/num_gene_counts_total.png', dpi=300)

    # Plot UMAP of values
    sc.pl.umap(
        adata[adata.obs['sample'] == sample],
        color=['pct_counts_mt', 'pct_counts_rb', 'doublet_score', 'total_counts', 'n_genes_by_counts'],
        size=2,
        ncols=3,
        cmap='viridis',
        show=False
    )
    plt.savefig(f'plots/{sample}/qc_umap.png', dpi=300)