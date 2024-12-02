import scvi
import scanpy as sc
import torch
import pandas as pd
import scipy
import numpy as np

scvi.settings.seed = 0
torch.set_float32_matmul_precision('high')

adata = sc.read_h5ad(snakemake.input.merged_rna_anndata)

adata.layers['counts'] = scipy.sparse.csr_matrix(adata.layers['counts'].copy())

sc.pp.highly_variable_genes(
    adata,
    flavor="seurat",
    batch_key="sample",
    subset=True,
    n_top_genes=5000
)

scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="sample")

model = scvi.model.SCVI(
    adata, 
    dispersion="gene-batch", 
    n_layers=2, 
    n_latent=30, 
    gene_likelihood="nb"
)

model.train(
    max_epochs=1000,
    accelerator='gpu',  
    early_stopping=True,
    early_stopping_patience=20
)

elbo = model.history['elbo_train']
elbo['elbo_validation'] = model.history['elbo_validation']

adata.obsm['X_scvi'] = model.get_latent_representation()
sc.pp.neighbors(adata, use_rep='X_scvi')
sc.tl.umap(adata, min_dist=0.3)
sc.tl.leiden(adata,  resolution=.5, key_added='leiden_05')

adata.write_h5ad(snakemake.output.merged_rna_anndata, compression='gzip')
elbo.to_csv(snakemake.output.model_history, index=False)
model.save(snakemake.params.model, overwrite=True)
