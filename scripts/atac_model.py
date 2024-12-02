import scvi
import scanpy as sc
import pandas as pd
import numpy as np
import scipy

scvi.settings.seed = 0

adata = sc.read_h5ad(filename=snakemake.input.merged_atac_anndata)

adata.layers['peaks'] = scipy.sparse.csr_matrix(adata.layers['peaks'].copy())

adata.X = adata.layers['peaks'].copy()

sc.pp.normalize_total(adata, layer='peaks')

# Logarithmize the data
sc.pp.log1p(adata, layer='peaks')

# Save the normalized-log data
adata.layers['log_normal_peaks']=adata.X.copy() 

sc.pp.highly_variable_genes(
    adata, 
    layer='log_normal_peaks', 
    flavor='seurat_v3', 
    subset=True, batch_key='sample',
    n_top_genes=250000
)

# Setup how the modeling will be done
scvi.external.POISSONVI.setup_anndata(
    adata, 
    layer="peaks", 
    batch_key="sample"
    )

# Initialize the model
model = scvi.external.POISSONVI(
    adata, 
    n_latent=30, 
    n_layers=2, 
    latent_distribution='normal'
    )

model.train(
    max_epochs=1000,
    early_stopping=True,
    accelerator='cpu',
    early_stopping_patience=40,
    )

elbo = model.history['elbo_train']
elbo['elbo_validation'] = model.history['elbo_validation']


model.save(snakemake.params.model, overwrite=True)
atac.obsm['X_poissonvi'] = model.get_latent_representation() 

sc.pp.neighbors(
    adata, n_neighbors=20,
    metric='cosine', use_rep='X_poissonvi'
    )

sc.tl.umap(adata, min_dist=0.3)

adata.write_h5ad(snakemake.output.merged_rna_anndata, compression='gzip')
elbo.to_csv(snakemake.output.model_history, index=False)