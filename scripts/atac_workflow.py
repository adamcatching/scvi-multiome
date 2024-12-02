import scvi
import scanpy as sc
import pickle
import torch
import scipy

scvi.settings.seed = 0
torch.set_float32_matmul_precision('high')

adata = sc.read_h5ad(filename='/data/CARD_singlecell/SN_atlas/data/atlas/03_filtered_anndata_atac.h5ad')

sc.pp.filter_genes(adata, min_cells=int(10000))

adata.layers['peaks'] = scipy.sparse.csr_matrix(adata.layers['peaks'].copy())

adata.X = adata.layers['peaks'].copy()

print("# regions before filtering:", adata.shape[-1])

# compute the threshold: 1% of the cells
min_cells = int(adata.shape[0] * 0.01)

# Filter the number of
sc.pp.filter_genes(adata, min_cells=min_cells)

print("# regions before filtering:", adata.shape[-1])

"""

sc.pp.normalize_total(adata)

# Logarithmize the data
sc.pp.log1p(adata)

# Save the normalized-log data
adata.layers['log_normal_peaks']=adata.X.copy() 

#adata.layers['fragments'] = adata.X.copy()

sc.pp.highly_variable_genes(
    adata, 
    flavor='seurat', 
    layer='log_normal_peaks',
    subset=True, 
    batch_key='batch',
    n_top_genes=200000
)"""

# Setup how the modeling will be done
scvi.external.POISSONVI.setup_anndata(adata, layer='peaks', batch_key='sample')

# Initialize the model
model = scvi.external.POISSONVI(
    adata, 
    n_latent=100, 
    n_layers=2, 
    )

# Train the model 
model.train(
    max_epochs=1000,
    accelerator='cpu',  
    early_stopping=True,
    early_stopping_patience=20
    )

# Save the ELBO of the model
elbo = model.history['elbo_train']
elbo['elbo_validation'] = model.history['elbo_validation']

# Export model parameters
model.save('data/CARD_singlecell/SN_atlas/data/models/atac', overwrite=True) # type: ignore
# Add the model parameters to the object
adata.obsm['X_poissonvi'] = model.get_latent_representation() # type: ignore
# Compute nearest neighbors from the model
sc.pp.neighbors(adata, use_rep='X_poissonvi')
# Cluster from nearest neighbors
sc.tl.leiden(adata, resolution=.5, key_added='leiden_05')
# Compute the UMAP projection
sc.tl.umap(adata, min_dist=0.3)

# Save the model and the object
adata.write_h5ad(filename='/data/CARD_singlecell/SN_atlas/data/atlas/04_modeled_anndata_atac.h5ad', compression='gzip')
elbo.to_csv('data/model_elbo/atac/model_history.csv', index=False) # type: ignore
