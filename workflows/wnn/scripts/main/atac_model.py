import scvi
import torch
import scipy
import scanpy


scvi.settings.seed = 0
torch.set_float32_matmul_precision('high')


with open(snakemake.input.pickle, 'rb') as file:
    arg_dict = pickle.load(file)

adata = scanpy.read_h5ad(snakemake.input.anndata) # type: ignore



scvi.external.POISSONVI.setup_anndata(adata, layer='counts', batch_key='sample') 

model = scvi.external.POISSONVI(adata, n_latent=arg_dict['model_args']['n_latent'], n_hidden=arg_dict['model_args']['n_hidden']) # type: ignore

model.train(
    max_epochs=1000,
    accelerator='gpu', 
    early_stopping=True,
    early_stopping_patience=20,
    devices=snakemake.resources.gpu, 
    plan_kwargs=arg_dict['train_args']['plan_kwargs']
)

elbo = model.history['elbo_train']
elbo['elbo_validation'] = model.history['elbo_validation']


model.save(snakemake.params.model_dir, overwrite=True) # type: ignore
adata.obsm[snakemake.params.latent_key] = model.get_latent_representation() # type: ignore


elbo.to_csv(snakemake.output.model_history, index=False) # type: ignore
adata.write_h5ad(filename=snakemake.output.anndata, compression='gzip') # type: ignore

