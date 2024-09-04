import ray
import scvi
import torch
import scanpy
import pickle
from scvi import autotune

scvi.settings.seed = 0
torch.set_float32_matmul_precision('high')


adata = scanpy.read_h5ad(snakemake.input.anndata)

model = scvi.external.POISSONVI

model.setup_anndata(adata, layer='counts', batch_key='sample')

tuner = autotune.ModelTuner(model)


search_space = {
    'reduce_lr_on_plateau': True,
    'lr': ray.tune.loguniform(1e-4, 1e-2),
    'n_hidden': ray.tune.choice([128, 256]),
    'lr_patience': ray.tune.choice([10, 15]),
    'lr_factor': ray.tune.loguniform(0.01, 0.9),
    'n_latent': ray.tune.choice([i for i in range(10, 101, 10)])
}


ray.init(log_to_driver=False, num_gpus=snakemake.resources.gpu)

results = tuner.fit(
    adata,
    metric='validation_loss',
    search_space=search_space,
    num_samples=50, max_epochs=50,
    resources={'cpu': snakemake.threads / snakemake.resources.gpu, 'gpu': 1},
)


arg_dict = {
    'model_args': results.model_kwargs, 'train_args': results.train_kwargs
}


with open(snakemake.output.pickle, 'wb') as f:
    pickle.dump(arg_dict, f)
