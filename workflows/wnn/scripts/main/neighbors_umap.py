import scanpy

scanpy.settings.n_jobs = snakemake.threads


adata = scanpy.read_h5ad(snakemake.input.anndata) # type: ignore


scanpy.pp.neighbors(
    adata, n_neighbors=20,
    metric='cosine', use_rep=snakemake.params.latent_key
) # type: ignore

scanpy.tl.umap(adata, min_dist=0.3)


adata.write_h5ad(filename=snakemake.output.anndata, compression='gzip') # type: ignore