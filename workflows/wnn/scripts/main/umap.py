import muon
import scanpy

scanpy.settings.n_jobs = snakemake.threads


mdata = muon.read_h5mu(snakemake.input.muon_object)


muon.tl.umap(mdata, neighbors_key='wnn', min_dist=0.3)

mdata.obsm['X_wnn_umap'] = mdata.obsm['X_umap']


mdata.write_h5mu(filename=snakemake.output.muon_object, compression='gzip') # type: ignore