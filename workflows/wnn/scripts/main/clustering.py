import muon
import scanpy

scanpy.settings.n_jobs = snakemake.threads


mdata = muon.read_h5mu(snakemake.input.muon_object)


scanpy.tl.leiden(mdata, neighbors_key='wnn', resolution=snakemake.params.resolution)

mdata.mod['rna'].obs['leiden'] = mdata.mod['atac'].obs['leiden'] = mdata.obs['leiden']


mdata.write_h5mu(filename=snakemake.output.muon_object, compression='gzip') # type: ignore