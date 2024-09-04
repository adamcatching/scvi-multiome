import muon
import scanpy


adata = muon.read_anndata(filename=snakemake.input.muon_object, mod=snakemake.params.assay) # type: ignore


scanpy.pp.filter_genes(adata, min_cells=int(adata.n_obs * snakemake.params.cutoff))


adata.write_h5ad(filename=snakemake.output.anndata, compression='gzip') # type: ignore