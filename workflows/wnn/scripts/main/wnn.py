import muon
import scanpy


mdata = muon.MuData(
    {
        'rna': scanpy.read_h5ad([mod for mod in snakemake.input.anndata if 'rna' in mod][0]), # type: ignore
        'atac': scanpy.read_h5ad([mod for mod in snakemake.input.anndata if 'atac' in mod][0]) # type: ignore
    }
)


muon.pp.neighbors(mdata, key_added='wnn', metric='cosine', low_memory=True)


mdata.write_h5mu(filename=snakemake.output.muon_object, compression='gzip') # type: ignore