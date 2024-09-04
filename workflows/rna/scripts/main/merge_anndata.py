import scanpy
import anndata


adata = anndata.concat(
    merge='same', index_unique='_', join='outer',
    adatas={[item for item in dataset.split('_') if 'ARC' in item][0]: 
            scanpy.read_h5ad(dataset) for dataset in snakemake.input.anndata} # type: ignore
)


adata.write_h5ad(filename=snakemake.output.anndata, compression='gzip') # type: ignore