import scanpy
import anndata


adata = anndata.concat(
    merge='same', index_unique='_', join='inner',
    adatas={[item for item in dataset.split('_') if 'ARC' in item][0]: 
            scanpy.read_h5ad(dataset) for dataset in snakemake.input.anndata} # type: ignore
)


peaks = adata.var_names.str.split(pat='-', expand=True).to_frame(index=False)

peaks.set_index(adata.var.index, inplace=True)
peaks.rename({0: 'chr', 1: 'start', 2: 'end'}, axis=1, inplace=True)

adata.var = peaks.copy()

adata.var_names = adata.var_names.str.replace(pat='-', repl=':', n=1) #type: ignore


adata.write_h5ad(filename=snakemake.output.anndata, compression='gzip') # type: ignore