import scanpy as sc
import anndata as ad


adata = ad.concat(
    merge='same', index_unique='_', join='outer',
    adatas={[item for item in dataset.split('_')][-4]: 
            sc.read_h5ad(dataset) for dataset in snakemake.input.atac_anndata}
)

adata.write_h5ad(filename=snakemake.output.merged_atac_anndata, compression='gzip') # type: ignore