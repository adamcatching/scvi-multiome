import episcanpy as epi
import scanpy as sc

# Load in data
adata = sc.read_h5ad(snakemake.input.atac_anndata) 

# Filter for quality control on general values

# Threshold below a given mitochondria percent
epi.pp.filter_cells(adata, min_features=1000)
epi.pp.filter_features(adata, min_cells=5)

adata.write(filename=snakemake.output.atac_anndata, compression='gzip')