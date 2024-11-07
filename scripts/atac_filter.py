import muon as mu
import scanpy as sc

# Load in data
adata = sc.read_h5ad(snakemake.input.atac_anndata) # type: ignore

# Filter for quality control on general values

# Threshold below a given mitochondria percent
mu.pp.filter_var(adata, 'n_cells_by_counts', lambda x: x >= snakemake.min_num_cell_by_counts)
mu.pp.filter_obs(adata, 'n_genes_by_counts', lambda x: x <= snakemake.min_peak_counts)

adata.write(filename=snakemake.output.atac_anndata, compression='gzip') # type: ignore