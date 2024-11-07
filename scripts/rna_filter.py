import muon as mu
import scanpy as sc

# Load in data
adata = sc.read_h5ad(snakemake.input.rna_anndata) # type: ignore

# Filter for quality control on general values

# Threshold below a given mitochondria percent
mu.pp.filter_obs(adata, 'pct_counts_mt', lambda x: x <= snakemake.mito_percent_thresh)
# Threshold below a given doublet score
mu.pp.filter_obs(adata, 'doublet_score', lambda x: x < snakemake.doublet_thresh)   
# Threshold above a given genes per cell threshold
mu.pp.filter_obs(adata, 'n_genes_by_counts', lambda x: x >= snakemake.min_genes_per_cell)   
# Threshold below a given ribosome threshold
mu.pp.filter_obs(adata, 'pct_counts_rb', lambda x: x <= snakemake.ribo_percent_thresh)

adata.write(filename=snakemake.output.rna_anndata, compression='gzip') # type: ignore