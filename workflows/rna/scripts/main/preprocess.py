import os
import pandas
import scanpy

# Read the samples table once
samples = pandas.read_csv(snakemake.input.input_table) # type: ignore

# Extract the metadata for the specific sample in one step
metadata = samples[samples['sample'] == snakemake.params.sample].iloc[0] # type: ignore

# Read the single-cell data
adata = scanpy.read_10x_h5(os.path.join(snakemake.params.data_root, metadata['cohort'], snakemake.params.sample, 'cellbender_gex_counts_filtered.h5')) # type: ignore

# Ensure unique variable names
adata.var_names_make_unique()

# Add mitochondrial and ribosomal markers
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['rb'] = adata.var_names.str.startswith(('RPL', 'RPS'))

# Calculate QC metrics
scanpy.pp.calculate_qc_metrics(adata, qc_vars=['rb', 'mt'], percent_top=None, log1p=False, inplace=True)

# Run scrublet to identify doublets
scanpy.external.pp.scrublet(adata, expected_doublet_rate=(adata.n_obs / 1000) * 0.008)
adata.obs.drop('predicted_doublet', axis=1, inplace=True)

# Add metadata to the AnnData object directly from the metadata dataframe
metadata_dict = {
    'batch': metadata['batch'],
    'cohort': metadata['cohort'],
    'sex': metadata['Sex'],
    'age': metadata['Age'],
    'pmi': metadata['PMI'],
    'ancestry': metadata['Ancestry'],
    'brain_bank': metadata['Brain_bank'],
    'homogenization': metadata['Homogenization'],
    'library': metadata['LibraryPrep'],
    'seq': metadata['Sequencing'],
    'sample': snakemake.params.sample
}

for key, value in metadata_dict.items():
    adata.obs[key] = value

# Save the AnnData object
adata.write_h5ad(filename=snakemake.output.anndata, compression='gzip') # type: ignore
