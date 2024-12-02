import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import scanpy as sc
import scipy
import muon as mu
import scvi
import anndata as ad

# Read the samples table once
samples = pd.read_csv(snakemake.input.input_table) # type: ignore

# Extract the metadata for the specific sample in one step
metadata = samples[samples['Sample_ID'] == snakemake.params.sample].iloc[0]

# Import all peaks to a muon object
mdata = mu.read_10x_h5(snakemake.input.atac_anndata)

# Split up the two modalities
atac = mdata['atac']
rna = mdata['rna']
rna.var_names_make_unique()

# Save the raw layer
atac.layers['peaks'] = atac.X.copy()
atac.raw = atac
atac.obs['cell_barcode'] = atac.obs_names

# Add metadata to the AnnData object directly from the metadata dataframe
metadata_dict = {
    'batch': metadata['Use_batch'],
    'sex': metadata['Sex'],
    'age': metadata['Age'],
    'pmi': metadata['PMI'],
    'ethnicity': metadata['Ethnicity'],
    'race': metadata['Race'],
    'brain_bank': metadata['Brain_bank'],
    'homogenization': metadata['Homogenizing_batch'],
    'library': metadata['Library_batch'],
    'seq': metadata['Sequencing_batch'],
    'sample': metadata['Sample_ID'],
    'short diagnosis': metadata['Short Diagnosis']
}

# Save the metadata into the ATAC object
for key, value in metadata_dict.items():
    atac.obs[key] = value

# Calculate QC metrics
sc.pp.calculate_qc_metrics(atac, percent_top=None, log1p=False, inplace=True)

atac.obs['NS']=1

# Calculate nucleosome signal
mu.atac.tl.nucleosome_signal(atac, n=1e6)

# Calculate translation start site enrichment
mu.atac.tl.get_gene_annotation_from_rna(mdata['rna']).head(3)
mu.atac.tl.tss_enrichment(mdata, n_tss=1000) 

# Save the AnnData object
atac.write(filename=snakemake.output.atac_anndata, compression='gzip')