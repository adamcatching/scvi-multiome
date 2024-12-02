import anndata as ad
import pandas as pd
import scanpy as sc



# Open the data containing the important diagnosis
important_diagnosis_df = pd.read_csv('input/SN_diagnoses2024 - control50_PD_DLB_only.csv')
important_diagnosis_df = important_diagnosis_df[['SampleID', 'PrimaryDiagnosis']].iloc[:-3]
important_diagnosis_df['PrimaryDiagnosis'].drop_duplicates()

# Change column names for easier merge
important_diagnosis_df.columns = ['sample', 'primary diagnosis']

# Remove the dataset
important_diagnosis_df['sample'] = [str(x) for x in important_diagnosis_df['sample']]
samples_to_keep = important_diagnosis_df['sample'].to_numpy()

# Concatenate with names of sample
adata = ad.concat(
    merge='same', index_unique='_', join='outer',
    adatas={[item for item in dataset.split('_')][0]: 
            sc.read_h5ad(dataset) for dataset in snakemake.input.rna_anndata} # type: ignore
)

# Filter for samples that are useful for the dataset
adata = adata[adata.obs['sample'].isin(samples_to_keep)].copy()

# Merge the primary disease dataset
adata.obs = adata.obs.merge(important_diagnosis_df, on='sample')

# Write out the unfiltered dataset
adata.write_h5ad(filename=snakemake.output.merged_rna_anndata, compression='gzip') # type: ignore