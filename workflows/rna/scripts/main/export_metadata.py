import scanpy
import pandas



covariates = ['sample', 'sex', 'age', 'ancestry', 'brain_bank', 'homogenization', 'library', 'seq', 'pmi', 'batch', 'cohort']

scanpy.read_h5ad(snakemake.input.anndata).obs.to_csv(snakemake.output.metadata, columns=covariates, index_label='cells')