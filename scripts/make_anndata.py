import os
import scanpy as sc
import pandas as pd
import muon as mu

dataset = snakemake.params.sample # type: ignore
input_table = pandas.read_csv(snakemake.input.samples) # type: ignore
#batch = pandas.read_csv(input_table, header=None).loc[:, 0].tolist()[1:] # type: ignore

raw_counts = os.path.join('/data/CARD_singlecell/SN_atlas/cellbender/data/sample/', dataset, 'filtered_feature_bc_matrix.h5')

mdata = mu.read_10x_h5(raw_counts)

adata.write_h5ad(filename=snakemake.output.rna, compression='gzip') # type: ignore
