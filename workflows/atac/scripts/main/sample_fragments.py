import os
import gzip
import pandas
import scanpy


adata = scanpy.read_h5ad(snakemake.input.anndata)



batch = adata.obs['batch'].unique()[0]
sample = adata.obs['sample'].unique()[0]
cohort = adata.obs['cohort'].unique()[0]

data_path = '/data/CARD_singlecell/Brain_atlas/'
frag_file = os.path.join(data_path, cohort.upper() + '_multiome', batch, 'Multiome', sample, 'outs/atac_fragments.tsv.gz')


with gzip.open(frag_file, 'rt') as fragments:
    df = pandas.read_csv(
        fragments, sep='\t', header=None, names=[
            'chr', 'start', 'end', 'barcode', 'count'
            ]
    )

df['sample'] = sample
df.dropna(inplace=True)
df = df[df['barcode'].isin(adata.obs_names)]

df['chr'] = df['chr'].astype('category')
df['start'] = df['start'].astype('int32')
df['end'] = df['end'].astype('int32')
df['barcode'] = df['barcode'].astype('category')
df['count'] = df['count'].astype('int32')



df.to_csv(snakemake.output.bed, sep='\t', index=False, header=False)








