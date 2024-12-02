import numpy as np
import scanpy as sc
import decoupler as dc
import scvi

# Open the RNA merged and filtered
adata = sc.read_h5ad(snakemake.output.merged_rna_anndata)

# Cell types
cell_types = ['Neurons', 'Astrocytes', 'Endothelial cells', 'Microglia', 'Oligodendrocyte progenitor cells', 'Oligodendrocytes', 'Pericytes']

markers = pd.read_csv('/data/CARD_singlecell/SN_atlas/input/PanglaoDB_markers_27_Mar_2020.tsv', delimiter='\t')
brain_markers = markers[
    (markers['cell type'].isin(cell_types)) & 
    (markers['canonical marker'] ==1) & 
    (markers['species'].isin(['Mm Hs', 'Hs']))]

sc.tl.leiden(adata, flavor='igraph', resolution=1, key_added='leiden')

dc.run_ora(
    mat=adata,
    net=brain_markers,
    source='cell type',
    target='official gene symbol',
    min_n=1,
    verbose=True,
    use_raw=False
)

acts = dc.get_acts(adata, obsm_key='ora_estimate')

acts_v = acts.X.ravel()
max_e = np.nanmax(acts_v[np.isfinite(acts_v)])
acts.X[~np.isfinite(acts.X)] = max_e

df = dc.rank_sources_groups(acts, groupby='leiden', reference='rest', method='t-test_overestim_var')

annotation_dict = df.groupby('group').head(1).set_index('group')['names'].to_dict()

adata.obs['cell_type'] = [annotation_dict[clust] for clust in adata.obs['leiden']]

adata.obs['cell_type'] = cell_assign

adata.write_h5ad(filename=snakemake.output.merged_rna_anndata, compression='gzip')