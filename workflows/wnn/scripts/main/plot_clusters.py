import muon
import scanpy


scanpy.settings.verbosity = 1
scanpy.settings.n_jobs = snakemake.threads
scanpy.settings.figdir = snakemake.params.plot_directory # type: ignore
scanpy.settings.set_figure_params(dpi=300, fontsize=10, dpi_save=300, format='pdf', figsize=('5', '4')) # type: ignore


mode = snakemake.params.modalities
cluster_mapper = snakemake.params.type_mapper


if mode == 'wnn':
    mdata = muon.read_h5mu(snakemake.input.muon_object)
else:
    mdata = muon.read_anndata(snakemake.input.muon_object, mod=mode)

cluster_to_celltype = {cluster: celltype for celltype, clusters in cluster_mapper.items() for cluster in clusters}


mdata.obs['celltype'] = mdata.obs['leiden'].map(cluster_to_celltype)

plot_name = ''.join(['_', mode, '_celltype.pdf'])


muon.pl.embedding(
    mdata, color='celltype', frameon=False, 
    basis='X_wnn_umap', save=plot_name, title=mode.toupper()
)

