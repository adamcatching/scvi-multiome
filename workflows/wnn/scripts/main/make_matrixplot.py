import muon
import scanpy
import pandas


scanpy.settings.verbosity = 1
scanpy.settings.n_jobs = snakemake.threads
scanpy.settings.figdir = snakemake.params.plot_directory # type: ignore
scanpy.settings.set_figure_params(dpi=300, fontsize=10, dpi_save=300, format='png', figsize=('6', '4')) # type: ignore


df = pandas.read_csv(snakemake.input.markers, index_col=0) # type: ignore
adata = muon.read_anndata(snakemake.input.muon_object, mod='rna') # type: ignore

markers = {
    cell_type: df[df[cell_type] == 1].index.tolist()
    for cell_type in df.columns
}

scanpy.tl.dendrogram(adata, groupby=snakemake.params.group, use_rep='X_scvi', cor_method='pearson')   # type: ignore


scanpy.pl.matrixplot(
    adata, var_names=markers, dendrogram=True, standard_scale='var', cmap='viridis', 
    save=''.join([snakemake.params.group, '.png']), groupby=snakemake.params.group, use_raw=True # type: ignore
)
