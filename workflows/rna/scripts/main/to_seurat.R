source(snakemake@params[['packages']])

sceasy::convertFormat(
    obj=snakemake@input[['anndata']],
    outFile=snakemake@output[['seurat']],
    assay='RNA', main_layer='counts',
    from='anndata', to='seurat'
)



object <- readRDS(snakemake@output[['seurat']])

object$n_genes_by_counts <- object$nFeaturess_RNA_by_counts
object$nFeaturess_RNA_by_counts <- NULL


metadata <- copy(object@meta.data)
setDT(metadata, keep.rownames='cells')

metadata[, counts_deviation := median_deviation_score(total_counts)]

counts_deviation_score <- metadata[, counts_deviation]; names(counts_deviation_score) <- metadata[, cells]

object <- object %>% AddMetaData(metadata=counts_deviation_score, col.name='counts_deviation_score')


saveRDS(object, snakemake@output[['seurat']])