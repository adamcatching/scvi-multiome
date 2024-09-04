source(snakemake@params[['packages']])

future::plan('multicore', workers=snakemake@threads)
options(future.globals.maxSize=50 * 1000 * 1024^2)


object <- readRDS(snakemake@input[['seurat']])
peaks.use <- readRDS(snakemake@input[['merged_peaks']])


annotations <- get_annotations()
    
dataset <- Project(object)
batch <- unique(object@meta.data$batch)
cohort <- unique(object@meta.data$cohort)

data_root <- paste0(data_path, toupper(cohort), '_multiome/', batch, '/Multiome/', dataset)


fragments <- Fragments(object[['ATAC']])[[1]]
cells.use <- GetFragmentData(fragments, slot='cells')               
fragpath <- paste0(data_root, '/outs/atac_fragments.tsv.gz')


object <- object %>% subset(cells=cells.use)
counts <- FeatureMatrix(cells=cells.use, features=peaks.use, fragments=fragments, process_n=5000)
object[['ATAC']] <- CreateChromatinAssay(counts=counts, fragments=fragpath, annotation=annotations)


saveRDS(object, snakemake@output[['seurat']])
