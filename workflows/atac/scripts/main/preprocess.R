source(snakemake@params[['packages']])

assay <- 'ATAC'

hub <- AnnotationHub::AnnotationHub()
reference <- hub[['AH75011']]

dataset <- snakemake@params[['sample']]
m <- fread(snakemake@input[['metadata']])[sample %chin% dataset][, cells := tstrsplit(cells, '_')[[1]]]

batch <- m[, unique(batch)]; cohort <- m[, unique(cohort)]

data_root <- paste0(data_path, toupper(cohort), '_multiome/', batch, '/Multiome/', dataset)


frag.file <- paste0(data_root, '/outs/atac_fragments.tsv.gz')
counts <- Read10X_h5(paste0(data_root, '/outs/raw_feature_bc_matrix.h5'))$Peaks


counts <- counts[, colnames(counts) %in% m[, cells]]

frag_stats <- CountFragments(fragments=frag.file, verbose=FALSE)
setDT(frag_stats); setnames(frag_stats, old='CB', new='cells')


metadata <- frag_stats[m, on='cells']
setDF(metadata, rownames=metadata[, cells]); metadata$cells <- NULL

 
object <- CreateSeuratObject(
    assay=assay, project=dataset, meta.data=metadata,
    counts=BuildChromiumAssay(counts=counts, enDB=reference, frag.file=frag.file)
)  


object <- object %>% 

            TSSEnrichment(assay=assay, verbose=FALSE) %>%
            NucleosomeSignal(assay=assay, verbose=FALSE) %>% 
            FRiP(assay=assay, total.fragments='frequency_count', col.name='frip', verbose=FALSE)


object$atac_counts_deviation_score <- median_deviation_score(object$nCount_ATAC)


saveRDS(object, snakemake@output[['seurat']])
