source(snakemake@params[['packages']])

future::plan('multicore', workers=snakemake@threads)
options(future.globals.maxSize=50 * 1000 * 1024^2)


metadata <- rbindlist(future.apply::future_lapply(snakemake@input[['seurat']], function(sample) {
    metadata <- copy(readRDS(sample)@meta.data); setDT(metadata, keep.rownames='cells')
}))

metadata <- na.omit(metadata)


counts_mixture <- mixtools::normalmixEM(metadata[, atac_counts_deviation_score], k=2)
counts_cutoff <- NormalMixCutoff(counts_mixture, lower=min(counts_mixture$mu), upper=max(counts_mixture$mu))


metadata[, `:=` (
    project=snakemake@params[['project_name']],
    atac_counts_qc_status=fifelse(atac_counts_deviation_score < counts_cutoff, 'PASS', 'FAIL')
)]


fwrite(x=metadata, file=snakemake@output[['metadata']])

