source(snakemake@params[['packages']])

future::plan('multicore', workers=snakemake@threads)
options(future.globals.maxSize=50 * 1000 * 1024^2)


metadata <- rbindlist(future.apply::future_lapply(snakemake@input[['seurat']], function(sample) {
    metadata <- copy(readRDS(sample)@meta.data); setDT(metadata, keep.rownames='cells')
}))


metadata <- na.omit(metadata)

doublet_mixture <- mixtools::normalmixEM(metadata[, doublet_score], k=2)
doublet_cutoff <- NormalMixCutoff(doublet_mixture, lower=min(doublet_mixture$mu), upper=max(doublet_mixture$mu))

counts_mixture <- mixtools::normalmixEM(metadata[, counts_deviation_score], k=2)
counts_cutoff <- NormalMixCutoff(counts_mixture, lower=min(counts_mixture$mu), upper=max(counts_mixture$mu))


metadata[, `:=` (
    project=snakemake@params[['project_name']],
    qc_status=fifelse(counts_deviation_score < counts_cutoff, 'PASS', 'FAIL'),
    predicted_gmm_doublets=fifelse(doublet_score < doublet_cutoff, 'singlet', 'doublet')
)]


fwrite(x=metadata, file=snakemake@output[['metadata']])

