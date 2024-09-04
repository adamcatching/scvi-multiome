source(snakemake@params[['packages']])

future::plan('multicore', workers=snakemake@threads)
options(future.globals.maxSize=50 * 1000 * 1024^2)


peaks.list <- future.apply::future_lapply(snakemake@input[['peaks']], readRDS)

peaks.use <- suppressWarnings(reduce(x=unlist(as(peaks.list, 'GRangesList'))))

saveRDS(peaks.use, snakemake@output[['merged_peaks']])