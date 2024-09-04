source(snakemake@params[['packages']])

object <- readRDS(snakemake@input[['seurat']])
        
object <- object %>% subset(
    
    subset=
            nFeature_ATAC > 300 & nucleosome_signal < 3 & frip > 0.1 &
            TSS.enrichment > 3 & nCount_ATAC %between% c(300, 50000)
    )

saveRDS(object, snakemake@output[['seurat']])

macs_path <- '/data/abbass2/miniforge3/envs/macs3/bin/macs3'

peaks <- object %>% 

            CallPeaks(assay='ATAC', macs2.path=macs_path) %>% 
            GenomeInfoDb::keepStandardChromosomes(pruning.mode='coarse') %>% 
            IRanges::subsetByOverlaps(ranges=blacklist_hg38_unified, invert=TRUE)


saveRDS(peaks, snakemake@output[['peaks']])

