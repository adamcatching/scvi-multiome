source(snakemake@params[['packages']])


object <- readRDS(snakemake@input[['seurat']])

assay <- 'RNA'
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


object <- object %>% 

    NormalizeData(verbose=FALSE) %>% 
    ScaleData(features=unique(unlist(cc.genes)), verbose=FALSE) %>% 
    CellCycleScoring(s.features=s.genes, g2m.features=g2m.genes) %>% 

    DietSeurat(assays=assay)


saveRDS(object, snakemake@output[['seurat']])