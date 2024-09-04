source(snakemake@params[['packages']])


object <- readRDS(snakemake@input[['seurat']])

Project(object) <- tstrsplit(Project(object), '_')[[2]]


object <- object %>% subset(
    
    subset=
            total_counts < 100000 & n_genes_by_counts > 300 & pct_counts_mt < 5 & pct_counts_rb < 5,
    cells=
            fread(snakemake@input[['metadata']])[
                sample %chin% Project(object) & 
                predicted_gmm_doublets %chin% 'singlet', cells
            ]  
    )


saveRDS(object, snakemake@output[['seurat']])
