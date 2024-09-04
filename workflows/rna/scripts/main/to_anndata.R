
sceasy::convertFormat(
    obj=readRDS(snakemake@input[['seurat']]),
    outFile=snakemake@output[['anndata']],
    assay='RNA', main_layer='counts', 
    from='seurat', to='anndata', 
    drop_single_values=FALSE
)

