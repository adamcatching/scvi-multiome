# scvi-multiome
Snakemake Workflow for using SCVI on scMultiome Datasets


### scVI

scVI (single-cell Variational Inference) is a scalable probabilistic framework designed to analyze single-cell RNA sequencing (scRNA-seq) data. Built on variational autoencoders, scVI offers powerful tools for batch effect correction, dimensionality reduction, and data integration across multiple samples or experiments. It enables efficient processing of large scRNA-seq datasets, supports flexible modeling of gene expression, and provides biologically interpretable latent spaces for downstream analysis such as clustering, differential expression, and cell type annotation.

### poissonVI

PoissonVI is an extension of the scVI framework specifically designed for single-cell ATAC sequencing (scATAC-seq) data. Unlike scVI, which models gene expression counts using a negative binomial distribution, PoissonVI leverages the Poisson distribution to more accurately model sparse, discrete peak accessibility counts characteristic of scATAC-seq data. This allows PoissonVI to provide more precise analysis of chromatin accessibility and better capture the underlying structure of epigenetic data for tasks such as clustering, dimensionality reduction, and identifying regulatory elements.

### WNN

WNN (Weighted Nearest Neighbors), implemented in the muon framework, is a powerful method for integrating multiple single-cell modalities such as RNA and ATAC sequencing. WNN leverages information from each modality by calculating separate nearest neighbor graphs for each data type and then weighting these graphs based on the relative contribution of each modality to a cell's overall identity. This approach allows for robust multimodal integration, enabling researchers to capture complementary information from both gene expression and chromatin accessibility, or other paired modalities, leading to improved cell type annotation, clustering, and interpretation of complex cellular states in single-cell multiome datasets.
