import snapatac2 as snap
import scvi
import anndata as ad
import scipy
import pandas as pd

# Import sample metadata
samples = pd.read_csv('/data/CARD_singlecell/SN_atlas/input/SNsamples.csv')

# From the metadata csv extract the ID of the sample and the batch folder it is saved under
sample_names = samples['Sample_ID'].to_list()
batches = samples['Use_batch'].to_list()
# Implement the file structure heirarchy here to use the correct sample and 
sample_locs = [f'/data/CARD_singlecell/Brain_atlas/SN_Multiome/batch{batches[i]}/Multiome/{sample_names[i]}-ARC/outs/' for i in range(len(batches))]

# Define the location of the fragment files to use for each snapATAC2 run
fragment_files = [f'{sample_locs}atac_fragments.tsv.gz' for sample_locs in sample_locs]

# Save the individual snapATAC2 file locations and their file names
outputs = []
names = []

# Create a list of the file names and the location to save the files
for fl in fragment_files:
    name = fl.split('/')[-3]
    name = name.split('-')[0]
    temp_dir = '/data/CARD_singlecell/SN_atlas/data/samples/' + name + '/01_atac_anndata.h5ad'
    names.append(name)
    outputs.append(temp_dir)

# Import the files as a list of AnnData objects (these are not saved in memory)
adatas = snap.pp.import_data(
    fragment_files, 
    file=outputs, 
    chrom_sizes=snap.genome.hg38.chrom_sizes,
    min_num_fragments=1000,
    sorted_by_barcode=False,
    n_jobs=60)

# Generate the equaly spaced tiles 
snap.pp.add_tile_matrix(adatas, bin_size=5000)

# Select variable features
snap.pp.select_features(adatas, n_features=250000)

# Create the AnnDataSet object
adataset = snap.AnnDataSet(
    adatas=[(names[i], f) for i, f in enumerate(adatas)],
    filename='/data/CARD_singlecell/SN_atlas/data/atlas/02_merged_anndata_atac.h5ad'
)

# Select variable features on the new combined AnnData object
snap.pp.select_features(adataset, n_features=250000)

""" This isn't quite working yet"""
#snap.pp.mnc_correct(adataset, batch='sample', key_added='X_spectral')

#snap.tl.spectral(adataset)
#snap.pp.knn(adataset)
#snap.tl.leiden(adataset)
#snap.tl.umap(adataset, n_comps=2)

# Convert to an AnnData atlas object (this takes a lot of memory, ~1.5 Tb on the SN atlas)
adata = adataset.to_adata()

# For redundancy, save the AnnData object
adata.write_h5ad('/data/CARD_singlecell/SN_atlas/data/atlas/01_merged_anndata_atac.h5ad', compression='gzip')
