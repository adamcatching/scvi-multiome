import muon
import scipy
import scanpy


rna = scanpy.read_h5ad(snakemake.input.rna)
atac = scanpy.read_h5ad(snakemake.input.atac)


rna.raw = rna; atac.raw = atac

rna.layers['counts'] = scipy.sparse.csr_matrix(rna.X.copy())
atac.layers['counts'] = scipy.sparse.csr_matrix(atac.X.copy())


mdata = muon.MuData({'rna': rna, 'atac': atac})


muon.pp.intersect_obs(mdata)

mdata.write_h5mu(filename=snakemake.output.muon_object, compression='gzip') # type: ignore

