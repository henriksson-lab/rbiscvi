
if(FALSE){
  ############# A scMetaG dataset ###################
  adata <- readRDS("/home/mahogny/github/rbiscvi/adata.RDS")
  
  fname <- "/home/mahogny/github/rbiscvi/counts.biscvi5"  #biscvi5 is the name of the file format
  fp <- BiscviOpen(fname)
  
  #df <- adata@meta.data
  #BiscviEntryExists(fp, "/layers")
  #BiscviStoreDataframe(fp, df, "/meta")
  
  BiscviStoreMeta(adata, fp)
  BiscviStoreCounts(adata, fp, assayName = "RNA")
  BiscviStoreReduction(adata, fp, redName = "umap")
  BiscviStoreReduction(adata, fp, redName = "pca")
  BiscviStoreReduction(adata, fp, redName = "kraken_umap")
  
  
  # /home/mahogny/github/rbiscvi/counts.biscvi5
  
  rhdf5::h5ls(fname, all=TRUE) # H5T_STRING
  
}

if(FALSE) {
  h5f <- rhdf5::H5Fopen(fname)
}
