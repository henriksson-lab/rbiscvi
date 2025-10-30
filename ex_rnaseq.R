


if(FALSE){
  ############# A scMetaG dataset ###################
  #adata <- readRDS("/home/mahogny/github/rbiscvi/adata.RDS")
  adata <- LoadSeuratRds("~/mystore/cartdata/data/publicDengDatIntegrated.rds")
  
  fname <- "/home/mahogny/github/rbiscvi/counts_deng.biscvi5"  #biscvi5 is the name of the file format
  fp <- BiscviOpen(fname)
  
  BiscviStoreMeta(adata, fp)
  BiscviStoreCounts(adata, fp, assayName = "RNA")
  BiscviStoreReduction(adata, fp, redName = "umap")
  BiscviStoreReduction(adata, fp, redName = "pca")
  
  # /home/mahogny/github/rbiscvi/counts.biscvi5
  
  rhdf5::h5ls(fname, all=TRUE) # H5T_STRING
  
}
