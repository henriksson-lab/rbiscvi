

# for dense arrays
# 'encoding-type': 'array' and 'encoding-version': '0.2.0'


# https://anndata.readthedocs.io/en/latest/fileformat-prose.html
# https://www.bioconductor.org/packages/devel/bioc/vignettes/rhdf5/inst/doc/rhdf5.html

######### Could support both Zarr and hdf5 with fairly little work



###############################################
#' @export
setClass("BiscviAnndata", slots=list(
  fname="character"
  )
) 




###############################################
#' 
#' 
BiscviOpen <- function(fname) {
  if(!file.exists(fname)) {
    rhdf5::h5createFile(fname)
  }
  
  new(
    "BiscviAnndata", 
    fname=fname
  )
}



###############################################
#' Delete previous object if it exists
BiscviDeleteIfExists <- function(fp, obname) {
  if(BiscviEntryExists(fp, obname)) {
    #print("Deleting previous object")
    rhdf5::h5delete(fp@fname, obname)
  }
}

###############################################
#' Check if entry exists. This can likely be made faster
BiscviEntryExists <- function(fp, obname) {
  contents <- rhdf5::h5ls(fp@fname)
  contents$group[contents$group=="/"] <- ""
  obname %in% paste0(contents$group, "/", contents$name)
}
  



###############################################
#' Store xx
#' 
BiscviStoreDataframe  <- function(fp, df, obname) {

  #Create the group  
  BiscviDeleteIfExists(fp, obname)
  gr <- rhdf5::h5createGroup(fp@fname, obname)
  attr(gr, "_index") <- "name of some col"
  attr(gr, "encoding-type") <- "dataframe"
  attr(gr, "encoding-version") <- "0.2.0"
  attr(gr, "column-order") <- colnames(df)

  #Store each row separately
  for(i in 1:ncol(df)) {
    col_obname <- paste0(obname,"/",colnames(df)[i])
    BiscviStoreVec(fp, df[,i], col_obname)
  }
}



###############################################
#' Store xx
#' 
BiscviStoreVec <- function(fp, dat, obname) {

  #Ensure that strings are always stored as factors
  if(is.character(dat)) {
    dat <- factor(dat)
  }
  
  if(is.factor(dat)) {
    ## Store a factor vector
    print("factor")
    print(head(dat))
    
    #Create the group
    gr <- rhdf5::h5createGroup(fp@fname, obname)
    attr(gr, "encoding-type") <- "categorical"
    attr(gr, "encoding-version") <- "0.2.0"
    attr(gr, "ordered") <- FALSE #not sure what for
    
    #Store the values
    towrite <- as.integer(dat)
    attr(towrite, "shape") <- length(towrite) #needed? and type?
    rhdf5::h5write(towrite, fp@fname, paste0(obname,"/codes"), level=0, write.attributes=TRUE)

    #Store the factor levels
    towrite <- levels(dat)
    print(towrite)
    attr(towrite, "shape") <- length(towrite)  #needed? and type?
    
    name_categories <- paste0(obname,"/categories")
    rhdf5::h5createDataset(
      fp@fname, name_categories, 
      dims = length(towrite), #chunk=1000, 
      storage.mode=storage.mode("character"), size=NULL
    )
    rhdf5::h5write(
      towrite, 
      fp@fname, 
      name_categories, 
      level=0, 
      write.attributes=TRUE
    )
    
  } else {
    ## Store a regular vector

    print("values")
    print(head(dat))
    
    attr(dat, "encoding-type") <- "array"
    attr(dat, "encoding-version") <- "0.2.0"
    rhdf5::h5write(dat, fp@fname, obname, level=0, write.attributes=TRUE)
    
  }
}



###############################################
#' Store metadata
#' 
#' @param adata A Seurat object
#' @param fname File name with .h5-ending
#' 
#' @export
BiscviStoreMeta  <- function(adata, fp) {
  
  BiscviStoreDataframe(fp, adata@meta.data, "/obs")
  #### cannot store data in var!!! different features per layer
  #rownames(adata) 
  #BiscviStoreDataframe(fp, adata@meta.data, "/var")
}





###############################################
#' Store counts
#' 
#' @param adata A Seurat object
#' @param fname File name with .h5-ending
#' @param assayName Name of assay to store
#' 
#' @export
BiscviStoreCounts <- function(adata, fp, assayName) {
  
  compress_algo <- "NONE"

  #Different treatment based on Seurat version
  the_assay <- adata@assays[[assayName]]
  #the_assay <- adata@assays[["RNA"]]
  #class(the_assay)
  #exists(the_assay@data)
  
  #Get the counts to store  
  if(inherits(the_assay,"Assay5")) {
    spmat <- the_assay@layers$data
  } else {
    spmat <- the_assay@data
  }
  
  #spmat <- adata@assays[[assayName]]@data

  #Ensure layers group exists
  if(!BiscviEntryExists(fp, "/counts")) {
    rhdf5::h5createGroup(fp@fname,"/counts")
  }

  #Create a new group for this type of counts
  groupname <- paste0("/counts/", assayName)
  BiscviDeleteIfExists(fp, groupname)
  gr <- rhdf5::h5createGroup(fp@fname,groupname)
  attr(gr, "encoding-type") <- "csr_matrix"  # or csc_matrix
  attr(gr, "encoding-version") <- "0.1.0"
  attr(gr, "shape") <- dim(spmat)
  
  #Store data without compression to enable fast subsetting
  name_data <- paste0(groupname,"/data")
  rhdf5::h5createDataset(fp@fname, name_data, dims = length(spmat@x), filter=compress_algo, chunk=1000, storage.mode=storage.mode(spmat@x))
  rhdf5::h5write(spmat@x, fp@fname, name_data, level=0, write.attributes=TRUE)
  
  name_indices <- paste0(groupname,"/indices")
  rhdf5::h5createDataset(fp@fname, name_indices, dims = length(spmat@i), filter=compress_algo, chunk=1000, storage.mode=storage.mode(spmat@i))
  rhdf5::h5write(spmat@i, fp@fname, name_indices, level=0, write.attributes=TRUE)
  
  name_indptr <- paste0(groupname,"/indptr")
  rhdf5::h5createDataset(
    fp@fname, name_indptr, 
    dims = length(spmat@p), 
    filter=compress_algo, 
    chunk=1000, 
    storage.mode=storage.mode(spmat@p)
  )
  rhdf5::h5write(spmat@p, fp@fname, name_indptr, level=0, write.attributes=TRUE)
  
  #Store the name of features. This replaces "var" in regular anndata
  feature_names <- rownames(the_assay)
  name_feature_names <- paste0(groupname,"/feature_names")
  rhdf5::h5createDataset(
    fp@fname, name_feature_names, 
    dims = length(feature_names), filter=compress_algo, chunk=1000, 
    storage.mode=storage.mode("character"), size=NULL
  )
  rhdf5::h5write(feature_names, fp@fname, name_feature_names, level=0, write.attributes=TRUE)
}


#BiscviStoreCounts(adata, fp, assayName = "RNA")


#library(rhdf5)
BiscviStoreReduction <- function(adata, fp, redName) {
  
  #Ensure reductions group exists
  if(!BiscviEntryExists(fp, "/reductions")) {
    rhdf5::h5createGroup(fp@fname,"/reductions")
  }
  
  
  #redName <- "kraken_umap"
  compress_algo <- "NONE"
  
  emb <- as.matrix(adata@reductions[[redName]]@cell.embeddings)

  #Prepare to overwrite if needed
  groupname <- paste0("/reductions/", redName)
  BiscviDeleteIfExists(fp, groupname)
    
  #name_feature_names <- paste0(groupname,"/feature_names")
  rhdf5::h5createDataset(fp@fname, groupname, dims = dim(emb), filter=compress_algo, chunk=c(1000,2), storage.mode=storage.mode(emb))
  rhdf5::h5write(emb, fp@fname, groupname, level=0)
  
}  


