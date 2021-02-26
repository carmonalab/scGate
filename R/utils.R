check_CTmarkers <- function(obj, markers.list, min.gene.frac=0.5, verbose=TRUE) {
  markers.list.exp <- unlist(lapply(markers.list,function(x){mean(x %in% rownames(obj))}))
  
  markers.list.pass <- markers.list[markers.list.exp >= min.gene.frac]
  markers.list.pass <- lapply(markers.list.pass,function(x){x[x %in% rownames(obj)]})
  
  skip <- names(markers.list.exp)[markers.list.exp < min.gene.frac]
  if (verbose==TRUE & length(skip) > 0) {
    skip <- paste(skip, collapse=", ")
    warning(paste0("Warning: not enough genes to evaluate the following signatures: ",skip))
  }
  return(markers.list.pass)
}

get_CTscores <- function(obj, markers.list, rm.existing=TRUE, bg=NULL) {
  if (rm.existing) {  #remove existing CTfilter scores - necessary because AddModuleScore won't overwrite them
    index.rm <- grep("_CTfilter",colnames(obj@meta.data), perl=T)
    if (length(index.rm)>0) {
      obj@meta.data <- subset(obj@meta.data, select = -index.rm)
    }
  }
  
  obj <- suppressWarnings(AddModuleScore(obj, features = markers.list, name="CTfilterScore"))
  sign.names <- paste0(names(markers.list),"_CTfilter")
  
  index.start <- ncol(obj@meta.data)-length(markers.list)+1
  index.end <- ncol(obj@meta.data)
  
  names(obj@meta.data)[index.start:index.end] <- sign.names
  
  #Convert to Z-score
  for (j in index.start:index.end) {
     sign <- colnames(obj@meta.data)[j]
     obj@meta.data[,j] <- (obj@meta.data[,j] - bg[sign,"mean"])/bg[sign,"sd"]
  }
  return(obj)
}

vectorize.parameters <- function(par, lgt=1) {
  vec <- vector(mode="numeric", lgt)
  for (i in 1:length(vec)) {
    if (i <= length(par)) {
      vec[i] <- par[i]
    } else {
      vec[i] <- par[i-1]
    }
  }
  return(vec)
}