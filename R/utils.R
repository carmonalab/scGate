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

get_CTscores <- function(obj, markers.list, rm.existing=TRUE, bg=NULL, raw.score=FALSE, chunk.size=1000, ncores=1) {
  
  if (rm.existing) {  #remove existing CTfilter scores - necessary because AddModuleScore won't overwrite them
    index.rm <- grep("_CTfilter",colnames(obj@meta.data), perl=T)
    if (length(index.rm)>0) {
      obj@meta.data <- subset(obj@meta.data, select = -index.rm)
    }
  }
  
  #obj <- suppressWarnings(AddModuleScore(obj, features = markers.list, name="CTfilterScore"))
  obj <- run_AUCell(obj, markers.list=markers.list, chunk.size=chunk.size, ncores=ncores)
  
  if (raw.score) {
    return(obj) 
  }
  
  #Convert to Z-score
  for (j in index.start:index.end) {
    sign <- colnames(obj@meta.data)[j]
    obj@meta.data[,j] <- (obj@meta.data[,j] - bg[sign,"mean"])/bg[sign,"sd"]
  }
  return(obj)
}

run_AUCell <- function(obj, markers.list, chunk.size=1000, ncores=1) {
  
   if (ncores>1) {
     require(future.apply)
   }
  
   assay <- DefaultAssay(obj)
   
   #Split into manageable chunks
   split.data <- split_data.matrix(matrix=obj@assays[[assay]]@data, chunk.size=chunk.size)
   
   #Parallelize?
   if (ncores>1) {
     plan(future::multisession(workers=ncores))
     
     meta.list <- future_lapply(
       X = 1:length(x = split.data),
       FUN = function(i) {
         cells_rankings <- AUCell_buildRankings(split.data[[i]], plotStats=F, verbose=F)
         cells_AUC <- AUCell_calcAUC(markers.list, cells_rankings)

         new.meta <- as.data.frame(t(getAUC(cells_AUC)))
         colnames(new.meta) <- paste0(colnames(new.meta),"_CTfilter")
         return(new.meta)
       },
       future.seed = 123
     )
     plan(strategy = "sequential")
   } else {
     meta.list <- lapply(
       X = 1:length(x = split.data),
       FUN = function(i) {
         cells_rankings <- AUCell_buildRankings(split.data[[i]], plotStats=F, verbose=F)
         cells_AUC <- AUCell_calcAUC(markers.list, cells_rankings)
         
         new.meta <- as.data.frame(t(getAUC(cells_AUC)))
         colnames(new.meta) <- paste0(colnames(new.meta),"_CTfilter")
         return(new.meta)
       } )
   }
   
   meta.merge <- Reduce(rbind, meta.list)
   
   obj <- AddMetaData(obj, meta.merge)
   return(obj)
}

split_data.matrix <- function(matrix, chunk.size=1000) {
   ncols <- dim(matrix)[2]
   nchunks <- ncols %/% chunk.size + 1
   
   split.data <- list()
   for (i in 1:nchunks) {
      min <- 1 + (i-1)*chunk.size
      max <- min(i*chunk.size, ncols)
      split.data[[i]] <- matrix[,min:max]
   }
   return(split.data)
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