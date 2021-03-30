
#Calculate CTfilter scores
get_CTscores <- function(obj, markers.list, rm.existing=TRUE, bg=NULL, z.score=FALSE,
                         method=c("UCell","AUCell","ModuleScore"), chunk.size=1000, ncores=1, maxRank=1500) {
  
  method.use <- method[1]
  
  if (rm.existing) {  #remove existing CTfilter scores
    index.rm <- grep("_CTfilter|_Zscore",colnames(obj@meta.data), perl=T)
    if (length(index.rm)>0) {
      obj@meta.data <- subset(obj@meta.data, select = -index.rm)
    }
  }
  
  if (method.use == "UCell") {
     obj <- UCell::AddModuleScore_UCell(obj, features=markers.list, chunk.size=chunk.size, ncores=ncores, maxRank=maxRank, name="_CTfilter")
  } else if (method.use == "AUCell") {
     obj <- AddModuleScore_AUCell(obj, features=markers.list, chunk.size=chunk.size, ncores=ncores)
  } else if (method.use == "ModuleScore") {
     obj <- suppressWarnings(AddModuleScore(obj, features = markers.list, name="CTfilterScore"))
  } else {
     stop("Please give a valid method for signature scoring (see 'method' parameter)")
  }

  if (!z.score) {
    return(obj) 
  }
  
  #Convert to Z-score
  cols <- paste0(names(markers.list),"_CTfilter")
  zcols <- paste0(names(markers.list),"_Zscore")
  
  for (i in seq_along(cols)) {
    sig <- cols[i]
    zsig <- zcols[i]
    obj@meta.data[,zsig] <- (obj@meta.data[,sig] - bg[sig,"mean"])/bg[sig,"sd"]
  }
  return(obj)
}




AddModuleScore_AUCell <- function(obj, features, chunk.size=1000, ncores=1) {
  
  require(AUCell)
  assay <- DefaultAssay(obj)
  
  #Split into manageable chunks
  split.data <- split_data.matrix(matrix=obj@assays[[assay]]@data, chunk.size=chunk.size)
  
  #Parallelize?
  if (ncores>1) {
    plan(future::multisession(workers=future_param_ncores))
    
    meta.list <- future_lapply(
      X = split.data,
      FUN = function(x) {
        cells_rankings <- AUCell_buildRankings(x, plotStats=F, verbose=F)
        cells_AUC <- AUCell_calcAUC(features, cells_rankings)
        
        new.meta <- as.data.frame(t(getAUC(cells_AUC)))
        colnames(new.meta) <- paste0(colnames(new.meta),"_CTfilter")
        return(new.meta)
      },
      future.seed = future_param_seed
    )
    plan(strategy = "sequential")
    
  } else {
    meta.list <- lapply(
      X = split.data,
      FUN = function(x) {
        cells_rankings <- AUCell_buildRankings(x, plotStats=F, verbose=F)
        cells_AUC <- AUCell_calcAUC(features, cells_rankings)
        
        new.meta <- as.data.frame(t(getAUC(cells_AUC)))
        colnames(new.meta) <- paste0(colnames(new.meta),"_CTfilter")
        return(new.meta)
      } )
  }
  
  meta.merge <- Reduce(rbind, meta.list)
  obj <- Seurat::AddMetaData(obj, as.data.frame(meta.merge))
  
  return(obj)
}

split_data.matrix <- function(matrix, chunk.size=1000) {
   ncols <- dim(matrix)[2]
   nchunks <- (ncols-1) %/% chunk.size + 1
   
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