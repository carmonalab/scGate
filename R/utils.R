
check_selected_celltypes <- function(celltype, db, autocomplete=T) {
  match <- vector()
  
  for (i in seq_along(celltype)) {
    regexp = paste0("^",celltype[i])
    if (!autocomplete) {  #exact match
      regexp = paste0(regexp, '$')   
    }
    
    this <- grep(regexp, db, value=T, perl=T)
    match <- c(match,this)
  }
  return(unique(match))
}

#Calculate scGate scores
get_CTscores <- function(obj, markers.list, rm.existing=TRUE, bg=NULL, z.score=FALSE,
                         method=c("UCell","ModuleScore"), chunk.size=1000, ncores=1, maxRank=1500) {
  
  method.use <- method[1]
  
  if (rm.existing) {  #remove existing scGate scores
    index.rm <- grep("_scGate|_Zscore",colnames(obj@meta.data), perl=T)
    if (length(index.rm)>0) {
      obj@meta.data <- subset(obj@meta.data, select = -index.rm)
    }
  }
  
  if (method.use == "UCell") {
     obj <- suppressWarnings(UCell::AddModuleScore_UCell(obj, features=markers.list, chunk.size=chunk.size, ncores=ncores, maxRank=maxRank, name="_scGate"))
  } else if (method.use == "AUCell") {
     obj <- suppressWarnings(AddModuleScore_AUCell(obj, features=markers.list, chunk.size=chunk.size, ncores=ncores, maxRank=maxRank))
  } else if (method.use == "ModuleScore") { ##TO DO: Rename output to make it compatible with other methods
     obj <- suppressWarnings(AddModuleScore(obj, features = markers.list, name="scGateScore"))
  } else {
     stop("Please give a valid method for signature scoring (see 'method' parameter)")
  }

  if (!z.score) {
    return(obj) 
  }
  
  #Convert to Z-score
  celltypes <- names(markers.list)
  
  for (ct in celltypes) {
    ct_z <- paste0(ct,"_Zscore")
    ct_s <- paste0(ct,"_scGate")
    obj@meta.data[,ct_z] <- (obj@meta.data[,ct_s] - bg[ct,"mean"])/bg[ct,"sd"]
  }
  return(obj)
}

AddModuleScore_AUCell <- function(obj, features, chunk.size=1000, ncores=1, maxRank=1500, name="_scGate") {
  
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
        cells_AUC <- AUCell_calcAUC(features, cells_rankings, aucMaxRank=maxRank)
        
        new.meta <- as.data.frame(t(getAUC(cells_AUC)))
        colnames(new.meta) <- paste0(colnames(new.meta), name)
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
        cells_AUC <- AUCell_calcAUC(features, cells_rankings, aucMaxRank=maxRank)
        
        new.meta <- as.data.frame(t(getAUC(cells_AUC)))
        colnames(new.meta) <- paste0(colnames(new.meta), name)
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