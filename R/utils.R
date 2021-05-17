#Match selected cell types with the DB of cell types
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

#Guess the species from the gene names of the query
detect_species <- function(query) {
  mm.genes <- unique(unlist(scGate_mouse_TIL_model@markers))
  hs.genes <- unique(unlist(scGate_human_TIL_model@markers))
  
  mm.intersect <- length(intersect(mm.genes, rownames(query)))/length(mm.genes)
  hs.intersect <- length(intersect(hs.genes, rownames(query)))/length(hs.genes)
  if (max(mm.intersect, hs.intersect)<0.2) {
    warning("More than 80% of genes not found in reference signatures...did you remove genes from the query data?")
  }
  if (mm.intersect>hs.intersect) {
    species <- "mouse"
  } else {
    species <- "human"
  }
  return(species)
}

#Calculate scGate scores
get_CTscores <- function(obj, markers.list, bg=NULL, z.score=FALSE, additional.signatures=NULL,
                         method=c("UCell","AUCell","ModuleScore"), chunk.size=1000, ncores=1, maxRank=1500) {
  
  method.use <- method[1]
  celltypes <- names(markers.list)
  
  if (is.list(additional.signatures) & length(additional.signatures)>0) {
     markers.list <- append(markers.list, additional.signatures)
  } 
  celltypes.add <- names(markers.list)
  
  #remove existing scGate scores
  index.rm <- grep("_scGate|_Zscore",colnames(obj@meta.data), perl=T)
  if (length(index.rm)>0) {
    obj@meta.data <- subset(obj@meta.data, select = -index.rm)
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
  #We only keep the new signature scores
  ct_s <- paste0(celltypes.add,"_scGate")
  scores <- obj@meta.data[,ct_s]
  
  if (!z.score) {
    return(scores) 
  }
  
  #Convert to Z-score
  for (ct in celltypes) {
    ct_z <- paste0(ct,"_Zscore")
    ct_s <- paste0(ct,"_scGate")
    scores[,ct_z] <- (scores[,ct_s] - bg[ct,"mean"])/bg[ct,"sd"]
  }
  return(scores)
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