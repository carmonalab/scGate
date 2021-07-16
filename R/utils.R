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
  mm.genes <- unique(unlist(scGate_DB$mouse$Tcell@markers))
  hs.genes <- unique(unlist(scGate_DB$human$Tcell@markers))
  
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
    if (verbose) {
      mess <- "Running on single core... You can speed-up computation by increasing number of cores ('ncores' parameter)"
      message(mess)
    }
    
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
      vec[i] <- vec[i-1]
    }
  }
  return(vec)
}

scGate_helper <- function(data, gating.model=NULL, max.impurity=0.5,
                   max.impurity.decay=0,        
                   ndim=30, resol=3, assay="RNA",
                   sd.in=3, sd.out=7, species="mouse",
                   chunk.size=1000, maxRank=1500,
                   max.iterations=10, min.cells=100, stop.iterations=0.01,
                   additional.signatures=NULL,
                   genes.blacklist="Tcell.blacklist", 
                   skip.normalize=FALSE, ncores=1,
                   return_signature_scores=TRUE, verbose=FALSE, quiet=FALSE) {
  
  #Get data stored in trained model
  if (! class(gating.model) == "scGate_Model") {
    stop("Background model must be a scGate_Model object")
  }
  markers <- gating.model@markers
  pos.celltypes <- gating.model@positive_celltypes
  model <- gating.model@bg_model
  method <- gating.model@scoring_method
  
  if (!is.null(genes.blacklist)) {
    if (length(genes.blacklist)==1 && genes.blacklist == "Tcell.blacklist") {  #Default
      if (species=="human") {
        genes.blacklist <- genes.blacklist.Hs
      } else {
        genes.blacklist <- genes.blacklist.Mm 
      }
    }  
    if (is.list(genes.blacklist)) {
      genes.blacklist <- unlist(genes.blacklist)
    }
    genes.blacklist <- unique(genes.blacklist)
  }
  
  #Check that markers and threshold names correspond
  marker.names.pass <- intersect(names(markers), rownames(model))
  
  if (length(marker.names.pass) ==0) {
    mess <- "Error. Could not match markers with threshold file"
    stop(mess)
  }
  markers <- markers[marker.names.pass]
  
  #Check selected cell type(s), and autoexpand if required
  sign.names <- names(markers)
  celltype.pass <- check_selected_celltypes(pos.celltypes, db=sign.names, autocomplete=F)
  celltype.nopass <- setdiff(sign.names, celltype.pass)
  
  if (length(celltype.pass)<1) {
    mess <- sprintf("Could not find selected cell types in marker list")
    stop(mess)
  }  
  if (!quiet) {
    mess <- paste(celltype.pass, collapse=", ")
    message(sprintf("--- Filtering for target cell type(s): %s", mess))
  }
  
  #Get Zscores
  scores <- get_CTscores(obj=data, markers.list=markers, method=method, chunk.size=chunk.size, ncores=ncores,
                         bg=model, z.score=TRUE, maxRank=maxRank, additional.signatures=additional.signatures)
  
  #START FILTERING
  all_ind <- seq_along(rownames(scores))
  tags <- rep("pass", nrow(scores))
  names(tags) <- rownames(scores)
  
  #First, only keep cells with minimum expression levels for selected signature(s)
  positive_select <- c()
  for (sig in celltype.pass) {
    sig.zscore <- paste0(sig,"_Zscore")
    
    sd.use <- sd.in
    min.sd.in <- model[sig,"mean"]/model[sig,"sd"]-1e-3 # minimum z-score for positive scores
    if (sd.use > min.sd.in) {
      sd.use <- min.sd.in
      if(verbose) message(sprintf("Readjusting sd.in for %s to %.4f to avoid negative thresholds", sig, sd.use)) 
    }
    positive_select <- c(positive_select, which(scores[,sig.zscore] >= -sd.use))  # Z.score threshold for desired cell types
  }
  positive_select <- unique(positive_select)
  to.filter <- setdiff(all_ind, positive_select)
  
  #Then, filter out cells with excessive levels on undesired signatures
  negative_select <- c()
  
  for (sig in celltype.nopass){
    sig.zscore <- paste0(sig,"_Zscore")
    if(! sig %in% celltype.pass) {
      negative_select <- c(negative_select, which(scores[,sig.zscore] > sd.out))   # Z.score threshold for contaminants
    }
  }
  negative_select <- unique(negative_select)
  
  sig.zscores.out <- paste0(celltype.nopass,"_Zscore") 
  if (length(negative_select)>0) {
    max.index <- apply(scores[negative_select, sig.zscores.out], MARGIN = 1, which.max)
    tags[negative_select] <- sig.zscores.out[max.index]
  }
  
  unknown <- setdiff(to.filter, negative_select)
  tags[unknown] <- "Unknown"
  
  to.filter.ID <- names(tags[!tags=="pass"])

  #The vector 'labs' will contain the identification as "Pure" or "Impure" cell
  labs <- rep("Pure", nrow(scores))
  names(labs) <- rownames(scores)
  labs.ann <- labs
  tot.cells <- length(labs)
  q <- data
  
  #Start iterations
  for (iter in 1:max.iterations) {
    
    filterCells.this <- to.filter.ID[to.filter.ID %in% Cells(q)]
    imp.thr <- max.impurity - (iter-1)*max.impurity.decay
    
    ndim.use <- ifelse(ncol(q)<300, 5, ndim)  #with very few cells, reduce dimensionality
    
    ##High resolution clustering to detect dense regions of undesired cell types
    if (!skip.normalize) {
      q <- NormalizeData(q, verbose = FALSE)
    }
    q <- FindVariableFeatures(q, selection.method = "vst", nfeatures = 500, verbose = FALSE)
    q@assays[[assay]]@var.features <- setdiff(q@assays[[assay]]@var.features, genes.blacklist)
    q <- ScaleData(q, verbose=FALSE)
    q <- RunPCA(q, features = q@assays[[assay]]@var.features, verbose = FALSE)
    q <- FindNeighbors(q, reduction = "pca", dims = 1:ndim.use, k.param = 5, verbose=FALSE)
    q  <- FindClusters(q, resolution = resol, verbose = FALSE)
    q$clusterCT <- q@active.ident
    
    m <- q@meta.data
    impure.freq <- tapply(row.names(m), m$clusterCT,function(x) {mean(x %in% filterCells.this)})
    
    filterCluster <- names(impure.freq)[impure.freq > imp.thr]
    n_rem <- sum(q$clusterCT %in% filterCluster)
    frac.to.rem <- n_rem/tot.cells
    
    mess <- sprintf("----- Iter %i - max.impurity=%.3f\n-- Detected %i non-pure cells for selected signatures (%.2f%% of total cells)",
                    iter, imp.thr, n_rem, 100*frac.to.rem)
    if (verbose) {
      message(mess)
    }
    
    for (cl in filterCluster) {
       ids <- Cells(q)[which(m$clusterCT == cl)]
       tab <- table(tags[ids])
       tab <- tab[names(tab) != "pass"]
       ctype.out <- names(sort((tab), decreasing = T))[1]
       labs.ann[ids] <- gsub("_Zscore",replacement="", ctype.out)
    }
    q$is.pure <- ifelse(q$clusterCT %in% filterCluster,"Impure","Pure")
    
    labs[colnames(q)] <- q$is.pure
    
    #Subset on pure cells
    if (frac.to.rem < 1) {
       q <- subset(q, subset=is.pure=="Pure")
    }
    if (frac.to.rem < stop.iterations | iter>=max.iterations | ncol(q)<min.cells | frac.to.rem==1) {
      #Return clusters and active idents for easy filtering
      n_rem <- sum(labs=="Impure")
      frac.to.rem <- n_rem/tot.cells
      mess <- sprintf("scGate: Detected %i non-pure cells for selected signatures - %.2f%% cells marked for removal (active.ident)",
                      n_rem, 100*frac.to.rem)
      if (!quiet) {
        message(mess)
      }
      if (n_rem == tot.cells) {
        message("Warning: all cells were removed by scGate. Check filtering thresholds or input data")
      }
      
      if (return_signature_scores) {
        data <- AddMetaData(data, scores)
      }
      
      data$is.pure <- labs
      data$scGate.annotation <- labs.ann
      
      Idents(data) <- data$is.pure
      return(data)
    }
  }  
  return()
}

