
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


#Calculate AUC as Mann–Whitney U statistic from a vector of ranks
# (len_sig*(len_sig+1))/2 is the total rank sum
# U1+U2=n1*n2 (n1,n2: sample sizes); U1/(n1*n2)=AUC
u_stat = function(rank_value, maxRank=1000){
  rank_value[rank_value==0] <- maxRank+1 #in case we are using sparse matrices
  insig <- rank_value > maxRank
  if(all(insig)) {
    return(0L)
  } else {
    rank_value[insig] <- maxRank+1
    rank_sum = sum(rank_value)
    len_sig <- length(rank_value)
    u_value = rank_sum - (len_sig * (len_sig + 1))/2
    auc = u_value/(len_sig * maxRank)
    auc = 1 - auc
    return(auc)
  }
}

#Calculate AUC for a list of signatures, from a ranks matrix
u_stat_signature_list = function(sig_list, ranks_matrix, maxRank=1000) {
  
  u_matrix <- sapply(sig_list, function(sig) {
    as.numeric(ranks_matrix[sig, lapply(.SD, function(x) u_stat(x,maxRank = maxRank)),.SDcols=-1, on="rn"])
  })
  
  rownames(u_matrix) <- colnames(ranks_matrix)[-1]
  return (u_matrix)
}


# Calculate feature ranks from expression data matrices
data_to_ranks_data_table = function(data) {
  dt <- as.data.table(data)
  rnaDT.ranks.dt <- dt[, lapply(.SD, function(x) frankv(x,ties.method="random",order=c(-1L)))]
  rnaDT.ranks.rownames <- rownames(data)
  rnaDT.ranks.dt.rn <- cbind(rn=rnaDT.ranks.rownames, rnaDT.ranks.dt)
  setkey(rnaDT.ranks.dt.rn, rn, physical = F)
  return(rnaDT.ranks.dt.rn)
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