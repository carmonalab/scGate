#Helper function to subset signatures on genes observed in the query set
check_CTmarkers <- function(obj, markers.list, min.gene.frac=0.5, verbose=TRUE) {
  
  markers.list.exp <- unlist(lapply(markers.list,function(x){mean(x %in% rownames(obj))}))
  
  markers.list.pass <- markers.list[markers.list.exp >= min.gene.frac]
  markers.list.pass <- lapply(markers.list.pass,function(x){x[x %in% rownames(obj)]})
  
  skip <- names(markers.list.exp)[markers.list.exp < min.gene.frac]
  
  missing <- 100*length(skip)/length(markers.list)
  if (missing > 50) {
     warning(sprintf("Warning: %s %% of signatures could not be evaluated, too many genes are missing from your data"))
  } 
  
  if (verbose==TRUE & length(skip) > 0) {
    skip <- paste(skip, collapse=", ")
    message <- paste0("Not enough genes to evaluate the following signatures:\n",skip,"\n")
    message <- paste0(message, "Probably this just means that these cells are not present in your dataset")
    warning(message)
  }
  
  return(markers.list.pass)
}

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
                         method=c("UCell","ModuleScore"), chunk.size=1000, ncores=1) {
  
  method.use <- method[1]
  
  if (rm.existing) {  #remove existing scGate scores
    index.rm <- grep("_scGate|_Zscore",colnames(obj@meta.data), perl=T)
    if (length(index.rm)>0) {
      obj@meta.data <- subset(obj@meta.data, select = -index.rm)
    }
  }
  
  if (method.use == "UCell") {
     obj <- UCell::AddModuleScore_UCell(obj, features=markers.list, chunk.size=chunk.size, ncores=ncores, name = "_scGate")
  } else if (method.use == "ModuleScore") {
     obj <- suppressWarnings(AddModuleScore(obj, features = markers.list, name="scGateScore"))
  } else {
     stop("Please give a valid method for signature scoring (see 'method' parameter)")
  }

  if (!z.score) {
    return(obj) 
  }
  
  #Convert to Z-score
  cols <- paste0(names(markers.list),"_scGate")
  zcols <- paste0(names(markers.list),"_Zscore")
  
  for (i in seq_along(cols)) {
    sig <- cols[i]
    zsig <- zcols[i]
    obj@meta.data[,zsig] <- (obj@meta.data[,sig] - bg[sig,"mean"])/bg[sig,"sd"]
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