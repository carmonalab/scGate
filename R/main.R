CTfilter <- function(query, celltype="T.cell", CT.thresholds=NULL, markers=NULL, genes.blacklist=NULL,
                     assay="RNA", min.gene.frac=0.5, rm.existing=TRUE, seed=1234,
                     ndim=5, resol=3, max.impurity=0.5, skip.normalize=FALSE, verbose=FALSE) {
  
  set.seed(seed)
  def.assay <- DefaultAssay(query)
  DefaultAssay(query) <- assay
  
  if (is.null(markers)) {
    markers.list <- readRDS("aux/MCA_Markers.rds")   #HERE change to internal package data
    names(markers.list) <- paste0(names(markers.list),"_CT")
  } else {
    markers.list <- readRDS(markers)
  }
  if (is.null(genes.blacklist)) {
    genes.blacklist <- "aux/genesBlacklist.Rds"   #HERE change to internal package data
  }  
  genes.bl <- readRDS(genes.blacklist)
  if (is.list(genes.bl)) {
    genes.bl <- unlist(genes.bl)
  }
  genes.bl <- unique(genes.bl)
  
  markers.list.pass <- check_CTmarkers(obj=query, markers.list=markers.list, min.gene.frac=min.gene.frac, verbose=verbose)
  celltype <- paste0(celltype,"_CT")
  
  if (!celltype %in% names(markers.list)) {
    message <- sprintf("Cell type provided (%s) not found in marker list", celltype)
    stop(message)
  }  
  if (!celltype %in% names(markers.list.pass)) {
     message <- sprintf("Fewer than %.1f%% of genes from selected signature %s in query data", 100*min.gene.frac, celltype)
     stop(message)
  }
  
  query <- get_CTscores(obj=query, markers.list=markers.list.pass, rm.existing=rm.existing)
  
  sign.names <- names(markers.list.pass)
  Celltype_score_th_ref <- ref@misc$CTFilter
  
  meta <- query@meta.data
  filterCells <- c()
  for (sig in sign.names){
    sig.meta <- paste0("CelltypeScore_",sig)
    if( sig == celltype ) {
       filterCells <- c(filterCells, which(meta[,sig.meta] < -1*CT.thresholds[sig.meta]))
    } else {
       filterCells <- c(filterCells, which(meta[,sig.meta] > CT.thresholds[sig.meta]))
    }
  }
  filterCells <- unique(filterCells)
  filterCellsId <- colnames(query)[filterCells]
  
  ##High resolution clustering to detect dense regions of undesired cell types
  if (!skip.normalize) {
     query <- NormalizeData(query, verbose = FALSE)
  }
  query <- FindVariableFeatures(query, selection.method = "vst", nfeatures = 500, verbose = FALSE)
  query@assays[[assay]]@var.features <- setdiff(query@assays[[assay]]@var.features, genes.bl)
  query <- ScaleData(query, verbose=FALSE)
  query <- RunPCA(query, features = query@assays[[assay]]@var.features, verbose = FALSE)
  query <- FindNeighbors(query, reduction = "pca", dims = 1:ndim, k.param = 5, verbose=FALSE)
  query  <- FindClusters(query, resolution = resol, verbose = FALSE)
  query$clusterCT <- query@active.ident
  
  m <- query@meta.data
  impure.freq <- tapply(row.names(m), m$clusterCT,function(x) {mean(x %in% row.names(m)[filterCells])})

  filterCluster <- names(impure.freq)[impure.freq > max.impurity]
  n_rem <- sum(query$clusterCT %in% filterCluster)
  message <- sprintf("Detected %i non-pure cells for signature %s - %.2f%% cells marked for removal (active.ident)",
              n_rem, celltype, n_rem*100/ncol(query))
  message(message)
  
  #Return clusters and active idents for easy filtering
  query$is.pure <- ifelse(query$clusterCT %in% filterCluster,"Impure","Pure")
  Idents(query) <- query$is.pure
  DefaultAssay(query) <- def.assay
  return(query)
}

calculate_thresholds_CTfilter <- function(ref, markers=NULL, sd.dev=7, quant=0.995, assay="RNA", min.gene.frac=0.5, rm.existing=TRUE) {
  
  def.assay <- DefaultAssay(ref) 
  DefaultAssay(ref) <- assay
  if (is.null(markers)) {
     markers.list <- readRDS("aux/MCA_Markers.rds")   #HERE change to internal package data
     names(markers.list) <- paste0(names(markers.list),"_CT")
  } else {
     markers.list <- readRDS(markers)
  }
  
  markers.list.pass <- check_CTmarkers(obj=ref, markers.list=markers.list, min.gene.frac=min.gene.frac, verbose=verbose)
  ref <- get_CTscores(obj=ref, markers.list=markers.list.pass, rm.existing=rm.existing)
  
  sign.names <- paste0("CelltypeScore_",names(markers.list.pass))
  ref_thr <- vector(length=length(sign.names))
  for(sig in sign.names){
    bulk <- ref@meta.data[,sig]
    bulk <- bulk[bulk < quantile(bulk,p=quant)]
    ref_thr[sig] <- mean(bulk) + sd.dev*sd(bulk)
 }
  
  ref@misc$CTfilter <- ref_thr
  DefaultAssay(ref) <- def.assay
  return(ref)
}

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

get_CTscores <- function(obj, markers.list, rm.existing=TRUE) {
  
  if (rm.existing) {  #remove existing CTfilter scores - necessary because AddModuleScore won't overwrite them
    index.rm <- grep("_CT$",colnames(obj@meta.data), perl=T)
    if (length(index.rm)>0) {
      obj@meta.data <- subset(obj@meta.data, select = -index.rm)
    }
  }
  
  obj <- suppressWarnings(AddModuleScore(obj, features = markers.list, name="CelltypeScore"))
  sign.names <- paste0("CelltypeScore_",names(markers.list))
  
  index.start <- ncol(obj@meta.data)-length(markers.list)+1
  index.end <- ncol(obj@meta.data)
  
  names(obj@meta.data)[index.start:index.end] <- sign.names
  return(obj)
}
