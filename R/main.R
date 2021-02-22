#' Filter single-cell data by cell type
#'
#' Apply CTfilter for a specific cell type to a query dataset
#'
#' @param query Seurat object containing a query data set - filtering will be applied to this object
#' @param celltype Cell type to preserve from the query data set (should be one of cell types in \code{names(markers)})
#' @param CT.thresholds A named list with thresholds for each cell type - see function \code{\link{calculate_thresholds_CTfilter}}
#' @param markers List of markers for each cell type, for example \code{CTfilter::MCA.markers.Mm}
#' @param max.impurity Maximum fraction of impurity allowed in clusters to flag cells as "pure". Can be either:
#' \itemize{
#'   \item{Single number between 0 and 1 - the same impurity threshold is applied to all iterations}
#'   \item{Vector of numbers between 0 and 1 - specific impurity thresholds for successive iterations e.g. \code{max.impurity=c(0.7,0.5,0.3)}}
#' }
#' @param ndim Number of dimensions for cluster analysis
#' @param resul Resolution for cluster analysis
#' @param assay Seurat assay to use
#' @param min.gene.frac Only consider signatures covered by this fraction of genes in query set
#' @param max.iterations Maximum number of iterations
#' @param stop.iterations Stop iterating if fewer than this fraction of cells were removed in the last iteration
#' @param min.cells Stop iterating if fewer than this number of cells is left
#' @param rm.existing Overwrite existing CTfilter scores in query object
#' @param genes.blacklist Genes blacklisted from variable features. The default loads the list of genes in \code{CTfilter::genes.blacklist.Mm};
#'     you may deactivate blacklisting by setting \code{genes.blacklist=NULL}
#' @param skip.normalize Skip data normalization
#' @param seed Random seed for cluster analysis
#' @param verbose Verbose output
#' @param quiet Suppress all output
#' @return A new metadata column \code{is.pure} is added to the query Seurat object, indicating which cells correspond to the desidered cell type.
#'     The \code{active.ident} is also set to this variable. Additionally, scores for all signatures in \code{markers} are added to the metadata of the Seurat object.
#' @examples
#' query <- CTfilter(query, celltype="T.cell")
#' DimPlot(query)
#' query <- subset(query, subset=is.pure=="Pure")
#' @seealso \code{\link{calculate_thresholds_CTfilter()}} to calculate celltype-specific thresholds
#' @export
CTfilter <- function(query, celltype="T.cell", CT.thresholds=NULL, markers=NULL, max.impurity=0.5, 
                     ndim=30, resol=3, assay="RNA", genes.blacklist="Tcell.blacklist", min.gene.frac=0.5, rm.existing=TRUE,
                     max.iterations=10, stop.iterations=0.01, min.cells=100,
                     seed=1234, skip.normalize=FALSE, verbose=FALSE, quiet=FALSE) {
  
  set.seed(seed)
  def.assay <- DefaultAssay(query)
  DefaultAssay(query) <- assay
  
  if (is.null(CT.thresholds)) {
    CT.thresholds <- Tcell.TIL.thr #Default
  }  
  if (is.null(markers)) {
    markers <- MCA.markers.Mm   #Default
  }
  if (!is.null(genes.blacklist)) {
    if (genes.blacklist == "Tcell.blacklist") {
       genes.blacklist <- genes.blacklist.Mm  #Default
    }  
    if (is.list(genes.blacklist)) {
      genes.blacklist <- unlist(genes.blacklist)
    }
    genes.blacklist <- unique(genes.blacklist)
  }
  #Allow setting different impurity levels in successive iterations
  max.impurity.vec <- vector(mode="numeric", max.iterations)
  for (i in 1:length(max.impurity.vec)) {
    if (i <= length(max.impurity)) {
      max.impurity.vec[i] <- max.impurity[i]
    } else {
      max.impurity.vec[i] <- max.impurity.vec[i-1]
    }
  }
  
  markers.list.pass <- check_CTmarkers(obj=query, markers.list=markers, min.gene.frac=min.gene.frac, verbose=verbose)
  
  if (!celltype %in% names(markers)) {
    mess <- sprintf("Cell type provided (%s) not found in marker list", celltype)
    stop(mess)
  }  
  if (!celltype %in% names(markers.list.pass)) {
    mess <- sprintf("Fewer than %.1f%% of genes from selected signature %s in query data", 100*min.gene.frac, celltype)
    stop(mess)
  }
  
  query <- get_CTscores(obj=query, markers.list=markers.list.pass, rm.existing=rm.existing)
  
  sign.names <- names(markers.list.pass)
  
  meta <- query@meta.data
  filterCells <- c()
  for (sig in sign.names){
    sig.meta <- paste0(sig,"_CTfilter")
    if( sig == celltype ) {
      filterCells <- c(filterCells, which(meta[,sig.meta] < -1*CT.thresholds[sig.meta]))
    } else {
      filterCells <- c(filterCells, which(meta[,sig.meta] > CT.thresholds[sig.meta]))
    }
  }
  filterCells <- unique(filterCells)
  filterCells.ID <- colnames(query)[filterCells]
  
  #The vector labs will hold the global cells status to return
  labs <- rep("Pure", dim(query)[2])
  names(labs) <- colnames(query)
  q <- query
  
  #Start iterations
  for (iter in 1:max.iterations) {
    
    filterCells.this <- filterCells.ID[filterCells.ID %in% colnames(q)]
    imp.thr <- max.impurity.vec[iter]
    
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
    frac.to.rem <- n_rem/ncol(q)
    mess <- sprintf("---- Iter %i - max.impurity=%.3f\n-- Detected %i non-pure cells for signature %s (%.2f%% of remaining pure cells)",
                     iter, imp.thr, n_rem, celltype, 100*frac.to.rem)
    if (verbose) {
        message(mess)
    }
  
    q$is.pure <- ifelse(q$clusterCT %in% filterCluster,"Impure","Pure")
    labs[colnames(q)] <- q$is.pure
    q <- subset(q, subset=is.pure=="Pure")
    
    if (frac.to.rem <= stop.iterations | iter>=max.iterations | ncol(q)<min.cells) {
      #Return clusters and active idents for easy filtering
      n_rem <- sum(labs=="Impure")
      frac.to.rem <- n_rem/length(labs)
      mess <- sprintf("CTfilter: Detected %i non-pure cells for signature %s - %.2f%% cells marked for removal (active.ident)",
                         n_rem, celltype, 100*frac.to.rem)
      if (!quiet) {
        message(mess)
      }
      query$is.pure <- labs
      Idents(query) <- query$is.pure
      DefaultAssay(query) <- def.assay
      return(query)
    }
  }  
  return()
}

#' Filter single-cell data by cell type
#'
#' Apply CTfilter for a specific cell type to a query dataset
#'
#' @param query Seurat object containing a query data set - filtering will be applied to this object
#' @param celltype Cell type to preserve from the query data set (should be one of cell types in \code{names(markers)})
#' @param CT.thresholds A named list with thresholds for each cell type - see function \code{\link{calculate_thresholds_CTfilter}}
#' @param markers List of markers for each cell type, for example \code{CTfilter::MCA.markers.Mm}
#' @param max.impurity Maximum fraction of impurity allowed in clusters to flag cells as "pure"
#' @param ndim Number of dimensions for cluster analysis
#' @param resul Resolution for cluster analysis
#' @param assay Seurat assay to use
#' @param min.gene.frac Only consider signatures covered by this fraction of genes in query set
#' @param rm.existing Overwrite existing CTfilter scores in query object
#' @param genes.blacklist Genes blacklisted from variable features
#' @param skip.normalize Skip data normalization
#' @param seed Random seed for cluster analysis
#' @param verbose Verbose output
#' @param quiet Suppress all output
#' @return A new metadata column \code{is.pure} is added to the query Seurat object, indicating which cells correspond to the desidered cell type.
#'     The \code{active.ident} is also set to this variable. Additionally, scores for all signatures in \code{markers} are added to the metadata of the Seurat object.
#' @examples
#' query <- CTfilter(query, celltype="T.cell")
#' DimPlot(query)
#' query <- subset(query, subset=is.pure=="Pure")
#' @seealso \code{\link{calculate_thresholds_CTfilter()}} to calculate celltype-specific thresholds
#' @export
CTfilter.single <- function(query, celltype="T.cell", CT.thresholds=NULL, markers=NULL, max.impurity=0.5, 
                     ndim=30, resol=3, assay="RNA", genes.blacklist=NULL, min.gene.frac=0.5, rm.existing=TRUE,
                     seed=1234, skip.normalize=FALSE, verbose=FALSE, quiet=FALSE) {
  
  set.seed(seed)
  def.assay <- DefaultAssay(query)
  DefaultAssay(query) <- assay
  
  if (is.null(CT.thresholds)) {
    CT.thresholds <- Tcell.TIL.thr #Default
  }  
  if (is.null(markers)) {
    markers <- MCA.markers.Mm   #Default
  } 
  if (is.null(genes.blacklist)) {
    genes.blacklist <- genes.blacklist.Mm  #Default
  }  
  if (is.list(genes.blacklist)) {
    genes.blacklist <- unlist(genes.blacklist)
  }
  genes.blacklist <- unique(genes.blacklist)
  
  markers.list.pass <- check_CTmarkers(obj=query, markers.list=markers, min.gene.frac=min.gene.frac, verbose=verbose)
  
  if (!celltype %in% names(markers)) {
    mess <- sprintf("Cell type provided (%s) not found in marker list", celltype)
    stop(mess)
  }  
  if (!celltype %in% names(markers.list.pass)) {
    mess <- sprintf("Fewer than %.1f%% of genes from selected signature %s in query data", 100*min.gene.frac, celltype)
    stop(mess)
  }
  
  query <- get_CTscores(obj=query, markers.list=markers.list.pass, rm.existing=rm.existing)
  
  sign.names <- names(markers.list.pass)
  
  meta <- query@meta.data
  filterCells <- c()
  for (sig in sign.names){
    sig.meta <- paste0(sig,"_CTfilter")
    if( sig == celltype ) {
      filterCells <- c(filterCells, which(meta[,sig.meta] < -1*CT.thresholds[sig.meta]))
    } else {
      filterCells <- c(filterCells, which(meta[,sig.meta] > CT.thresholds[sig.meta]))
    }
  }
  filterCells <- unique(filterCells)
  filterCellsId <- colnames(query)[filterCells]
  
  ##High resolution clustering to detect dense regions of undesired cell types
  q <- query
  if (!skip.normalize) {
    q <- NormalizeData(q, verbose = FALSE)
  }
  q <- FindVariableFeatures(q, selection.method = "vst", nfeatures = 500, verbose = FALSE)
  q@assays[[assay]]@var.features <- setdiff(q@assays[[assay]]@var.features, genes.blacklist)
  q <- ScaleData(q, verbose=FALSE)
  q <- RunPCA(q, features = q@assays[[assay]]@var.features, verbose = FALSE)
  q <- FindNeighbors(q, reduction = "pca", dims = 1:ndim, k.param = 5, verbose=FALSE)
  q  <- FindClusters(q, resolution = resol, verbose = FALSE)
  q$clusterCT <- q@active.ident
  
  m <- q@meta.data
  impure.freq <- tapply(row.names(m), m$clusterCT,function(x) {mean(x %in% row.names(m)[filterCells])})
  
  filterCluster <- names(impure.freq)[impure.freq > max.impurity]
  n_rem <- sum(q$clusterCT %in% filterCluster)
  mess <- sprintf("Detected %i non-pure cells for signature %s - %.2f%% cells marked for removal (active.ident)",
                     n_rem, celltype, n_rem*100/ncol(q))
  if (!quiet) {
    message(mess)
  }
  
  #Return clusters and active idents for easy filtering
  query$is.pure <- ifelse(q$clusterCT %in% filterCluster,"Impure","Pure")
  Idents(query) <- query$is.pure
  DefaultAssay(query) <- def.assay
  return(query)
}


#' Calculate thresholds for CTfilter
#'
#' Given a reference set, calculate thresholds that identify outliers from the reference
#'
#' @param ref Seurat object containing the reference data set
#' @param markers List of markers for each cell type, for example \code{CTfilter::MCA.markers.Mm}
#' @param sd.dev Maximum standard deviations from mean to identify outliers
#' @param quant Quantile cutoff for score distribution
#' @param assay Seurat assay to use
#' @param min.gene.frac Only consider signatures covered by this fraction of genes in query set
#' @param level Annotation level - thresholds are saved in list element \code{ref@@misc$CTfilter[[level]]}
#' @param rm.existing Overwrite existing CTfilter scores in query object
#' @param verbose Verbose output
#' @return Return the \code{ref} reference object, with celltype-specific thresholds in the field \code{ref@@misc$CTfilter}. Scores for individual signatures are
#'     returned as metadata in the Seurat object
#' @examples
#' library(ProjecTILs)
#' ref <- load.reference.map()
#' ref <- calculate_thresholds_CTfilter(ref)
#' ref@@misc$CTfilter
#' @seealso \code{\link{CTfilter()}} to apply signatures on a query dataset and filter on a specific cell type
#' @export
calculate_thresholds_CTfilter <- function(ref, markers=NULL, sd.dev=7, quant=0.995, assay="RNA", min.gene.frac=0.5,
                                          level=1, rm.existing=TRUE, verbose=TRUE) {
  
  def.assay <- DefaultAssay(ref) 
  DefaultAssay(ref) <- assay
  if (is.null(markers)) {
     markers <- MCA.markers.Mm   #Default
  } 
  markers.list.pass <- check_CTmarkers(obj=ref, markers.list=markers, min.gene.frac=min.gene.frac, verbose=verbose)
  ref <- get_CTscores(obj=ref, markers.list=markers.list.pass, rm.existing=rm.existing)
  
  sign.names <- paste0(names(markers.list.pass),"_CTfilter")
  ref_thr <- vector(length=length(sign.names))
  names(ref_thr) <- sign.names
  for(sig in sign.names){
    bulk <- ref@meta.data[,sig]
    bulk <- bulk[bulk < quantile(bulk,p=quant)]
    ref_thr[sig] <- mean(bulk) + sd.dev*sd(bulk)
 }

  if (!is.list(ref@misc$CTfilter)) {
     ref@misc$CTfilter <- list()
  } 
  ref@misc$CTfilter[[level]] <- ref_thr
  
  DefaultAssay(ref) <- def.assay
  message("Cell type thresholds available in ref@misc$CTfilter")
  return(ref)
}