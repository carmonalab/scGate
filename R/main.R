#' Filter single-cell data by cell type
#'
#' Apply CTfilter for a specific cell type to a query dataset
#'
#' @param query Seurat object containing a query data set - filtering will be applied to this object
#' @param celltype Cell type to preserve from the query data set (should be one of cell types in `names(markers)`)
#' @param CT.thresholds A named list with thresholds for each cell type - see function `calculate_thresholds_CTfilter`
#' @param markers List of markers for each cell type, for example `CTfilter::MCA.markers.Mm`
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
#' @return A new metadata column `is.pure` is added to the query Seurat object, indicating which cells correspond to the desidered cell type.
#'     The `active.ident` is also set to this variable. Additionally, scores for all signatures in `markers` are added to the metadata of the Seurat object.
#' @examples
#' query <- CTfilter(query, celltype="T.cell")
#' DimPlot(query)
#' query <- subset(query, subset=is.pure=="Pure")
#' @seealso [calculate_thresholds_CTfilter()] to calculate celltype-specific thresholds
#' @export
CTfilter <- function(query, celltype="T.cell", CT.thresholds=NULL, markers=NULL, max.impurity=0.5,
                     ndim=5, resol=3, assay="RNA", genes.blacklist=NULL, min.gene.frac=0.5, rm.existing=TRUE,
                     seed=1234, skip.normalize=FALSE, verbose=FALSE) {
  
  set.seed(seed)
  def.assay <- DefaultAssay(query)
  DefaultAssay(query) <- assay
  
  if (is.null(CT.thresholds)) {
    CT.thresholds <- Tcell.TIL.thr #Default
  }  
  if (is.null(markers)) {
    markers.list <- MCA.markers.Mm   #Default
  } 
  if (is.null(genes.blacklist)) {
    genes.blacklist <- genes.blacklist.Mm  #Default
  }  
  if (is.list(genes.blacklist)) {
    genes.blacklist <- unlist(genes.blacklist)
  }
  genes.blacklist <- unique(genes.blacklist)
  
  markers.list.pass <- check_CTmarkers(obj=query, markers.list=markers.list, min.gene.frac=min.gene.frac, verbose=verbose)
  
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
  message <- sprintf("Detected %i non-pure cells for signature %s - %.2f%% cells marked for removal (active.ident)",
              n_rem, celltype, n_rem*100/ncol(q))
  message(message)
  
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
#' @param markers List of markers for each cell type, for example `CTfilter::MCA.markers.Mm`
#' @param sd.dev Maximum standard deviations from mean to identify outliers
#' @param quant Quantile cutoff for score distribution
#' @param assay Seurat assay to use
#' @param min.gene.frac Only consider signatures covered by this fraction of genes in query set
#' @param rm.existing Overwrite existing CTfilter scores in query object
#' @param verbose Verbose output
#' @return A new metadata column `is.pure` is added to the query Seurat object, indicating which cells correspond to the desidered cell type.
#'     The `active.ident` is also set to this variable. Additionally, scores for all signatures in `markers` are added to the metadata of the Seurat object.
#' @examples
#' library(ProjecTILs)
#' ref <- load.reference.map()
#' ref <- calculate_thresholds_CTfilter(ref)
#' ref@@misc$CTfilter
#' @seealso [CTfilter()] to apply signatures on a query dataset and filter on a specific cell type
#' @export
calculate_thresholds_CTfilter <- function(ref, markers=NULL, sd.dev=7, quant=0.995, assay="RNA", min.gene.frac=0.5,
                                          rm.existing=TRUE, verbose=TRUE) {
  
  def.assay <- DefaultAssay(ref) 
  DefaultAssay(ref) <- assay
  if (is.null(markers)) {
     markers.list <- MCA.markers.Mm   #Default
  } 
  markers.list.pass <- check_CTmarkers(obj=ref, markers.list=markers.list, min.gene.frac=min.gene.frac, verbose=verbose)
  ref <- get_CTscores(obj=ref, markers.list=markers.list.pass, rm.existing=rm.existing)
  
  sign.names <- paste0(names(markers.list.pass),"_CTfilter")
  ref_thr <- vector(length=length(sign.names))
  names(ref_thr) <- sign.names
  for(sig in sign.names){
    bulk <- ref@meta.data[,sig]
    bulk <- bulk[bulk < quantile(bulk,p=quant)]
    ref_thr[sig] <- mean(bulk) + sd.dev*sd(bulk)
 }
  
  ref@misc$CTfilter <- ref_thr
  DefaultAssay(ref) <- def.assay
  message("Cell type thresholds available in ref@misc$CTfilter")
  return(ref)
}