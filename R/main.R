#' Filter single-cell data by cell type
#'
#' Apply CTfilter for a specific cell type to a query dataset. This function expands \code{\link{CTfilter}} to allow multi-level signatures to be used
#' for filtering; e.g. Tcells at level 1, and CD8.Tcells at level 2. To prepare markers and signatures in the correct format, store them in a list where
#' the list index corresponds to the signature level. Then select the desired cell at each level with the celltype parameter
#' (e.g. \code{celltype=c("T.cell","Tcell.CD8")}). Note that for the parameters listed below, you can also specify different values for different levels, by
#' setting them as vectors (e.g. sd.out=c(7,4)).
#'
#' @param query Seurat object containing a query data set - filtering will be applied to this object
#' @param celltype Cell type to preserve from the query data set (should be one of cell types in \code{names(markers)}). For multi-level filtering,
#'     provide nested subtypes as a vector (e.g. \code{celltype=c("T.cell","Tcell.CD8")})
#' @param CT.thresholds A named list with thresholds for each cell type - see function \code{\link{calculate_thresholds_CTfilter}}
#' @param markers List of markers for each cell type, for example \code{CTfilter::MCA.markers.Mm}
#' @param max.impurity Maximum fraction of impurity allowed in clusters to flag cells as "pure"
#' @param sd.in Maximum standard deviations from mean (Z-score) to identify outliers for selected signature
#' @param sd.out Maximum standard deviations from mean (Z-score) to identify outliers for all other signatures
#' @param min.gene.frac Only consider signatures covered by this fraction of genes in query set
#' @param assay Seurat assay to use
#' @param seed Random seed for cluster analysis
#' @param quiet Suppress all output
#' @param ... Additional parameters for \code{\link{CTfilter}}
#' @return A new metadata column \code{is.pure} is added to the query Seurat object, indicating which cells correspond to the desidered cell type.
#'     The \code{active.ident} is also set to this variable.
#' @examples
#' query <- CTfilter.multilevel(query, celltype=c("T.cell","Tcell.CD4"))
#' DimPlot(query)
#' query <- subset(query, subset=is.pure=="Pure")
#' @seealso \code{\link{calculate_thresholds_CTfilter()}} to calculate celltype-specific thresholds
#' @export
CTfilter.multilevel <- function(query, celltype="T.cell", CT.thresholds=NULL, markers=NULL, max.impurity=0.5, 
                     assay="RNA", min.gene.frac=0.5, sd.in=3, sd.out=7, seed=1234, quiet=FALSE, ...) {

  set.seed(seed)
  def.assay <- DefaultAssay(query)
  DefaultAssay(query) <- assay
  
  n.levels <- length(celltype)
  
  #Single-level run
  if (n.levels==1) {
     query <- CTfilter(query=query, celltype=celltype, CT.thresholds=CT.thresholds, markers=markers, max.impurity=max.impurity,
                       min.gene.frac=min.gene.frac, sd.in=sd.in, sd.out=sd.out, assay = assay, seed=seed, quiet=quiet, ...)
     return(query)
  }
  
  #Multi-level run
  if (!is.list(CT.thresholds) || length(CT.thresholds) < length(celltype)) {
     stop(sprintf("Please provide a list of thresholds (CT.thresholds) for each desidered celltype level (%s levels)", n.levels))
  }
  
  parameters <- list("max.impurity"=max.impurity,
                     "min.gene.frac"=min.gene.frac,
                     "sd.in"=sd.in,
                     "sd.out"=sd.out)
  parameters <- lapply(parameters, vectorize.parameters, lgt=n.levels)
  
  query$is.pure <- "Impure"
  sub <- query
  for (lev in 1:n.levels) {
     message(sprintf("--- Running CTfilter for level %i: %s celltype...",lev, celltype[lev]))
     
     sub <- CTfilter(query=sub, celltype=celltype[lev],
                      CT.thresholds=CT.thresholds[[lev]],
                      markers=markers[[lev]],
                      max.impurity=parameters$max.impurity[lev],
                      min.gene.frac=parameters$min.gene.frac[lev],
                      sd.in=parameters$sd.in[lev],
                      sd.out=parameters$sd.out[lev],
                      assay = assay, seed=seed, quiet=quiet, ...)
    
     sub <- subset(sub, subset=is.pure=="Pure")
  }
  pure.cells <- colnames(sub)
  query@meta.data[pure.cells, "is.pure"] <- "Pure"
  Idents(query) <- query$is.pure
  
  return(query)

}  
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
#' @param sd.in Maximum standard deviations from mean (Z-score) to identify outliers for selected signature
#' @param sd.out Maximum standard deviations from mean (Z-score) to identify outliers for all other signatures
#' @param ndim Number of dimensions for cluster analysis
#' @param resol Resolution for cluster analysis
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
#'     The \code{active.ident} is also set to this variable. Additionally, Z-scores for all signatures in \code{markers} are added to the metadata of the Seurat object.
#' @examples
#' query <- CTfilter(query, celltype="T.cell")
#' DimPlot(query)
#' query <- subset(query, subset=is.pure=="Pure")
#' @seealso \code{\link{calculate_thresholds_CTfilter()}} to calculate celltype-specific thresholds
#' @export
CTfilter <- function(query, celltype="T.cell", CT.thresholds=NULL, markers=NULL, max.impurity=0.5, 
                     ndim=30, resol=3, assay="RNA", genes.blacklist="Tcell.blacklist", min.gene.frac=0.5, 
                     sd.in=3, sd.out=7, rm.existing=TRUE,
                     max.iterations=10, stop.iterations=0.01, min.cells=100,
                     seed=1234, skip.normalize=FALSE, verbose=FALSE, quiet=FALSE) {
  
  set.seed(seed)
  def.assay <- DefaultAssay(query)
  DefaultAssay(query) <- assay
  celltype_CT <- paste0(celltype, "_CTfilter")
  
  #First guess the species from the gene names of the query
  mm.genes <- unique(unlist(MCA.markers.Mm))
  hs.genes <- unique(unlist(HCA.markers.Hs))
  
  mm.intersect <- length(intersect(mm.genes, rownames(query)))/length(mm.genes)
  hs.intersect <- length(intersect(hs.genes, rownames(query)))/length(hs.genes)
  if (max(mm.intersect, hs.intersect)<0.2) {
    warning("More than 80% of genes not found in reference signatures...did you remove genes from the query data?")
  }
  if (mm.intersect>hs.intersect) {
    species <- "Mouse"
  } else {
    species <- "Human"
  }
  
  if (is.null(markers)) { #Default
    if (species=="Human") {
      markers <- HCA.markers.Hs
    } else {
      markers <-  MCA.markers.Mm
    }
  }
  
  if (is.null(CT.thresholds)) { #Default
    if (species=="Human") {
      CT.thresholds <- Tcell.Hs.thr
    } else {
      CT.thresholds <- Tcell.TIL.thr 
    }
  }  

  if (!is.null(genes.blacklist)) {
    if (length(genes.blacklist)==1 && genes.blacklist == "Tcell.blacklist") {  #Default
      if (species=="Human") {
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
  
  #Allow setting different impurity levels in successive iterations
  max.impurity.vec <- vector(mode="numeric", max.iterations)
  for (i in 1:length(max.impurity.vec)) {
    if (i <= length(max.impurity)) {
      max.impurity.vec[i] <- max.impurity[i]
    } else {
      max.impurity.vec[i] <- max.impurity.vec[i-1]
    }
  }
  
  #Check that markers and threshold names correspond
  marker.cols <- paste0(names(markers), "_CTfilter")
  names(marker.cols) <- names(markers)
  marker.cols.pass <- marker.cols[marker.cols %in% rownames(CT.thresholds)]
  if (length(marker.cols.pass) ==0) {
    mess <- "Error. Could not match markers with threshold file"
    stop(mess)
  }
  markers.list.pass <- markers[names(marker.cols.pass)]
  markers.list.pass <- check_CTmarkers(obj=query, markers.list=markers.list.pass, min.gene.frac=min.gene.frac, verbose=verbose)
  
  if (!celltype %in% names(markers)) {
    mess <- sprintf("Cell type provided (%s) not found in marker list", celltype)
    stop(mess)
  }  
  if (!celltype_CT %in% rownames(CT.thresholds)) {
    mess <- sprintf("Cell type provided (%s) not found in thresholds file", celltype)
    stop(mess)
  } 
  if (!celltype %in% names(markers.list.pass)) {
    mess <- sprintf("Fewer than %.1f%% of genes from selected signature %s in query data", 100*min.gene.frac, celltype)
    stop(mess)
  }
  
  #Get Zscores
  query <- get_CTscores(obj=query, markers.list=markers.list.pass, rm.existing=rm.existing, bg=CT.thresholds, raw.score=F)
  
  sign.names <- names(markers.list.pass)
  
  meta <- query@meta.data
  filterCells <- c()
  for (sig in sign.names){
    sig.meta <- paste0(sig,"_CTfilter")
    if( sig.meta == celltype_CT ) {
      filterCells <- c(filterCells, which(meta[,sig.meta] < -sd.in))  # Z.score threshold for desired cell type
    } else {
      filterCells <- c(filterCells, which(meta[,sig.meta] > sd.out))   # Z.score threshold for contaminants
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
    mess <- sprintf("---- Iter %i - max.impurity=%.3f\n-- Detected %i non-pure cells for signature %s (%.2f%% of remaining cells)",
                     iter, imp.thr, n_rem, celltype, 100*frac.to.rem)
    if (verbose) {
        message(mess)
    }
  
    q$is.pure <- ifelse(q$clusterCT %in% filterCluster,"Impure","Pure")
    labs[colnames(q)] <- q$is.pure
    q <- subset(q, subset=is.pure=="Pure")
    
    if (frac.to.rem < stop.iterations | iter>=max.iterations | ncol(q)<min.cells) {
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


#' Calculate thresholds for CTfilter
#'
#' Given a reference set, calculate thresholds that identify outliers from the reference
#'
#' @param ref Seurat object containing the reference data set
#' @param markers List of markers for each cell type, for example \code{CTfilter::MCA.markers.Mm}
#' @param quant Quantile cutoff for score distribution
#' @param assay Seurat assay to use
#' @param min.gene.frac Only consider signatures covered by this fraction of genes in query set
#' @param min.sd Minimum value for standard deviation - set to this value if calculated standard deviation is lower 
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
calculate_thresholds_CTfilter <- function(ref, markers=NULL, quant=0.995, assay="RNA", min.gene.frac=0.5,
                                          min.sd=0.05, level=1, rm.existing=TRUE, chunk.size=1000, ncores=1, verbose=TRUE) {
  
  def.assay <- DefaultAssay(ref) 
  DefaultAssay(ref) <- assay
  if (is.null(markers)) {
     markers <- MCA.markers.Mm   #Default
  } 
  markers.list.pass <- check_CTmarkers(obj=ref, markers.list=markers, min.gene.frac=min.gene.frac, verbose=verbose)
  
  ref <- get_CTscores(obj=ref, markers.list=markers.list.pass, rm.existing=rm.existing, raw.score=TRUE,
                      chunk.size=chunk.size, ncores=ncores)
  
  
  sign.names <- paste0(names(markers.list.pass),"_CTfilter")
  
  ref_thr <- matrix(nrow=length(sign.names), ncol=2)
  rownames(ref_thr) <- sign.names
  colnames(ref_thr) <- c("mean","sd")
  
  for(sig in sign.names){
    #bulk <- ref@meta.data[,sig]
    bulk <- as.numeric(ref@meta.data[,sig])
    
    bulk <- bulk[bulk < quantile(bulk,p=quant)]
    ref_thr[sig,1] <- mean(bulk)
    ref_thr[sig,2] <- ifelse(sd(bulk) > min.sd, sd(bulk), min.sd)
 }

  if (!is.list(ref@misc$CTfilter)) {
     ref@misc$CTfilter <- list()
  } 
  ref@misc$CTfilter[[level]] <- ref_thr
  
  if (!is.list(ref@misc$CTfilter.markers)) {
    ref@misc$CTfilter.markers <- list()
  } 
  ref@misc$CTfilter.markers[[level]] <- markers.list.pass
  
  DefaultAssay(ref) <- def.assay
  message("Cell type thresholds available in ref@misc$CTfilter")
  message("Marker list available in ref@misc$CTfilter.markers")
  return(ref)
}

#' Calculate thresholds for CTfilter
#'
#' Calculate some stats after a CTfilter run.
#'
#' @param query A query object in Seurat format returned by CTfilter
#' @param celltype The desidered cell type used in CTfilter
#' @param sd.in Maximum standard deviations from mean (Z-score) to identify outliers for selected signature
#' @param sd.out Maximum standard deviations from mean (Z-score) to identify outliers for all other signatures
#' @param min.cells Minimum number of cells to report a specific cell type. Everything else will be grouped in 'Others'
#' @return Returns the query object with an additional column, reporting the cell type with the highest Z-score compared to the reference
#'     thresholds. In the slot \code{query@@misc$CTfilter.counts}, a dataframe reports the number of cells for each signature that exceed the Z-score threshold. Note that individual
#'     cells can be outliers for several signatures, so the total count does not correspond to the number of removed cells.
#' @examples
#' query <- CTfilter(query, celltype="Tcell")
#' query <- CTfilter.stats(query, celltype="Tcell")
#' Dimplot(query, group.by="top.zscore")
#' head(query@@misc$CTfilter.counts)
#' @seealso \code{\link{CTfilter()}} to apply signatures on a query dataset and filter on a specific cell type
#' @export


CTfilter.stats <- function(query, celltype="T.cell", sd.in=3, sd.out=7, min.cells=50) {
    
    celltype_CT <- paste0(celltype,"_CTfilter")
    query$top.zscore <- celltype_CT
    
    imp <- subset(query, is.pure=="Impure")
    
    cols <- grep("_CTfilter$", colnames(imp@meta.data),  perl=T, value=T)
    meta <- imp@meta.data[,cols]
    this.i <- which(colnames(meta) == celltype_CT)
    if (length(this.i) == 0) {
      stop(paste0("Celltype %s not found in query object",celltype))
    }  
    
    #Get raw counts of outlier from background distribution
    counts <- rep(0, length(cols))
    names(counts) <- cols
    for (type in cols) {
       if (type==celltype_CT) {
          counts[type] <- sum(meta[,type] < -sd.in)
       } else {
          counts[type] <- sum(meta[,type] > sd.out)
       }
    }
    counts <- sort(counts[counts>0], decreasing = T)
    counts <- as.data.frame(counts)
    
    ind <- which(rownames(counts) == celltype_CT)
    if (length(ind) > 0) {
       rownames(counts)[ind] <- paste0("NON_",celltype_CT)
    }
    rownames(counts) = substr(rownames(counts), 1, nchar(rownames(counts))-9)
    
    
    top.imp <- vector(length=nrow(meta))
    names(top.imp) <- rownames(meta)
    #Now assign cell type with top z-score
    meta.out <- meta[,-this.i]
    max.sd <- apply(meta.out, 1, max)
    max.which <- apply(meta.out, 1, which.max)
    max.type <- colnames(meta.out)[max.which]
    
    meta.in <- meta[,this.i]
    
    for (i in 1:length(top.imp)) {
       if (max.sd[i] > sd.out) {
          top.imp[i] <- max.type[i]
       } else if (meta.in[i] < -sd.in) {
          top.imp[i] <- paste0("NON_",celltype_CT)
       } else {
          top.imp[i] <- celltype_CT
       }
    }
    query$top.zscore[names(top.imp)] <- top.imp
    
    query$top.zscore = substr(query$top.zscore, 1, nchar(query$top.zscore)-9)
    
    #Group cell types with few cells
    tab <- table(query$top.zscore)
    join <- names(tab[tab<min.cells])
    query$top.zscore[query$top.zscore %in% join] <- "Others"
    
    query@misc$CTfilter.counts <- counts
    return(query)
}
  
  
  
