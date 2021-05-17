#' Filter single-cell data by cell type (multi-level)
#'
#' Apply scGate for a specific cell type to a query dataset. This function expands \code{\link{scGate}} to allow multi-level signatures to be used
#' for filtering; e.g. Tcells at level 1, and CD8.Tcells at level 2. To prepare markers and signatures in the correct format, store them in a list where
#' the list index corresponds to the signature level. Then select the desired cell at each level with the `celltype` parameter, separating nested subtypes 
#' using the underscore: e.g. \code{celltype="T.cell_Tcell.CD8"}.
#' Note that for the parameters listed below, you can also specify different values for different levels, by
#' setting them as vectors (e.g. sd.out=c(7,4)).
#'
#' @param data Seurat object containing a query data set - filtering will be applied to this object
#' @param gating.model The background expected mean and SD for each cell type - see function \code{\link{train_scGate}}
#'     provide nested subtypes as a vector (e.g. \code{celltype=c("T.cell","Tcell.CD8")})
#' @param gating.model The background expected mean and SD for each cell type - see function \code{\link{train_scGate}}
#' @param max.impurity Maximum fraction of impurity allowed in clusters to flag cells as "pure"
#' @param sd.in Maximum standard deviations from mean (Z-score) to identify outliers for selected signature
#' @param sd.out Maximum standard deviations from mean (Z-score) to identify outliers for all other signatures
#' @param assay Seurat assay to use
#' @param seed Integer seed for random number generator
#' @param quiet Suppress all output
#' @param ... Additional parameters for \code{\link{scGate}}
#' @return A new metadata column \code{is.pure} is added to the query Seurat object, indicating which cells correspond to the desidered cell type.
#'     The \code{active.ident} is also set to this variable.
#' @examples
#' query <- scGate.multilevel(query)
#' DimPlot(query)
#' query <- subset(query, subset=is.pure=="Pure")
#' @seealso \code{\link{train_scGate()}} to generate celltype-specific models
#' @export
scGate.multilevel <- function(data, gating.model=NULL, max.impurity=0.5, 
                     assay="RNA", sd.in=3, sd.out=7, seed=123, quiet=FALSE, ...) {

  
  set.seed(seed)
  def.assay <- DefaultAssay(data)
  DefaultAssay(data) <- assay
  
  species <- detect_species(data)
  
  if (is.null(gating.model)) { #Default
    if (species=="human") {
      gating.model <- scGate_human_TIL_model
    } else {
      gating.model <- scGate_mouse_TIL_model 
    }
  }

  if (!is.list(gating.model)) {
     gating.model[[1]] <- gating.model 
  }
  n.levels <- length(gating.model)
  
  #Single-level run
  if (n.levels==1) {
     data <- scGate(data, gating.model=gating.model[[1]],
                     max.impurity=max.impurity,
                     sd.in=sd.in, sd.out=sd.out,
                     assay = assay, seed=seed, quiet=quiet, ...)
     return(data)
  }
  
  #Multi-level run
  parameters <- list("max.impurity"=max.impurity,
                     "sd.in"=sd.in,
                     "sd.out"=sd.out)
  parameters <- lapply(parameters, vectorize.parameters, lgt=n.levels)
  
  data$is.pure <- "Impure"
  sub <- data
  for (lev in 1:n.levels) {
     message(sprintf("--- Running scGate for level %i: ",lev))
    
     sub <- scGate(data=sub,
                      gating.model=gating.model[[lev]],
                      max.impurity=parameters$max.impurity[[lev]],
                      sd.in=parameters$sd.in[lev],
                      sd.out=parameters$sd.out[lev],
                      assay = assay, 
                      seed=seed, quiet=quiet, ...)
     
     meta.cols <- grep("_scGate|_Zscore",colnames(sub@meta.data), perl=T, value = T)
     data <- AddMetaData(data, metadata = sub@meta.data[,meta.cols])  #NB: this will generated NAs for 2nd+ level signatures
     
     sub <- subset(sub, subset=is.pure=="Pure")
  }
  pure.cells <- colnames(sub)
  data@meta.data[pure.cells, "is.pure"] <- "Pure"

  Idents(data) <- data$is.pure
  
  return(data)
}  
#' Filter single-cell data by cell type
#'
#' Apply scGate to filter specific cell types in a query dataset
#'
#' @param data Seurat object containing a query data set - filtering will be applied to this object
#' @param gating.model The background expected mean and SD for each cell type - see function \code{\link{train_scGate}}
#' @param max.impurity Maximum fraction of impurity allowed in clusters to flag cells as "pure". Can be either:
#' \itemize{
#'   \item{Single number between 0 and 1 - the same impurity threshold is applied to all iterations}
#'   \item{Vector of numbers between 0 and 1 - specific impurity thresholds for successive iterations e.g. \code{max.impurity=c(0.7,0.5,0.3)}}
#' }
#' @param sd.in Maximum standard deviations from mean (Z-score) to identify outliers for selected signatures
#' @param sd.out Maximum standard deviations from mean (Z-score) to identify outliers for all other signatures
#' @param ndim Number of dimensions for cluster analysis
#' @param resol Resolution for cluster analysis
#' @param assay Seurat assay to use
#' @param chunk.size Number of cells per batch to be scored by the method
#' @param maxRank Maximum number of genes to rank per cell; above this rank, a given gene is considered as not expressed (used when 'method' is 'UCell' or 'AUCell')
#' @param ncores Number of processors for parallel processing (requires \code{future.apply})
#' @param max.iterations Maximum number of iterations
#' @param stop.iterations Stop iterating if fewer than this fraction of cells were removed in the last iteration
#' @param min.cells Stop iterating if fewer than this number of cells is left
#' @param additional.signatures A list of additional signatures, not included in the model, to be evaluated (e.g. a cycling signature). The scores for this
#'     list of signatures will be returned but not used for filtering.
#' @param genes.blacklist Genes blacklisted from variable features. The default loads the list of genes in \code{scGate::genes.blacklist.Mm};
#'     you may deactivate blacklisting by setting \code{genes.blacklist=NULL}
#' @param skip.normalize Skip data normalization
#' @param return_signature_scores Add signature scores and Z-scores to object metadata
#' @param seed Integer seed for random number generator
#' @param verbose Verbose output
#' @param quiet Suppress all output
#' @return A new metadata column \code{is.pure} is added to the query Seurat object, indicating which cells correspond to the desidered cell type(s).
#'     The \code{active.ident} is also set to this variable. Additionally, Z-scores for all signatures in \code{markers} are added to the metadata of the Seurat object.
#' @examples
#' query <- scGate(query)
#' DimPlot(query)
#' query <- subset(query, subset=is.pure=="Pure")
#' @seealso \code{\link{train_scGate()}} to generate celltype-specific models
#' @export

scGate <- function(data, gating.model=NULL, max.impurity=0.5, 
                     ndim=30, resol=3, assay="RNA", 
                     sd.in=3, sd.out=7,
                     chunk.size=1000, ncores=1, maxRank=1500,
                     max.iterations=10, stop.iterations=0.01, min.cells=100,
                     additional.signatures=NULL,
                     genes.blacklist="Tcell.blacklist", 
                     seed=123, skip.normalize=FALSE, 
                     return_signature_scores=TRUE, verbose=FALSE, quiet=FALSE) {
  
  set.seed(seed)
  if (ncores>1) {
     require(future.apply)
     future_param_seed <<- seed
     future_param_ncores <<- ncores
  }
  
  def.assay <- DefaultAssay(data)
  DefaultAssay(data) <- assay
  
  species <- detect_species(data)
  
  if (is.null(gating.model)) { #Default
    if (species=="human") {
      gating.model <- scGate_human_TIL_model
    } else {
      gating.model <- scGate_mouse_TIL_model 
    }
  }
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
  marker.names.pass <- intersect(names(markers), rownames(model))
  
  if (length(marker.names.pass) ==0) {
    mess <- "Error. Could not match markers with threshold file"
    stop(mess)
  }
  markers <- markers[marker.names.pass]
  
  #Check selected cell type(s), and autoexpand if required
  sign.names <- names(markers)
  celltype.pass <- check_selected_celltypes(pos.celltypes, db=sign.names, autocomplete=F)
  
  if (length(celltype.pass)<1) {
    mess <- sprintf("Could not find selected cell types in marker list")
    stop(mess)
  }  
  if (!quiet) {
    mess <- paste(celltype.pass, collapse=", ")
    message(sprintf("Running scGate for selected signatures: %s", mess))
  }
  
  #Get Zscores
  scores <- get_CTscores(obj=data, markers.list=markers, method=method, chunk.size=chunk.size, ncores=ncores,
                        bg=model, z.score=TRUE, maxRank=maxRank, additional.signatures=additional.signatures)
  
  #First filter on cells with minimum levels of selected signature(s)
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

  if (length(positive_select)<1) {
    stop("All cells were removed by scGate. Check filtering thresholds or input data")
  }
  
  #Second, filter on cells with eccessive levels on undesired signatures
  negative_select <- seq_along(rownames(scores))[-positive_select]

  for (sig in sign.names){
    sig.zscore <- paste0(sig,"_Zscore")
    
    if(! sig %in% celltype.pass) {
      negative_select <- c(negative_select, which(scores[,sig.zscore] > sd.out))   # Z.score threshold for contaminants
    }
  }
  negative_select <- unique(negative_select)
  negative_select.ID <- rownames(scores)[negative_select]
  
  #The vector labs will hold the global cells status to return
  labs <- rep("Pure", dim(data)[2])
  names(labs) <- Cells(data)
  tot.cells <- length(labs)
  q <- data
  
  #Start iterations
  for (iter in 1:max.iterations) {
    
    filterCells.this <- negative_select.ID[negative_select.ID %in% Cells(q)]
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
    frac.to.rem <- n_rem/tot.cells
    mess <- sprintf("---- Iter %i - max.impurity=%.3f\n-- Detected %i non-pure cells for selected signatures (%.2f%% of total cells)",
                     iter, imp.thr, n_rem, 100*frac.to.rem)
    if (verbose) {
        message(mess)
    }
  
    q$is.pure <- ifelse(q$clusterCT %in% filterCluster,"Impure","Pure")
    labs[colnames(q)] <- q$is.pure
    q <- subset(q, subset=is.pure=="Pure")
    
    if (frac.to.rem < stop.iterations | iter>=max.iterations | ncol(q)<min.cells) {
      #Return clusters and active idents for easy filtering
      n_rem <- sum(labs=="Impure")
      frac.to.rem <- n_rem/tot.cells
      mess <- sprintf("scGate: Detected %i non-pure cells for selected signatures - %.2f%% cells marked for removal (active.ident)",
                         n_rem, 100*frac.to.rem)
      if (!quiet) {
        message(mess)
      }
      if (return_signature_scores) {
         data <- AddMetaData(data, scores)
      }
      
      data$is.pure <- labs
      Idents(data) <- data$is.pure
      DefaultAssay(data) <- def.assay
      return(data)
    }
  }  
  return()
}


#' Generate a scGate model
#'
#' Given a reference set, calculate expected mean and SD for all cell types in the reference set. Deviations from this expected distribution (Z-score)
#' can then be used to gate for specific cell types in any query dataset (see function \code{scGate})
#'
#' @param ref Seurat object containing the reference data set
#' @param markers List of markers for each cell type, for example \code{scGate::MCA.markers.Mm}
#' @param positive_celltypes List of celltypes included in the reference set. These must be one or more from the list of markers (`\code{markers} parameter).
#'     Only the `positive_celltypes` will be gated when the model is applied to a new query dataset
#' @param autocomplete Automatically autocomplete cell types in \code{positive_celltypes} that start with same prefix
#'     (e.g. B.cell gates for B.cell.1, B.cell.2 and B.cell.Plasmocyte signatures) 
#' @param quant Quantile cutoff for score distribution
#' @param assay Seurat assay to use
#' @param method Scoring method for cell signatures (default \code{UCell})
#' @param min.sd Minimum value for standard deviation - if calculated standard deviation is lower than min.sd, it is set to min.sd
#' @param chunk.size Number of cells per batch to be scored by the method
#' @param maxRank Maximum number of genes to rank per cell; above this rank, a given gene is considered as not expressed (used when 'method' is 'UCell' or 'AUCell')
#' @param ncores Number of processors for parallel processing (requires \code{future.apply})
#' @param seed Integer seed for random number generator
#' @param verbose Verbose output
#' @param quiet Suppress all STDOUT output
#' @return Returns a \code{scGate_Model} object, containing the expected signature score distribution (background model) for all evaluated signatures.
#' @examples
#' library(ProjecTILs)
#' ref <- load.reference.map()
#' scGate.model <- train_scGate(ref, positive_celltypes='T.cell')
#' head(scGate.model@@bg_model)
#' @seealso \code{\link{scGate()}} to apply signatures on a query dataset and filter on specific cell type(s)
#' @export
train_scGate <- function(ref, markers=NULL, positive_celltypes=NULL, autocomplete=T,
                                          maxRank=1500, min.sd=0.02,
                                          method=c("UCell","AUCell","ModuleScore"), quant=0.995, assay="RNA",
                                          chunk.size=1000, ncores=1, seed=123, verbose=TRUE, quiet=FALSE) {
  
  set.seed(seed)
  if (ncores>1) {
    require(future.apply)
    future_param_seed <<- seed
    future_param_ncores <<- ncores
  }
  DefaultAssay(ref) <- assay
  method <- method[1]
  
  species <- detect_species(ref)
  
  if (is.null(markers)) { #Default
    if (species=="human") {
      markers <- scGate_human_TIL_model@markers
    } else {
      markers <-  scGate_mouse_TIL_model@markers
    }
  }
  
  #Check selected cell type(s), and autoexpand if required
  sign <- names(markers)
  sign.names <- paste0(sign,"_scGate")
  names(sign.names) <- sign
  celltype.pass <- check_selected_celltypes(positive_celltypes, db=sign, autocomplete=autocomplete)
  
  if (length(celltype.pass)<1) {
    mess <- sprintf("Could not find selected cell types in marker list. Check your input to 'positive_celltypes' option")
    stop(mess)
  }  
  if (!quiet) {
    mess <- paste(celltype.pass, collapse=", ")
    message(sprintf("Training %s scGate model - positive set of signatures: %s", species, mess))
  }
  
  scores <- get_CTscores(obj=ref, markers.list=markers, z.score=FALSE,
                      method=method, chunk.size=chunk.size, ncores=ncores, maxRank=maxRank)
  

  ref_thr <- matrix(nrow=length(sign), ncol=2)
  rownames(ref_thr) <- sign
  colnames(ref_thr) <- c("mean","sd")
  
  for(s in sign){
    s.score <- sign.names[s]
    bulk <- scores[,s.score]
    bulk <- bulk[bulk >= quantile(bulk, p=1-quant) & bulk <= quantile(bulk,p=quant)]
    ref_thr[s,1] <- mean(bulk)
    ref_thr[s,2] <- ifelse(sd(bulk) > min.sd, sd(bulk), min.sd)
  }
  
  obj <- new("scGate_Model",
             species=species,
             positive_celltypes=celltype.pass,
             markers=markers,
             bg_model=as.data.frame(ref_thr),
             scoring_method=method)
  
  return(obj)
}


#' scGate_Model Class
#'
#' The scGate_Model class is use to store background models generated with the `scGate` package. It contains all the information needed
#' to apply the model to a query dataset and filter on the cell type(s) of interest stored in the `positive_celltypes` slot.
#'
#' @slot species The species of the training data
#' @slot positive_celltypes List of celltypes included in the reference set. These must be one or more from the list of markers (`markers` parameter).
#'     Only the `positive_celltypes` will be gated when the model is applied to a new query dataset
#' @slot markers List of markers for all celltypes included in the model
#' @slot bg_model A dataframe containing mean and SD calculated on training set for all cell types in `markers` list
#' @slot scoring_method The algorithm used to calculate signature scores     
#' @name scGate_Model-class
#' @rdname scGate_Model-class
#' @concept objects
#' @exportClass scGate_Model
#'
scGate_Model <- setClass(
  Class = "scGate_Model",
  slots = list(
    species = "character",
    positive_celltypes = "vector",
    markers = "list",
    bg_model = "data.frame",
    scoring_method = "character"
  )
)
