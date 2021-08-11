#' Filter single-cell data by cell type
#'
#' Apply scGate to filter specific cell types in a query dataset
#'
#' @param data Seurat object containing a query data set - filtering will be applied to this object
#' @param gating.model The background expected mean and SD for each cell type - see function \code{\link{train_scGate}}
#' @param max.impurity Maximum fraction of impurity allowed in clusters to flag cells as "pure".
#' @param max.impurity.decay Decrease by this fraction the max.impurity allowed at each iteration.
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
#'     The \code{active.ident} is also set to this variable. \cr
#'     The metadata field \code{scGate.annotation} will report the suggested identity of filtered cells based on the consensus maximum Zscores within clusters.
#'     If \code{return_signature_scores=TRUE}, Z-scores for all signatures in \code{markers} are added to the metadata of the Seurat object.
#' @note The parameters \code{max.impurity}, \code{sd.in} and \code{sd.out} also accept a vector of values. This allows specifying different values for this
#'     parameters for signatures of different levels
#' @examples
#' query <- scGate(query)
#' DimPlot(query)
#' DimPlot(query, group.by="scGate.annotation", label=T)
#' query <- subset(query, subset=is.pure=="Pure")
#' @seealso \code{\link{train_scGate()}} to generate celltype-specific models
#' @export

scGate <- function(data, gating.model=NULL, max.impurity=0.5,
                   max.impurity.decay=0,
                   ndim=30, resol=3, assay="RNA", 
                   sd.in=4, sd.out=7,
                   chunk.size=1000, ncores=1, maxRank=1500,
                   max.iterations=10, stop.iterations=0.01, min.cells=100,
                   additional.signatures=NULL,
                   genes.blacklist="Tcell.blacklist", 
                   seed=123, skip.normalize=FALSE, 
                   return_signature_scores=TRUE, verbose=FALSE, quiet=FALSE, compute_scores=T) {
  
  set.seed(seed)
  def.assay <- DefaultAssay(data)
  DefaultAssay(data) <- assay
  
  if (ncores>1) {
    require(future.apply)
    future_param_seed <<- seed
    future_param_ncores <<- ncores
  }
  
  species <- detect_species(data)
  
  if (is.null(gating.model)) { #Default
    gating.model <- scGate_DB[[species]]$Tcell
  }
  
  gating.model.list <- list()
  if (is.list(gating.model)) {
     gating.model.list <- gating.model
  } else { 
     gating.model.list[[1]] <- gating.model 
  }
  n.levels <- length(gating.model.list)
  
  #Vectorize parameters, to be able to assign different values to different levels
  parameters <- list("max.impurity"=max.impurity,
                     "sd.in"=sd.in,
                     "sd.out"=sd.out)
  parameters <- lapply(parameters, vectorize.parameters, lgt=n.levels)
  
  data$is.pure <- "Impure"
  data$scGate.annotation <- NA
  sub <- data
  for (lev in 1:n.levels) {
     message(sprintf("- Running scGate for level %i: ", lev))
    
     sub <- scGate_helper(data=sub,
                      gating.model=gating.model.list[[lev]],
                      max.impurity=parameters$max.impurity[lev],
                      max.impurity.decay = max.impurity.decay,
                      sd.in=parameters$sd.in[lev],
                      sd.out=parameters$sd.out[lev],
                      ndim=ndim, resol=resol, assay=assay,
                      chunk.size=chunk.size, species=species,
                      maxRank=maxRank, max.iterations=max.iterations,
                      stop.iterations=stop.iterations, min.cells=min.cells,
                      additional.signatures=additional.signatures,
                      genes.blacklist=genes.blacklist, skip.normalize=skip.normalize,
                      return_signature_scores=return_signature_scores,
                      ncores=ncores, quiet=quiet, verbose=verbose,compute_scores=compute_scores)
     
     if (return_signature_scores) {
        meta.cols <- grep("_scGate|_Zscore",colnames(sub@meta.data), perl=T, value = T)
        data <- AddMetaData(data, metadata = sub@meta.data[,meta.cols])  #NB: this will generate NAs for 2nd+ level signatures
     }
     data$scGate.annotation[colnames(sub)] <- sub$scGate.annotation
     if (sum(sub$is.pure=="Pure")==0) {
        sub <- NULL
        break   #all cells were removed
     } else {
        sub <- subset(sub, subset=is.pure=="Pure")
     }
  }
  pure.cells <- colnames(sub)
  data@meta.data[pure.cells, "is.pure"] <- "Pure"
  
  Idents(data) <- data$is.pure
  DefaultAssay(data) <- def.assay
  
  return(data)
} 

#' Generate a scGate model
#'
#' Given a reference set, calculate expected mean and SD for all cell types in the reference set. Deviations from this expected distribution (Z-score)
#' can then be used to gate for specific cell types in any query dataset (see function \code{scGate})
#'
#' @param ref Seurat object containing the reference data set
#' @param markers List of markers for each cell type, for example \code{scGate_DB$human$Tcell@@markers}
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
      markers <- scGate_DB[[species]]$Tcell@markers
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
    message(sprintf("Training %s scGate model - selected cell types: %s", species, mess))
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
