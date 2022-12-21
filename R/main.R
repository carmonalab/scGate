#' Filter single-cell data by cell type
#'
#' Apply scGate to filter specific cell types in a query dataset
#'
#' @param data Seurat object containing a query data set - filtering will be applied to this object
#' @param model A single scGate model, or a list of scGate models. See Details for this format
#' @param pca.dim Number of dimensions for cluster analysis
#' @param assay Seurat assay to use
#' @param slot Data slot in Seurat object
#' @param pos.thr Minimum UCell score value for positive signatures
#' @param neg.thr Maximum UCell score value for negative signatures
#' @param maxRank Maximum number of genes that UCell will rank per cell
#' @param nfeatures Number of variable genes for dimensionality reduction
#' @param k.param Number of nearest neighbors for knn smoothing
#' @param reduction Dimensionality reduction to use for knn smoothing. By default, calculates a new reduction
#'     based on the given \code{assay}; otherwise you may specify a precalculated dimensionality reduction (e.g.
#'     in the case of an integrated dataset after batch-effect correction)
#' @param pca.dim Number of principal components for dimensionality reduction
#' @param param_decay Controls decrease in parameter complexity at each iteration, between 0 and 1.
#'     \code{param_decay == 0} gives no decay, increasingly higher \code{param_decay} gives increasingly stronger decay
#' @param ncores Number of processors for parallel processing
#' @param output.col.name Column name with 'pure/impure' annotation
#' @param min.cells Minimum number of cells to cluster or define cell types
#' @param additional.signatures A list of additional signatures, not included in the model, to be evaluated (e.g. a cycling signature). The scores for this
#'     list of signatures will be returned but not used for filtering.
#' @param save.levels Whether to save in metadata the filtering output for each gating model level
#' @param keep.ranks Store UCell rankings in Seurat object. This will speed up calculations if the same object is applied again with new signatures.
#' @param genes.blacklist Genes blacklisted from variable features. The default loads the list of genes in \code{scGate::genes.blacklist.default};
#'     you may deactivate blacklisting by setting \code{genes.blacklist=NULL}
#' @param multi.asNA How to label cells that are "Pure" for multiple annotations: "Multi" (FALSE) or NA (TRUE)
#' @param seed Integer seed for random number generator
#' @param verbose Verbose output

#' @return A new metadata column \code{is.pure} is added to the query Seurat object, indicating which cells passed the scGate filter.
#'     The \code{active.ident} is also set to this variable.
#' @details Models for scGate are data frames where each line is a signature for a given filtering level.
#'     A database of models can be downloaded using the function \code{get_scGateDB}.
#'     You may directly use the models from the database, or edit one of these models to generate your own custom gating model.
#'     
#'     Multiple models can also be evaluated at once, by running scGate with a list of models. Gating for each individual model is
#'     returned as metadata, with a consensus annotation stored in \code{scGate_multi} metadata field. This allows using scGate as a
#'     multi-class classifier, where only cells that are "Pure" for a single model are assigned a label, cells that are "Pure" for
#'     more than one gating model are labeled as "Multi", all others cells are annotated as NA.
#' @examples
#' ### Test using a small toy set
#' data(query.seurat)
#' # Define basic gating model for B cells
#' my_scGate_model <- gating_model(name = "Bcell", signature = c("MS4A1")) 
#' query.seurat <- scGate(query.seurat, model = my_scGate_model, reduction="pca")
#' table(query.seurat$is.pure)
#' \donttest{
#' ### Test with larger datasets
#' library(Seurat)
#' testing.datasets <- get_testing_data(version = 'hsa.latest')
#' seurat_object <- testing.datasets[["JerbyArnon"]]
#' # Download pre-defined models
#' models <- get_scGateDB()
#' seurat_object <- scGate(seurat_object, model=models$human$generic$PanBcell)
#' DimPlot(seurat_object)
#' seurat_object_filtered <- subset(seurat_object, subset=is.pure=="Pure")
#'
#' ### Run multiple models at once
#' models <- get_scGateDB()
#' model.list <- list("Bcell" = models$human$generic$Bcell,
#'                    "Tcell" = models$human$generic$Tcell)
#' seurat_object <- scGate(seurat_object, model=model.list)
#' DimPlot(seurat_object, group.by = "scGate_multi")
#' }
#' @seealso \code{\link{load_scGate_model}} \code{\link{get_scGateDB}} \code{\link{plot_tree}}
#' @import Seurat
#' @import ggplot2
#' @importFrom dplyr %>% distinct bind_rows
#' @importFrom UCell AddModuleScore_UCell SmoothKNN
#' @import BiocParallel
#' @export

scGate <- function(data,
                   model,
                   pos.thr=0.2,
                   neg.thr=0.2,
                   assay=NULL,
                   slot="data",
                   ncores=1,
                   seed=123,
                   keep.ranks=FALSE,
                   reduction=c("calculate","pca","umap","harmony","Liors_elephant"),
                   min.cells=30,
                   nfeatures=2000,
                   pca.dim=30,
                   param_decay=0.25,
                   maxRank=1500,
                   output.col.name='is.pure',
                   k.param=30,
                   genes.blacklist="default",
                   multi.asNA = FALSE,
                   additional.signatures=NULL,
                   save.levels=FALSE,
                   verbose=FALSE) {
  
  set.seed(seed)
  
  if (!is.null(assay)) {
    DefaultAssay(data) <- assay
  }
  assay <- DefaultAssay(data)
  
  if (assay == "integrated") { #UCell should not run on integrated assay
    if ('RNA' %in% Assays(data)) {
      assay.ucell <- 'RNA'
    } else if ('SCT' %in% Assays(data)) {
      assay.ucell <- 'SCT'
    } else {
      stop("Cannot find assays with unintegrated data in this Seurat object")
    }
  } else {
    assay.ucell <- assay
  }
  
  reduction <- reduction[1]
  if (is.null(reduction) || tolower(reduction)=="calculate") {
    reduction = "calculate"
  } else {
    if (!reduction %in% Reductions(data)) {
      stop(sprintf("Could not find reduction %s in this object. Set reduction='calculate' to compute a new dimred", reduction))
    }
    pca.dim <- ncol(data@reductions[[reduction]])
  }
    
  #check gene blacklist
  if (!is.null(genes.blacklist)) {
    if (length(genes.blacklist)==1 && genes.blacklist == "default") {  #Default
       genes.blacklist = scGate::genes.blacklist.default
    }  
    if (is.list(genes.blacklist)) {
      genes.blacklist <- unlist(genes.blacklist)
    }
    genes.blacklist <- unique(genes.blacklist)
  }
  
  #With single model, make a list of one element
  if (!inherits(model, "list")) {
    model <- list("Target" = model)
  }
  if (is.null(names(model))) {
    names(model) <- paste(output.col.name, seq_along(model), sep = ".")
  }
  
  if (ncores>1) {
    bpp <- MulticoreParam(workers=ncores)
  } else {
    bpp <- SerialParam()
  }
  
  # compute signature scores using UCell
  if (verbose) {
    message(sprintf("Computing UCell scores for all signatures using %s assay...\n", assay.ucell))
  }
  data <- score.computing.for.scGate(data, model, bpp=bpp, assay=assay.ucell,
                                     slot=slot, maxRank=maxRank, 
                                     keep.ranks=keep.ranks,
                                     add.sign=additional.signatures)
  
  for (m in names(model)) {
    
    col.id <- paste0(output.col.name, "_", m)

    data <- run_scGate_singlemodel(data, model=model[[m]], k.param=k.param,
                           param_decay=param_decay, pca.dim=pca.dim,
                           nfeatures=nfeatures, min.cells=min.cells, bpp=bpp,
                           assay=assay, slot=slot, genes.blacklist=genes.blacklist,
                           pos.thr=pos.thr, neg.thr=neg.thr, verbose=verbose,
                           reduction=reduction, colname=col.id, save.levels=save.levels)
    
    Idents(data) <- col.id
    n_pure <- sum(data[[col.id]]=="Pure")
    frac.to.keep <- n_pure/ncol(data)
    mess <- sprintf("\n### Detected a total of %i pure '%s' cells (%.2f%% of total)",
                    n_pure, m, 100*frac.to.keep)
    message(mess)
  }
 
  #Combine results from multiple model into single cell type annotation 
  data <- combine_scGate_multiclass(data, prefix=paste0(output.col.name,"_"),
                            scGate_classes = names(m), multi.asNA = multi.asNA,
                            min_cells=min.cells, out_column = "scGate_multi")

  #Back-compatibility with previous versions
  if (names(model)[1] == 'Target') {
    cn <- paste0(output.col.name, "_Target")
    data@meta.data[,output.col.name] <- data@meta.data[,cn]
    data@meta.data[,cn] <- NULL
    
    if (save.levels) {
      for (l in unique(model[[1]]$levels)) {
        cn <- paste0(output.col.name, "_Target.",l)
        data@meta.data[,paste0(output.col.name,".",l)] <- data@meta.data[,cn]
        data@meta.data[,cn] <- NULL
      }
    }
  }
  return(data)
}


#' Plot model tree
#'
#' View scGate model as a decision tree (require ggparty package)
#'
#' @param model A scGate model to be visualized
#' @param box.size Box size
#' @param edge.text.size Edge text size
#' @return A plot of the model as a decision tree. At each level, green boxes
#'     indicate the 'positive' (accepted) cell types, red boxed indicate the
#'     'negative' cell types (filtered out). The final Pure population is the
#'     bottom right subset in the tree.
#' @examples
#' library(ggparty)
#' models <- get_scGateDB()
#' plot_tree(models$human$generic$Tcell)
#' @export


plot_tree <- function(model, box.size = 8, edge.text.size = 4) {
  
  if (!requireNamespace('ggparty', quietly = TRUE)) {  #check whether ggparty is available
    stop("Please install and load package 'ggparty'")
  }
  nlev <- length(unique(model$levels))
  if(nlev <= 1){
    stop("your model must contain at least two levels to be ploted as a tree")
  }
  #restructure data for visualization
  level.list <- list()
  for (i in 1:nlev) {
    level.list[[i]] <- list()
    sub <- subset(model, tolower(model$levels)==paste0("level",i))
    
    level.list[[i]][["positive"]] <- sub[sub$use_as=="positive","name"]
    level.list[[i]][["negative"]] <- sub[sub$use_as=="negative","name"]
  }
  
  #Initialize dataframe for tree
  df <- data.frame(matrix(ncol=nlev+1, nrow=nlev+1, data = 0))
  colnames(df) <- c(paste0("Level_", 1:nlev), "Pure")
  
  for (i in 2:nlev) {
    for (j in 1:(i-1)) {
      df[i,j] <- 1   
    }
  }
  df[nlev+1,] <- 1
  
  ##Construct tree structure
  pn <- list()
  #bottom level
  
  pn[[nlev]] <- partykit::partynode(nlev+1,
                          split = partykit::partysplit(nlev, index=1:2, breaks = 0),
                          kids = list(partykit::partynode(nlev+2),
                                      partykit::partynode(nlev+3))) 
  
  for (i in (nlev-1):1) {
    pn[[i]] <- partykit::partynode(i,
                         split = partykit::partysplit(i, index=1:2, breaks=0),
                         kids = list(partykit::partynode(i+1),
                                     pn[[i+1]]))
  }
  
  #first element in list has complete structure
  py <- partykit::party(pn[[1]], df)
  
  
  sign.annot <- vector(length=2*nlev+1)
  is.pos <- vector(length=2*nlev+1)
  sign.annot[1] <- "Root"
  is.pos[1] <- NA
  
  for (i in 1:nlev) {
    sign.annot[2*i] <- paste0(level.list[[i]]$negative, collapse = "\n")
    sign.annot[2*i+1] <- paste0(level.list[[i]]$positive, collapse = "\n")
    
    is.pos[2*i] <- "Negative"
    is.pos[2*i+1] <- "Positive"
  }
  
  gg <- ggparty::ggparty(py)
  gg$data$info <- sign.annot
  gg$data$p.value <- is.pos
  
  
  gg$data$breaks_label[grep("<=", gg$data$breaks_label)] <- "Negative"
  gg$data$breaks_label[grep(">", gg$data$breaks_label)] <- "Positive"
  
  gg <- gg + ggparty::geom_edge() +
    ggparty::geom_edge_label(size = edge.text.size) +
    ggparty::geom_node_label(ids = "inner",
                    mapping = aes(col = .data$p.value),
                    line_list = list(aes(label= .data$info)),
                    line_gpar = list(list(size = box.size)))  +
    ggparty::geom_node_label(ids = "terminal",
                    mapping = aes(col = .data$p.value),
                    nudge_y=0.01,
                    line_list = list(aes(label= .data$info)),
                    line_gpar = list(list(size = box.size))) +
    scale_color_manual(values=c("#f60a0a", "#00ae60")) +
    theme(legend.position = "none", plot.margin = unit(c(1,1,1,1), "cm")) 
  
  return(gg)
}

#' Load a single scGate model
#'
#' Loads a custom scGate model into R. For the format of these models, have a
#' look or edit one of the default models obtained with \code{\link{get_scGateDB}}
#'
#' @param model_file scGate model file, in .tsv format.
#' @param master.table File name of the master table (in repo_path folder) that contains cell type signatures.
#' @return A scGate model in dataframe format, which can given as input to the \code{\link{scGate}} function.
#' @examples
#' dir <- tempdir() # this may also be set to your working directory
#' models <- get_scGateDB(destination=dir)
#' # Original or edited model
#' model.path <- paste0(dir,"/scGate_models-master/human/generic/Bcell_scGate_Model.tsv")
#' master.path <- paste0(dir,"/scGate_models-master/human/generic/master_table.tsv")
#' my.model <- load_scGate_model(model.path, master.path)
#' my.model
#' @seealso \code{\link{scGate}} \code{\link{get_scGateDB}} 
#' @importFrom utils read.table
#' @export

load_scGate_model <- function(model_file, master.table = "master_table.tsv"){
  
  model <- read.table(model_file, sep ="\t",header =TRUE)
  model <- use_master_table(model, master.table = master.table)
  
  return(model)
}


#' Load scGate model database
#'
#' Download, update or load local version of the scGate model database. These are stored in a GitHub repository, from where you can download specific 
#' versions of the database.
#' 
#' @param destination Destination path for storing the DB. The default is tempdir(); if you wish to edit locally the models and
#'    link them to the current project, set this parameter to a new directory name, e.g. scGateDB
#' @param force_update  Whether to update an existing database.
#' @param version Specify the version of the scGate_models database (e.g. 'v0.1'). By default downloads the latest available version.
#' @param repo_url  URL path to scGate model repository database
#' @param branch  branch of the scGate model repository, either 'master' (default) or 'dev' for the latest models 
#' @param verbose  display progress messages
#' @return A list of models, organized according to the folder structure of the database. See the examples below.
#' @details Models for scGate are dataframes where each line is a signature for a given filtering level. A database of models can be downloaded using the function
#'     \code{get_scGateDB}. You may directly use the models from the database, or edit one of these models to generate your own custom gating model.  
#' @examples
#' scGate.model.db <- get_scGateDB()
#' # To see a specific model, browse the list of models:
#' scGate.model.db$human$generic$Myeloid
#' # Apply scGate with this model
#' data(query.seurat)
#' query <- scGate(query.seurat, model=scGate.model.db$human$generic$Myeloid, reduction="pca")
#' @seealso \code{\link{scGate}} \code{\link{load_scGate_model}}
#' @importFrom dplyr %>%  
#' @importFrom utils download.file unzip read.table
#' @export

get_scGateDB <- function(destination = tempdir(),
                         force_update = FALSE,
                         version = "latest",
                         branch=c("master","dev"), 
                         verbose=FALSE,
                         repo_url = "https://github.com/carmonalab/scGate_models"){
  
  branch = branch[1]
  if (version == "latest") {
    repo_url_zip = sprintf("%s/archive/%s.zip", repo_url,branch)
    repo.name <- paste0("scGate_models-",branch)
    repo.name.v <- repo.name
  } else {
    repo_url_zip = sprintf("%s/archive/refs/tags/%s.zip", repo_url, version)
    repo.name = sprintf("scGate_models-%s", version)
    #for some reason GitHub remove the 'v' from repo name after unzipping
    repo.name.v <- sprintf("scGate_models-%s", gsub("^v","",version, perl=TRUE)) 
  }
  destination <- normalizePath(destination, winslash = "/")
  repo_path = file.path(destination,repo.name)
  repo_path.v = file.path(destination,repo.name.v)
  temp <- tempfile()
  
  if(!dir.exists(repo_path)){
    if(!dir.exists(destination)) {
      dir.create(destination)
    }
    download.file(repo_url_zip,temp)
    unzip(temp,exdir = destination)
    unlink(temp)
  }else if(force_update){
    download.file(repo_url_zip,temp)
    system(sprintf("rm -r %s",repo_path))  # this ensure that db would be completely overwritten and old model will not persist. 
    unzip(temp,exdir = destination, overwrite = force_update)
    unlink(temp)
  }else{
    message(sprintf("Using local version of repo %s. If you want update it, set option force_update = TRUE",repo.name))
  }
  
  #Now load the models into a list structure
  allfiles <- list.files(repo_path.v, recursive = TRUE)
  modelfiles <- grep("scGate_Model.tsv", allfiles, value = TRUE)
  uniq_dirs <- sort(unique(dirname(modelfiles)))
  
  model_db <- list()
  for (dir in uniq_dirs) {
    sub <- strsplit(dir, split="/")[[1]]
    model_path <- file.path(repo_path.v, dir)
    
    if (length(sub)==0) {
      stop("Error in scGate DB format")
    } else if (length(sub)==1) {
      if(verbose) message(paste("loading ",model_path))
      model_db[[sub[1]]] <- load.model.helper(model_path,verbose=verbose)
    } else if (length(sub)==2) {
      if(verbose) message(paste("loading ",model_path))
      model_db[[sub[1]]][[sub[2]]] <- load.model.helper(model_path,verbose=verbose)
    } else if (length(sub)==2) {
      if(verbose) message(paste("loading ",model_path))
      model_db[[sub[1]]][[sub[2]]][[sub[[3]]]] <- load.model.helper(model_path,verbose=verbose)
    } else {
      message(sprintf("Warning: max depth for scGate models is 3. Skipping folder %s", model_path))
    } 
  }
  return(model_db)
}

#' Performance metrics
#'
#' Evaluate model performance for binary tasks
#' 
#' @param actual Logical or numeric binary vector giving the actual cell labels. 
#' @param pred  Logical or numeric binary vector giving the predicted cell labels. 
#' @param return_contingency  Logical indicating if contingency table must be returned.
#' @return Prediction performance metrics (Precision, Recall, MCC) between actual and
#'    predicted cell type labels.
#' @examples
#' results <- performance.metrics(actual= sample(c(1,0),20,replace=TRUE),
#'     pred =  sample(c(1,0),20,replace=TRUE,prob = c(0.65,0.35) ) )
#' @export

performance.metrics <- function(actual,pred,return_contingency=FALSE){
  actual <- as.numeric(actual +0)  
  pred <- as.numeric(pred +0)  
  tp <- sum(actual&pred)
  tn <- sum((!actual)&(!pred))
  fn <- sum(actual&(!pred))
  fp <- sum((!actual)&pred)  
  
  PREC <- tp/(tp +fp)
  REC <- tp/(tp + fn)
  #sqrt_ <- sqrt((tp + fp)*(tp+fn)*(tn+fp)*(tn+fn))
  sqrt_ <- exp(0.5* sum(log(c(tp+fp, tp+fn, tn+fp, tn+fn))) )
  MCC <- (tp*tn - fp*fn) / sqrt_
  
  
  if(!return_contingency){  
    res.Summary <- c(PREC,REC,MCC); names(res.Summary) <- c("PREC","REC","MCC")
    return(res.Summary)
  }else{
    ct <- table(actual,pred)
    ## ordering contingency table, but avoiding errors when all predictions (or all actual cells) are equals
    nam.act <- unique(actual)%>%sort(decreasing = TRUE)%>%as.character()  # 
    nam.pred <- unique(pred)%>%sort(decreasing = TRUE)%>%as.character()
    ct <- ct[nam.act,nam.pred]  
    return(list('counting' = ct,
                'summary' = res.Summary ))
  }
  
}

#' Test your model
#'
#' Wrapper for fast model testing on 3 sampled datasets 
#' 
#' @param model scGate model in data.frame format 
#' @param testing.version  Character indicating the version of testing tatasets
#'     to be used. By default "hsa-latest" will be used. It will be ignored if
#'     a custom dataset is provided (in Seurat format). 
#' @param custom.dataset  Seurat object to be used as a testing dataset. For
#'     testing purposes, metadata seurat object must contain a column named
#'     'cell_type' to be used as a gold standard. Also a set of positive
#'     targets must be provided in the target variable. 
#' @param target Positive target cell types. If default testing version is used
#'     this variable must be a character indicating one of the available target
#'     models ('immune','Lymphoid','Myeloid','Tcell','Bcell','CD8T','CD4T',
#'     'NK','MoMacDC','Plasma_cell','PanBcell'). 
#'     If a custom dataset is provided in Seurat format, this variable must be
#'     a vector of positive cell types in your data. The last case also require
#'     that such labels were named as in your cell_type meta.data column. 
#' @param plot Whether to return plots to device
#' @return Returns performance metrics for the benchmarking datasets, and optionally
#'     plots of the predicted cell type labels in reduced dimensionality space.  
#' @examples
#' \donttest{
#' scGate.model.db <- get_scGateDB()
#' # Browse the list of models and select one:
#' model.panBcell <-  scGate.model.db$human$generic$PanBcell
#' # Test the model with available testing datasets
#' panBcell.performance <- test_my_model(model.panBcell, target = "PanBcell")
#' model.Myeloid <-  scGate.model.db$human$generic$Myeloid
#' myeloid.performance <- test_my_model(model.Myeloid, target = "Myeloid")
#' }     
#' @importFrom utils download.file
#' @importFrom methods is
#' @importFrom patchwork wrap_plots
#' @export

test_my_model <- function(model, testing.version = 'hsa.latest',
                          custom.dataset = NULL,target = NULL,
                          plot = TRUE){
  
  performance.computation  <- ifelse (is.null(target), FALSE, TRUE)
  
  if (is(custom.dataset, "Seurat")){
    testing.datasets <- list()
    testing.datasets$user.dataset <- custom.dataset
    custom <- TRUE
  } else { 
    custom <- FALSE
  }

  if(!custom){
    targets <- c('immune','Lymphoid','Myeloid','Tcell','Bcell','CD8T','CD4T','NK','MoMacDC','Plasma_cell','PanBcell')
    
    if(is.null(target)){
      message("warning: target cell_type not provided. Avoiding performance computation")  
      performance.computation <- FALSE
    }else if(!target %in% targets){
      stop(sprintf("target must be one of %s; or NULL for avoiding performance computation",paste(targets,collapse = "';'")))
    }
    
    ## check dataset version
    available.datasets = c("hsa.latest")
    if(!testing.version %in% available.datasets){
      stop("Please provide a valid testing.version parameter or provide a custom.dataset in seurat format")
    }
    
    # load testing datasets
    if(testing.version == "hsa.latest"){
      testing.datasets <- get_testing_data(version = testing.version)
    }
  }  
  
  if(custom){
    if(!"cell_type" %in% colnames(custom.dataset@meta.data)){
      stop("please, provide a 'cell_type' column to be used as reference cell type")
    }
    
    if(is.null(target)){
      message("warning: target cell_type not provided. Avoiding performance computation")  
      performance.computation <- FALSE
    }else if(any(!target %in% custom.dataset$cell_type)){
      stop("all target celltypes must be included in cell_type metadata field. Otherwise, set target = NULL for avoiding performance computation")
    }
  }
  
  plt.out <- list()
  perf.out <- list()
  output <- list()
  # Drop is.pure cols if exists
  for(dset in names(testing.datasets)){
    obj <- testing.datasets[[dset]]
    plt <- list()
    cols <- colnames(obj@meta.data)
    dropcols = grep("^is.pure",cols,value =TRUE) %>% unique()
    if(length(dropcols)>0){
      for(col in dropcols){
        obj[[col]] <- NULL   
      }
    }
    
    ## scGate filtering
    obj <- scGate(obj, model = model, assay = DefaultAssay(obj))

    # add annotation plot
    nname <- sprintf("%s manual annot",dset)
    plt <- DimPlot(obj, group.by = "cell_type", label = TRUE,
                   repel =TRUE, label.size = 3) + 
      ggtitle(nname) + NoLegend() +  theme(aspect.ratio = 1)

    # add one DimPlot by model level
    pure.plot <- DimPlot(obj, group.by = "is.pure", cols = list("Pure"="green","Impure"="gray")) +
      theme(aspect.ratio = 1)
    plt <- list("Annotation"=plt, "Gating"=pure.plot)
    
    #reserve plots of this dset
    plt.out[[dset]] <- patchwork::wrap_plots(plt,ncol = length(plt))
    
    if(performance.computation){
      if(!custom){    
        performance = scGate::performance.metrics(actual = obj@meta.data[,target],
                                                  pred = obj$`is.pure`== "Pure")
      }else{
        performance = scGate::performance.metrics(actual = obj@cell_type %in% target,
                                                  pred = obj$`is.pure`== "Pure")
      }
      perf.out[[dset]] <- performance 
    }
    output[[dset]] <- obj
    
  }
  
  if(performance.computation){
    perf <- Reduce(rbind,perf.out)
    rownames(perf) <- names(perf.out)
  }
  
  if(plot) {
    print(patchwork::wrap_plots(plt.out, ncol = 1))
  }

  if(performance.computation){
    return(list(performance = perf, plots = plt.out, objects = output))
  }else{
    return(list(plots = plt.out, objects = output))
  }
}

#' Plot scGate filtering results by level
#'
#' Fast plotting of gating results over each model level.
#' 
#' @param obj Gated Seurat object output of scGate filtering function
#' @param pure.col Color code for pure category 
#' @param impure.col Color code for impure category
#' @return UMAP plots with 'Pure'/'Impure' labels for each level of the scGate model
#' @examples
#' scGate.model.db <- get_scGateDB()
#' model <- scGate.model.db$human$generic$Myeloid
#' # Apply scGate with this model
#' data(query.seurat)
#' query.seurat <- scGate(query.seurat, model=model,
#'     reduction="pca", save.levels=TRUE)
#' library(patchwork)     
#' pll <- plot_levels(query.seurat)
#' wrap_plots(pll)
#' @importFrom Seurat DimPlot
#' @export

plot_levels <- function(obj, pure.col = "green" ,impure.col = "gray"){
  myCols <- grep("^is.pure.", colnames(obj@meta.data),value = TRUE)
  plots <- list()
  for (myCol in myCols){
    plots[[myCol]] <- DimPlot(obj, group.by = myCol, 
                              cols = list(Pure = pure.col,Impure = impure.col)) +
      theme(aspect.ratio = 1)
  }
  return(plots)
}

#' Plot UCell scores by level
#'
#' Show distribution of UCell scores for each level of a given scGate model
#' 
#' @param obj Gated Seurat object (output of scGate)
#' @param model scGate model used to identify a target population in obj
#' @param pos.thr Threshold for positive signatures used in scGate model (set to NULL to disable)
#' @param neg.thr Threshold for negative signatures used in scGate model (set to NULL to disable) 
#' @param overlay Degree of overlay for ggridges
#' @param ncol Number of columns in output object (passed to wrap_plots)
#' @param combine Whether to combine plots into a single object, or to return a list of plots
#' @return Returns a density plot of UCell scores for the signatures in the scGate model,
#'     for each level of the model  
#' @examples
#' scGate.model.db <- get_scGateDB()
#' model <- scGate.model.db$human$generic$Tcell
#' # Apply scGate with this model
#' data(query.seurat)
#' query.seurat <- scGate(query.seurat, model=model,
#'     reduction="pca", save.levels=TRUE)
#' # View UCell score distribution
#' plot_UCell_scores(query.seurat, model)
#' @return Either a plot combined by patchwork (combine=T) or a list of plots (combine=F)
#' @importFrom reshape2 melt
#' @importFrom ggridges geom_density_ridges
#' @importFrom patchwork wrap_plots
#' @export
#' 
plot_UCell_scores <- function(obj, model, overlay=5, pos.thr=0.2,
                              neg.thr=0.2, ncol=NULL, combine=TRUE) {
  
  u_cols <- grep('_UCell', colnames(obj@meta.data), value = TRUE)
  
  levs <- unique(model$levels)
  pll <- list()
  
  palette <- c("#00fd0c","#f4340e")
  names(palette) <- c("Positive","Negative")
  
  if (sum(grepl("is.pure.level", colnames(obj@meta.data)))==0) {
     obj$is.pure.level1 <- obj$is.pure
     
     if (length(levs)>1) {
       warning("scGate levels were not stored in this object. Showing results only for top level.")
       levs <- "level1"
     }
  }
  
  for (l in seq_along(levs)) {
    
    lev.name <- levs[l]
    
    sigs <- model[model$levels == lev.name, c("use_as","name")]
    
    col <- sprintf("%s_UCell", sigs$name)
    col <- col[col %in% u_cols]
    
    meta <- obj@meta.data
    if (l>1) {
      meta <- meta[meta[sprintf("is.pure.level%i",l-1)]=="Pure",]
    }
    ncells <- nrow(meta)
    stat <- table(meta[,sprintf("is.pure.level%i",l)])
    
    to.plot <- meta[,col, drop=FALSE]
    colnames(to.plot) <- gsub("_UCell","",colnames(to.plot))
    
    to.plot <- reshape2::melt(to.plot, id=NULL)
    colnames(to.plot) <- c("Signature","Score")
    
    to.plot$Class <- "Positive"
    to.plot$Class[to.plot$Signature %in% sigs[sigs$use_as =="negative","name"]] <- "Negative"
    
    #vertical lines (thresholds)
    to.plot$Thr <- NA
    if (!is.null(pos.thr)) {
       to.plot[to.plot$Class=="Positive","Thr"] <- pos.thr
    }
    if (!is.null(neg.thr)) {
      to.plot[to.plot$Class=="Negative","Thr"] <- neg.thr
    }
    
    #Make ggridges distribution plot
    pll[[l]] <- ggplot(to.plot, aes(x =.data$Score, y =.data$Signature, fill=.data$Class)) + 
      geom_density_ridges(scale = overlay) +
      scale_fill_manual(values = palette) + theme_minimal() +
      theme(axis.title.y=element_blank()) + ggtitle(sprintf("%s - %i/%i pure ",lev.name, stat["Pure"], ncells)) 
    
    #Add threshold lines
    if (!is.null(pos.thr) |  !is.null(neg.thr)) {  
      pll[[l]] <- pll[[l]] + geom_segment(aes(x = .data$Thr, xend = .data$Thr,
                                              y = as.numeric(.data$Signature), 
                                              yend = as.numeric(.data$Signature)+0.9), linetype = "dashed")
    }
  }
  #Return combined plot or list of plots
  if (combine) {
    return(wrap_plots(pll, ncol=ncol))
  } else {
    return(pll)
  }
}

#' Model creation and editing
#'
#' Generate an scGate model from scratch or edit an existing one
#' 
#' @param model scGate model to be modified. When is NULL (default) a new model will be initialized.   
#' @param level integer. It refers to the hierarchical level of the model tree in which the signature will be added (level=1 by default)    
#' @param name Arbitrary signature name (i.e. Immune, Tcell, NK etc).   
#' @param signature character vector indicating gene symbols to be included in the signature (e.g. CD3D). If a minus sign is placed to the end of a gene name (e.g. "CD3D-"), this gene will be used as negative in UCell computing. See UCell documentation for details    
#' @param positive Logical indicating if the signature must be used as a positive signature in those model level. Default is TRUE. 
#' @param negative Same as `positive` but negated (negative=TRUE equals to positive=FALSE)
#' @param remove Whether to remove the given signature from the model
#' @return A scGate model that can be used by \code{\link{scGate}} to filter target cell types.
#' @examples
#' # create a simple gating model
#' my_model <- gating_model(level = 1, name = "immune", signature = c("PTPRC"))
#' my_model <- gating_model(model = my_model, level = 1, positive = FALSE,
#'     name = "Epithelial", signature = c("CDH1","FLT1") )
#' # Remove an existing signature
#' dropped_model <- gating_model(model = my_model, remove =TRUE, level = 1, name = "Epithelial")
#' @importFrom stats setNames
#' @export

gating_model <- function(model=NULL, level= 1, name, signature,
                         positive = TRUE, negative = FALSE, remove = FALSE){
  
  template <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("levels","use_as", "name", "signature"))
  
  if(negative){
    positive <- FALSE
  }
  
  if(is.null(model)){
    model <- template 
  }
  
  if(!remove){
    new.signature <- data.frame(levels = paste0("level",level),
                                use_as = ifelse(positive, "positive","negative"),
                                name = name,
                                signature = ifelse(length(signature) >1, paste(signature,collapse = ";") ,signature))
    model <- rbind(model,new.signature)
  }else{
    lev <- paste0("level",level)
    
    model <- model[!((model$levels == lev) & (model$name == name)),]
  }
  return(model)
}

#' Download sample data
#'
#' Helper function to obtain some sample data
#' 
#' @param version Which sample dataset   
#' @param destination Save to this directory
#' @return A list of datasets that can be used to test scGate    
#' @examples
#' \donttest{
#' testing.datasets <- get_testing_data(version = 'hsa.latest')
#' }
#' @export

get_testing_data <- function(version = 'hsa.latest', destination = tempdir()){
  data.folder = file.path(destination,"testing.data")
  if(!dir.exists(data.folder)){
    dir.create(data.folder,recursive = TRUE)
  }
  if(version == 'hsa.latest'){
    testing.data.url = "https://figshare.com/ndownloader/files/31114669?private_link=75b1193bd4c705ffb50b"
    testing.data.path = file.path(data.folder,"testing.dataset.2k.rds")
  }
  if(!file.exists(testing.data.path)){
      download.file(url = testing.data.url,destfile = testing.data.path)
  }
  
  testing.data <- readRDS(testing.data.path)
  return(testing.data)
}

#' Combine scGate annotations
#'
#' If a single-cell dataset has precomputed results for multiple scGate models, combined them in multi-class annotation
#' 
#' @param obj Seurat object with scGate results for multiple models stored as metadata   
#' @param prefix Prefix in metadata column names for scGate result models
#' @param scGate_classes Vector of scGate model names. If NULL, use all columns that start with "prefix" above.
#' @param min_cells Minimum number of cells for a cell label to be considered
#' @param multi.asNA How to label cells that are "Pure" for multiple annotations: "Multi" (FALSE) or NA (TRUE)
#' @param out_column The name of the metadata column where to store the multi-class cell labels
#' @return A Seurat object with multi-class annotations based on the combination of multiple models. A new
#'      column (by default "scGate_multi") is added to the metadata of the Seurat object.
#' @import Seurat
#' @examples
#' \donttest{
#' # Define gating models
#' model.B <- gating_model(name = "Bcell", signature = c("MS4A1")) 
#' model.T <- gating_model(name = "Tcell", signature = c("CD2","CD3D","CD3E"))
#' # Apply scGate with these models
#' data(query.seurat)
#' query.seurat <- scGate(query.seurat, model=model.T,
#'     reduction="pca", output.col.name = "is.pure_Tcell")
#' query.seurat <- scGate(query.seurat, model=model.B,
#'     reduction="pca", output.col.name = "is.pure_Bcell")
#' query.seurat <- combine_scGate_multiclass(query.seurat, scGate_class=c("Tcell","Bcell"))      
#' table(query.seurat$scGate_multi)
#' }
#' @export

combine_scGate_multiclass <- function(obj,
                                  prefix="is.pure_",
                                  scGate_classes=NULL,
                                  min_cells=20,
                                  multi.asNA = FALSE,
                                  out_column="scGate_multi"
){
  #Use all columns with given prefix
  if (is.null(scGate_classes)){  
    cols <- grep(prefix, colnames(obj@meta.data), value = TRUE)
    cols <- grep("\\.level\\d+$", cols, invert=TRUE, perl=TRUE, value=TRUE)
  } else {
    cols <- paste0(prefix, scGate_classes)
    cols <- cols[cols %in% colnames(obj@meta.data)]
  }
  if (is.null(cols)) {
    stop("Could not find scGate annotations in this object metadata.")
  }
  
  meta <- obj@meta.data[,cols, drop=FALSE]
  meta[is.na(meta)] <- "Impure"  #Avoid NAs
  obj.logical <- meta=="Pure"
  
  label.sums <- apply(obj.logical,1,sum)
  
  obj.single <- obj.logical[label.sums==1, , drop=FALSE]
  obj.single.labels <- apply(obj.single,1,function(x) names(x)[x])
  #remove prefix
  if (!is.null(prefix)) {
    obj.single.labels <- gsub(prefix, "", obj.single.labels)
  }
  
  #Assign labels to uniquely identified cells
  labs <- rep(NA, ncol(obj))
  names(labs) <- colnames(obj)
  
  labs[names(obj.single.labels)] <- obj.single.labels
 
  #Set to NA classes with too few cells
  tt <- table(labs, useNA = "always")
  labs[labs %in% names(tt)[tt<min_cells]] <- NA
  
  if (multi.asNA) {
    labs[names(label.sums[label.sums>1])] <- NA
  } else { 
    labs[names(label.sums[label.sums>1])] <- "Multi"
  }
  
  obj@meta.data[[out_column]] <- labs
  obj
}

#' Blocklist of genes for dimensionality reduction
#' 
#' A list of signatures, for mouse and human. These include cell cycling,
#' heat-shock genes, mitochondrial genes, and other genes classes, that may
#' confound the identification of cell types. These are used internally by scGate
#' and excluded from the calculation of dimensional reductions (PCA).
#' 
#' @name genes.blacklist.default
#' @docType data
#' @format A list of signatures
NULL

#' Toy dataset to test the package
#' 
#' A downsampled version (300 cells) of the single-cell dataset by Zilionis et
#' al. (2019) <doi:10.1016/j.immuni.2019.03.009>, with precalculated PCA and
#' UMAP reductions.
#' 
#' @name query.seurat
#' @docType data
#' @format A Seurat object
NULL
