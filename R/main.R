#' Filter single-cell data by cell type
#'
#' Apply scGate to filter specific cell types in a query dataset
#'
#' @param data Seurat object containing a query data set - filtering will be applied to this object
#' @param gating.model A tabular model with scGate signatures. See Details for this format
#' @param pca.dim Number of dimensions for cluster analysis
#' @param assay Seurat assay to use
#' @param pos.thr Minimum UCell score value for positive signatures
#' @param neg.thr Maximum UCell score value for negative signatures
#' @param nfeatures Number of variable genes for dimensionality reduction
#' @param k.param Number of nearest neighbors for knn smoothing
#' @param pca.dim Number of principal components for dimensionality reduction
#' @param resol Resolution for cluster analysis (if \code{by.knn=FALSE})
#' @param ncores Number of processors for parallel processing (requires \code{\link{future.apply}})
#' @param output.col.name Column name with 'pure/impure' annotation
#' @param min.cells Stop iterating if fewer than this number of cells is left
#' @param additional.signatures A list of additional signatures, not included in the model, to be evaluated (e.g. a cycling signature). The scores for this
#'     list of signatures will be returned but not used for filtering.
#' @param by.knn Perform k-nearest neighbor smoothing for UCell scores
#' @param keep.ranks Store UCell rankings in Seurat object. This will speed up calculations if the same object is applied again with new signatures.
#' @param genes.blacklist Genes blacklisted from variable features. The default loads the list of genes in \code{scGate::genes.blacklist.default};
#'     you may deactivate blacklisting by setting \code{genes.blacklist=NULL}
#' @param seed Integer seed for random number generator
#' @param verbose Verbose output

#' @return A new metadata column \code{is.pure} is added to the query Seurat object, indicating which cells correspond to the desidered cell type(s).
#'     The \code{active.ident} is also set to this variable.
#' @details Models for scGate are dataframes where each line is a signature for a given filtering level. A database of models can be downloaded using the function
#'     \code{get_scGateDB()}. You may directly use the models from the database, or edit one of these models to generate your own custom gating model.   
#' @examples
#' library(scGate)
#' models <- get_scGateDB()
#' query <- scGate(query, model=models$human$generic$PanBcell)
#' DimPlot(query)
#' query.filtered <- subset(query, subset=is.pure=="Pure")
#' 
#' ############
#' # Using a custom model
#' 
#' my.model <- load_scGate_model("custom_model.tsv")
#' plot_tree(my.model)
#' query <- scGate(query, model=my.model)
#' query.filtered <- subset(query, subset=is.pure=="Pure")
#' 
#' @seealso \code{\link{load_scGate_model}} \code{\link{get_scGateDB}} \code{\link{plot_tree}} 
#' @export

scGate <- function(data, model, pos.thr=0.2, neg.thr=0.2, assay="RNA", ncores=1, seed=123, keep.ranks=FALSE,
                   min.cells=30, nfeatures=2000, pca.dim=30, resol=3, output.col.name = 'is.pure',
                   by.knn = TRUE, k.param=10, genes.blacklist="default", additional.signatures=NULL, verbose=TRUE) {
  
  set.seed(seed)
  #check gene blacklist
  if (!is.null(genes.blacklist)) {
    if (length(genes.blacklist)==1 && genes.blacklist == "default") {  #Default
       genes.blacklist = genes.blacklist.default
    }  
    if (is.list(genes.blacklist)) {
      genes.blacklist <- unlist(genes.blacklist)
    }
    genes.blacklist <- unique(genes.blacklist)
  }
  
  #if model is NULL, use a default model (ADD)
  
  list.model <- table.to.model(model)
  # compute signature scores using UCell
  if (verbose) {
    message("Computing UCell scores for signatures...\n")
  }
  data <- score.computing.for.scGate(data, model, assay, ncores=ncores, keep.ranks=keep.ranks)
  
  #Iterate over levels
  q <- data  #local copy to progressively remove cells
  tot.cells <- dim(q)[2]
  
  ## prepare output object (one pure/impure flag by level)
  output_by_level <- rep("Impure",length(list.model)*tot.cells)
  dim(output_by_level) <- c(tot.cells,length(list.model))
  colnames(output_by_level) <- names(list.model)
  output_by_level <- data.frame(output_by_level)
  rownames(output_by_level) <- data@meta.data%>%rownames()
  
  for (lev in 1:length(list.model)) {
    if (verbose) {
      message(sprintf("Running scGate on level %i...", lev))
    }
    #We can introduce iterations over same level until convergence (optional)
    pos.names <- sprintf("%s_UCell", names(list.model[[lev]]$positive))
    neg.names <- sprintf("%s_UCell", names(list.model[[lev]]$negative))
    
    ##Reduce parameter complexity at each iteration
    pca.use <- round((3/4)**(lev-1) * pca.dim)
    nfeat.use <- round((3/4)**(lev-1) * nfeatures)
    res.use <- round((3/4)**(lev-1) * resol)
    
    q <- find.nn(q, by.knn=by.knn, assay=assay, min.cells=min.cells, nfeatures=nfeat.use, npca=pca.use, 
                 k.param=k.param, genes.blacklist=genes.blacklist)
    
    if(!by.knn){
      q  <- FindClusters(q, resolution = res.use, verbose = FALSE)
      q$clusterCT <- q@active.ident
    }
    
    q <- filter_bymean(q, positive=pos.names, negative=neg.names, pos.thr=pos.thr, assay=assay,
                       min.cells=min.cells, neg.thr=neg.thr, by.knn = by.knn)
    
    n_rem <- sum(q$is.pure=="Impure")
    frac.to.rem <- n_rem/tot.cells
    mess <- sprintf("scGate: Detected %i non-pure cells at level %i", n_rem, lev)
    if (verbose) { message(mess) }
    
    ## retain pure cells will give us info in case of all cell where filtered
    retain_pure_cells <- q$is.pure=="Pure"
    
    if(any(retain_pure_cells)){
      output_by_level[names(retain_pure_cells[retain_pure_cells==T]),lev] <- "Pure"  # save layer output
      q <- subset(q, is.pure=="Pure")
    }else{
      break  # in case of all cells became filtered, we do not continue with the next layer
    }
  }
  if(any(retain_pure_cells)){  
    pure.cells <- colnames(q)
    data <- AddMetaData(data,col.name = output.col.name,metadata = rep("Impure",tot.cells))
    data@meta.data[pure.cells, output.col.name] <- "Pure"
    
    # save output by level
    for(name.lev in names(list.model)){
      data <- AddMetaData(data,col.name = paste0(output.col.name,".",name.lev),metadata = output_by_level[[name.lev]])
    }
  }else{
    message(sprintf("Warning, all cells were removed at level %i. Consider reviewing signatures or model layout...", lev))
    data <- AddMetaData(data,col.name = output.col.name,metadata = rep("Impure",tot.cells))
    # save output by level
    for(name.lev in names(list.model)){
      data <- AddMetaData(data,col.name = paste0(output.col.name,".",name.lev),metadata = output_by_level[[name.lev]])
    }
  }
  
  Idents(data) <- output.col.name
  
  n_rem <- sum(data[[output.col.name]]=="Impure")
  frac.to.rem <- n_rem/tot.cells
  mess <- sprintf("\n### Detected a total of %i non-pure cells for selected signatures - %.2f%% cells marked for removal (active.ident)",
                  n_rem, 100*frac.to.rem)
  message(mess)  #verbose?
  
  return(data)
}


#' Plot model tree
#'
#' View scGate model as a decision tree
#'
#' @param model A scGate model to be visualized
#' @param box.size Box size
#' @param edge.text.size Edge text size
#' @return A plot of the model as a decision tree. At each level, green boxes indicate the 'positive' (accepted) cell types, red boxed indicate the
#'     'negative' cell types (filtered out). The final Pure population is the bottom right subset in the tree.
#' @examples
#' plot_tree(model)
#' @export


plot_tree <- function(model, box.size = 8, edge.text.size = 4) {
  
  require(ggparty)
  nlev <- length(unique(model$levels))
  
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
  
  pn[[nlev]] <- partynode(nlev+1, split = partysplit(nlev, index=1:2, breaks = 0),
                          kids = list(partynode(nlev+2),
                                      partynode(nlev+3))) 
  
  for (i in (nlev-1):1) {
    pn[[i]] <- partynode(i, split = partysplit(i, index=1:2, breaks=0),
                         kids = list(partynode(i+1),
                                     pn[[i+1]]))
  }
  
  #first element in list has complete structure
  py <- party(pn[[1]], df)
  
  
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
  
  gg <- ggparty(py)
  gg$data$info <- sign.annot
  gg$data$p.value <- is.pos
  
  
  gg$data$breaks_label[grep("<=", gg$data$breaks_label)] <- "Negative"
  gg$data$breaks_label[grep(">", gg$data$breaks_label)] <- "Positive"
  
  gg <- gg + geom_edge() +
    geom_edge_label(size = edge.text.size) +
    geom_node_label(ids = "inner",
                    mapping = aes(col = p.value),
                    line_list = list(aes(label=info)),
                    line_gpar = list(list(size = box.size)))  +
    geom_node_label(ids = "terminal",
                    mapping = aes(col = p.value),
                    nudge_y=0.01,
                    line_list = list(aes(label=info)),
                    line_gpar = list(list(size = box.size))) +
    scale_color_manual(values=c("#f60a0a", "#00ae60")) +
    theme(legend.position = "none", plot.margin = unit(c(1,1,1,1), "cm")) 
  
  return(gg)
}


#' Load a single scGate model
#'
#' Loads a custom scGate model into R. For the format of these models, have a look or edit one of the default models obtained with \code{\link{get_scGateDB()}}
#'
#' @param model_file scGate model file, in .tsv format.
#' @param master.table File name of the master table (in repo_path folder) that contains cell type signatures.
#' @examples
#' library(scGate)
#' my.model <- load_scGate_model("my_model.tsv")
#' plot_tree(my.model)
#' query <- scGate(query, model=my.model)
#' @seealso \code{\link{scGate}} \code{\link{get_scGateDB}} 
#' @export

load_scGate_model <- function(model_file, master.table = "master_table.tsv"){
  
  model <- read.table(model_file, sep ="\t",header =T)
  model <- use_master_table(model, master.table = master.table)
  
  return(model)
}


#' Load scGate model database
#'
#' Download, update or load local version of the scGate model database. These are stored in a GitHub repository, from where you can download specific 
#' versions of the database.
#' 
#' @param destination destination path for db storage. The default is current location. 
#' @param force_update  Whether to update an existing database. WARNING: note that setting it TRUE, the current models folder will be overwritten . 
#' @param version Specify the version of the scGate_models database (e.g. 'v0.1'). By default downloads the latest available version.
#' @param repo_url  URL path to scGate model repository database
#' @return A list of models, organized according to the folder structure of the database. See the examples below.
#' @details Models for scGate are dataframes where each line is a signature for a given filtering level. A database of models can be downloaded using the function
#'     \code{get_scGateDB()}. You may directly use the models from the database, or edit one of these models to generate your own custom gating model.  
#' @examples 
#' library(scGate)
#' scGate.model.db <- get_scGateDB()
#' # To see a specific model, browse the list of models:
#' scGate.model.db$human$generic$Myeloid
#' # Apply scGate with this model
#' query <- scGate(query, model=scGate.model.db$human$generic$Myeloid)
#' @seealso \code{\link{scGate}} \code{\link{load_scGate_model}}  
#' @export

get_scGateDB <- function(destination = ".",
                         force_update = FALSE,
                         version = "latest",
                         repo_url = "https://github.com/carmonalab/scGate_models"){
  
  if (version=="latest") {
    repo_url_zip = sprintf("%s/archive/master.zip", repo_url)
    repo.name <- "scGate_models-master"
    repo.name.v <- repo.name
  } else {
    repo_url_zip = sprintf("%s/archive/refs/tags/%s.zip", repo_url, version)
    repo.name = sprintf("scGate_models-%s", version)
    repo.name.v <- sprintf("scGate_models-%s", gsub("^v","",version, perl = T))  #for some reason GitHub remove the 'v' from repo name after unzipping
  }
  
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
  allfiles <- list.files(repo_path.v, recursive = T, full.names = T)
  modelfiles <- grep("scGate_Model.tsv", allfiles, value = T)
  uniq_dirs <- sort(unique(dirname(modelfiles)))
  
  model_db <- list()
  for (dir in uniq_dirs) {
    pattern <- paste0(repo_path.v, '\\/?')
    tmp <- gsub(pattern, "", dir, perl=T)
    sub <- strsplit(tmp, split="/")[[1]]
    
    if (length(sub)==1) {
      model_db[[sub[1]]] <- load.model.helper(dir)
    } else if (length(sub)==2) {
      model_db[[sub[1]]][[sub[2]]] <- load.model.helper(dir)
    } else if (length(sub)==2) {
      model_db[[sub[1]]][[sub[2]]][[sub[[3]]]] <- load.model.helper(dir)
    } else {
      message(sprintf("Warning: max depth for scGate models is 3. Skipping folder %s", dir))
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
#' @param return_contingency  Logical indicating if contingency table must be returned. Default is FALSE
#' @examples
#' results <- performance.metrics(actual= sample(c(1,0),20,replace =T),pred =  sample(c(1,0),20,replace = T,prob = c(0.65,0.35) ) )
#' @export

performance.metrics <- function(actual,pred,return_contingency =F){
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
    nam.act <- unique(actual)%>%sort(decreasing = T)%>%as.character()  # 
    nam.pred <- unique(pred)%>%sort(decreasing = T)%>%as.character()
    ct <- ct[nam.act,nam.pred]  
    return(list('conting' = ct,
                'summary' = res.Summary ))
  }
  
}

