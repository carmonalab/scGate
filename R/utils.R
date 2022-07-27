run_scGate_singlemodel <- function(data, model, pos.thr=0.2, neg.thr=0.2, assay=NULL, slot="data",
                                   reduction="calculate", nfeatures=2000, pca.dim=30, resol=3,
                                   param_decay=0.25, min.cells=30, by.knn = TRUE, k.param=10, 
                                   genes.blacklist="default", verbose=FALSE,
                                   colname="is.pure", save.levels=FALSE) {
  
  if (!inherits(model, "data.frame")) {
    stop("Invalid scGate model. Please check the format of your model")
  }
  
  list.model <- table.to.model(model)
  
  q <- data  #local copy to progressively remove cells
  tot.cells <- ncol(q)
  
  ## prepare output object (one pure/impure flag by level)
  output_by_level <- rep("Impure",length(list.model)*tot.cells)
  dim(output_by_level) <- c(tot.cells,length(list.model))
  colnames(output_by_level) <- names(list.model)
  output_by_level <- data.frame(output_by_level)
  rownames(output_by_level) <- data@meta.data %>% rownames()
  
  for (lev in 1:length(list.model)) {
    if (verbose) {
      message(sprintf("Running scGate on level %i...", lev))
    }

    pos.names <- sprintf("%s_UCell", names(list.model[[lev]]$positive))
    neg.names <- sprintf("%s_UCell", names(list.model[[lev]]$negative))
    
    ##Reduce parameter complexity at each iteration
    if (param_decay < 0 | param_decay > 1) {
      stop("Parameter param_decay must be a number between 0 and 1")
    }
    
    if (reduction=="calculate") {
      pca.use <- round((1-param_decay)**(lev-1) * pca.dim)
      nfeat.use <- round((1-param_decay)**(lev-1) * nfeatures)
      res.use <- round((1-param_decay)**(lev-1) * resol)
    } else {
      pca.use <- pca.dim
    }
    q <- find.nn(q, by.knn=by.knn, assay=assay, slot=slot, min.cells=min.cells, nfeatures=nfeat.use, 
                 reduction=reduction, npca=pca.use, k.param=k.param, genes.blacklist=genes.blacklist)
    
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
      q <- subset(q, subset=`is.pure`=="Pure")
    } else {
      break  # in case of all cells became filtered, we do not continue with the next layer
    }
  }
  
  #Add 'pure' labels to metadata
  data <- AddMetaData(data,col.name = colname, metadata = rep("Impure",tot.cells))
  
  if(any(retain_pure_cells)){  
    pure.cells <- colnames(q)
    data@meta.data[pure.cells, colname] <- "Pure"
  } else {
    message(sprintf("Warning, all cells were removed at level %i. Consider reviewing signatures or model layout...", lev))
  }
  
  data@meta.data[,colname] <- factor(data@meta.data[,colname], levels=c("Pure","Impure"))
  
  # Save output by level
  if (save.levels) {
    for(name.lev in names(list.model)){
      combname <- paste0(colname,".",name.lev)
      data <- AddMetaData(data,col.name = combname, metadata = output_by_level[[name.lev]])
      data@meta.data[,combname] <- factor(data@meta.data[,combname], levels=c("Pure","Impure"))
    }
  }
  return(data)
}


find.nn <- function(q, assay = "RNA", slot="data", npca=30, nfeatures=2000, k.param=10,
                    min.cells=30, by.knn = F, reduction="calculate", genes.blacklist=NULL) {
  
  DefaultAssay(q) <- assay
  ncells <- length(Cells(q))
  ngenes <- nrow(q)
  
  if (reduction=="calculate") {
    if(ncells < min.cells){
      q$clusterCT <- 0    #with very few cells, consider them as a single cluster
      return(q)
    }  
    if (ngenes < nfeatures) {
      nfeatures <- ngenes
    }
    if (ngenes/2 < npca) {
      npca <- ngenes/2
    }
    
    if (slot=="counts") { 
      q <- NormalizeData(q, verbose = FALSE)
    }
    if (ngenes > 200) {  #only perform this for high-dim data
      q <- FindVariableFeatures(q, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
    } else {
      q@assays[[assay]]@var.features <- rownames(q)
    }
    
    q@assays[[assay]]@var.features <- setdiff(q@assays[[assay]]@var.features, genes.blacklist)
    
    q <- ScaleData(q, verbose=FALSE)
    q <- RunPCA(q, features = q@assays[[assay]]@var.features, npcs=npca, verbose = FALSE, reduction.key = "knnPCA_")
    
    red.use <- 'pca'
  } else {
    red.use <- reduction
  }
  
  q <- suppressMessages(FindNeighbors(q, reduction = red.use, dims = 1:npca, k.param = k.param, verbose=FALSE,
                                      return.neighbor = by.knn))
  return(q)
  
}

## Filter by mean
filter_bymean <- function(q, positive, negative, pos.thr=0.1, neg.thr=0.2,  min.cells=30,
                          assay="RNA", return_object = T, by.knn = F) {
  
  DefaultAssay(q) <- assay
  ncells <- dim(q)[2]
  notfound <- c(positive[!positive %in% colnames(q@meta.data)], negative[!negative %in% colnames(q@meta.data)])
  
  if (length(notfound)>0) {
    message(paste0("Warning: signatures not found: ", notfound))
  }
  positive <- positive[positive %in% colnames(q@meta.data)]
  negative <- negative[negative %in% colnames(q@meta.data)]
  
  cols <- c(positive, negative)
  means <- list()
  
  if(!by.knn){
    for (col in cols) {
      means[[col]] <- sapply(levels(q$clusterCT), function(x) {
        mean(q@meta.data[q$clusterCT == x, col])
      })
    }
    meds <- Reduce(rbind, means)
    if(is.null(dim(meds))){
      dim(meds) <- c(1,length(meds))
    }
    rownames(meds) <- cols
    
    pos <- vector(length=dim(meds)[2])
    neg <- vector(length=dim(meds)[2])
    for (j in 1:dim(meds)[2]) {
      pos[j] <- max(meds[positive,j])
      neg[j] <- max(meds[negative,j])
    }
    
    if(length(negative)>0){
      indices <- intersect(which(pos > pos.thr), which(neg < neg.thr))
    }else{
      indices <- which(pos > pos.thr) # case without negative signature
    }
    select.pures <- colnames(meds)[indices]
    ispure <- ifelse(q$clusterCT %in% select.pures,"Pure","Impure")
    
    
  } else{
    for (col in cols) {
      meta.nn <- sprintf("%s.nn", assay)
      if (ncells < min.cells) {   #very small dataset. Use all cells together
        neigs <- t(matrix(data = rep(1:ncells,ncells), nrow = ncells, ncol = ncells))
      } else {
        neigs <- q@neighbors[[meta.nn]]@nn.idx
      }
      m <- q[[col]][[1]][neigs]
      dim(m) <- dim(neigs)
      means[[col]] <- apply(m,1,mean)
    }
    meds <- Reduce(rbind, means)
    if(is.null(dim(meds))){
      dim(meds) <- c(1,length(meds))
    }
    rownames(meds) <- cols
    
    if(length(positive)>1){
      pos <- meds[positive,]%>%apply(2,max)
    }else{
      pos <- meds[positive,]
    }
    if(length(negative)>1){
      neg <- meds[negative,]%>%apply(2,max)
    }else{
      neg<- meds[negative,]
    }
    ispure <- rep("Impure",dim(q)[2])
    if(length(negative)>0){
      ispure[(pos > pos.thr)&(neg < neg.thr)] <- "Pure"
    }else{
      ispure[pos > pos.thr] <- "Pure"  # case without negative signature
    }
      
  }
  
  q$is.pure <- ispure
  if(return_object) {return(q)}
  
  return(ispure)
  
}

score.computing.for.scGate <- function(data, model, ncores=1, assay="RNA", slot="data",
                                       add.sign=NULL, keep.ranks=FALSE, maxRank=1500) {
  
  comb <- dplyr::bind_rows(model, .id = "Model_ID")
  # extract unique signatures
  model.uniq <- comb %>% dplyr::distinct(.data$name, .data$signature, .keep_all = T) 
  
  #Stop if there are signature with same name but different genes
  t <- table(model.uniq$name)
  dup <- names(t[t>1])
  if (length(dup)>0) {
    s <- paste(dup,collapse=", ")
    stop(sprintf("Different gene sets have been assigned to signature with same name: %s", s))
  }
  
  ## generate list object to be used in computing stage
  signatures <- model.uniq$signature %>% strsplit("[,; ]+") %>% lapply(unlist)   #also allow comma or space
  names(signatures) <- model.uniq$name
  
  if (!is.null(add.sign)) {
    if (!inherits(add.sign, "list")) {
      add.sign <- list("Additional_signature"=add.sign)
    }
    signatures <- append(signatures, add.sign)
  }
  
  data <- UCell::AddModuleScore_UCell(data, features = signatures, assay=assay, slot=slot,
                                      ncores=ncores, storeRanks = keep.ranks, maxRank = maxRank)
  
  return(data)
}

model.to.table <- function(scGate.model){
  tab <- data.frame(levels=NULL,use_as = NULL,name = NULL,signature = NULL)
  
  ### Define first column: Levels
  levels <- scGate.model%>%names()
  lev <- rep(levels,rep(2,length(levels)))
  
  len_pos_neg <- lapply(scGate.model,function(x){ 
    res = lapply(x,function(y){
      length(y)
    })
    return(res)
  })%>%unlist()
  
  extended.levels <- rep(lev,len_pos_neg)
  
  # second column: "use_as"
  useas <- rep(c("positive","negative"),length(levels))
  useas <- rep(useas , len_pos_neg)
  
  #third column: name
  signature.names <- lapply(scGate.model,function(x){ 
    res = lapply(x,function(y){
      names(y)
    })
    return(res)
  })%>%unlist()
  
  ## Four column: signature
  signatures <- lapply(scGate.model,function(x){ 
    res = lapply(x,function(y){
      lapply(y, function(z){
        paste(z,collapse = ";")
      })
    })
    return(res)
  })%>%unlist()
  
  tab <- data.frame("levels"=extended.levels,"use_as" = useas, "name" = signature.names,"signature" = signatures)
  return(tab)
}


table.to.model <- function(scGate.table){
  mod <- list()
  for(i in 1:nrow(scGate.table)){ 
    lev <- scGate.table$levels[i] 
    useas <- tolower(scGate.table$use_as[i])
    if(!useas %in% c("positive","negative")){
      message(sprintf("Error: row %i do not contain neither, 'positive' or 'negative' strings in 'use_as' column",i))
      return(NULL)
    }
    
    sign <- scGate.table$signature[i]
    
    name <- scGate.table$name[i]
    mod[[lev]][[useas]][[name]] <- strsplit(sign, "[,; ]+") %>% unlist()
  }
  return(mod)
}


load.model.helper <- function(models_path, master.table = "master_table.tsv") {

  df.models.toimpute <- list()
  files.to.impute <- list.files(file.path(models_path),"_scGate_Model.tsv")
  if(length(files.to.impute)==0){
    stop("Please, provide some model table files in your 'model folder' or set models_path = NULL for using the default ones")
  }
  # load models to impute
  for(f in files.to.impute){
    model.name <- strsplit(f,"_scGate_Model.tsv")[[1]][1]
    df.models.toimpute[[model.name]] <- read.table(file.path(models_path,f),sep ="\t",header =T)
  }
  # signature imputing
  imputed.models <-  lapply(df.models.toimpute,function(df){
    use_master_table(df.model = df, master.table = file.path(models_path, master.table))
  })
  model.list <- imputed.models
      

  return(model.list)
}

## This function allows to complete signatures in a table based model by using the name signature and a provided master.table of signatures
# the master.table must be a two column data.frame with two columns : 1) name: contains the signature names and 
# 2)signature: this column contain the genes present in each signature (separated with a semicolon) 

use_master_table <- function(df.model, master.table, name = "name",descript = "signature"){
  
  ## First, check if descript field exists in input table
  if(!descript %in% colnames(df.model)){
    df.model[descript] <- ""
  }
  
  ## second, check if there is something to do (or exit)
  input.sign <- df.model[[descript]]
  input.names <- df.model[[name]]
  complete.from.master.table <- as.vector(is.null(input.sign)|input.sign == "" | is.na(input.sign))
  
  if(!any(complete.from.master.table)){
    message("nothing to do, the model signatures are already provided")
    return(df.model)
  }
  
  #Does the master table exist?
  if(!file.exists(master.table)){
    stop("master_table.tsv file must be present in your 'model folder' unless signatures are completely specified in the models")
  }
  master.table <- read.table(master.table, sep ="\t", header =T) 
  
  # sanity check:
  warn <- setdiff(input.names[complete.from.master.table],master.table[[name]])
  if(length(warn)>0){
    stop(sprintf("signatures '%s' are not present in the provided master.signature table",paste(warn,collapse = " ; " )))
  }
  
  merged <- merge(df.model[complete.from.master.table,],master.table,by = name,suffixes = c("","_from.master"),all.x = T)
  vect <- merged[[paste0(descript,"_from.master")]]
  names(vect) <- merged[[name]]
  if(merged%>%nrow > sum(complete.from.master.table)){
    stop("please, master.table must not contain duplicated signature names")
  }
  
  # replace ensuring correct order (notice that the merge can change row order)
  nms <- df.model[complete.from.master.table,name]
  df.model[complete.from.master.table,descript] <- vect[nms]
  #  df.model[descript] <- output.sign
  return(df.model)
}
