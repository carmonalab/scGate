run_scGate_singlemodel <- function(data, model, pos.thr=0.2, neg.thr=0.2, assay=NULL, slot="data",
                                   reduction="calculate", nfeatures=2000, pca.dim=30,
                                   param_decay=0.25, min.cells=30, k.param=30,
                                   smooth.decay=0.1, smooth.up.only=FALSE,
                                   genes.blacklist="default", verbose=FALSE,
                                   colname="is.pure", save.levels=FALSE,
                                   return.as.meta=TRUE) {
  
  if (!inherits(model, "data.frame")) {
    stop("Invalid scGate model. Please check the format of your model")
  }
  
  list.model <- table.to.model(model)
  
  q <- data  #local copy to progressively remove cells
  tot.cells <- ncol(q)
  meta.cols <- c()
  
  ## prepare output object (one pure/impure flag by level)
  output_by_level <- rep("Impure",length(list.model)*tot.cells)
  dim(output_by_level) <- c(tot.cells,length(list.model))
  colnames(output_by_level) <- names(list.model)
  output_by_level <- data.frame(output_by_level)
  rownames(output_by_level) <- rownames(q@meta.data)
  
  for (lev in 1:length(list.model)) {
    
    if (ncol(q) < 2)  break  # if at any level we reach a number of cells below this threshold,
                                # we skip computation, considering 'Impure' by default

    if (verbose) {
      message(sprintf("Running scGate on level %i...", lev))
    }

    pos.names <- sprintf("%s_UCell", names(list.model[[lev]]$positive))
    neg.names <- sprintf("%s_UCell", names(list.model[[lev]]$negative))
    all.names <- c(pos.names, neg.names)
    
    ##Reduce parameter complexity at each iteration
    if (param_decay < 0 | param_decay > 1) {
      stop("Parameter param_decay must be a number between 0 and 1")
    }
    
    k.use <- round((1-param_decay)**(lev-1) * k.param)
    
    smooth.decay.use <- 1-(1-smooth.decay)*(1-param_decay*(lev-1))
    
    if (reduction=="calculate") {
      pca.use <- round((1-param_decay)**(lev-1) * pca.dim)
      nfeat.use <- round((1-param_decay)**(lev-1) * nfeatures)
    } else {
      pca.use <- pca.dim
    }
    q <- find.nn(q, assay=assay, slot=slot, signatures=all.names,min.cells=min.cells,
                 nfeatures=nfeat.use, reduction=reduction, npca=pca.use, k.param=k.use,
                 smooth.decay=smooth.decay.use, smooth.up.only=smooth.up.only,
                 genes.blacklist=genes.blacklist)
    
    q <- filter_bymean(q, positive=pos.names, negative=neg.names, assay=assay,
                       pos.thr=pos.thr, neg.thr=neg.thr)
    
    Idents(q) <- "is.pure"
    n_rem <- sum(Idents(q)=="Impure")
    frac.to.rem <- n_rem/tot.cells
    mess <- sprintf("scGate: Detected %i non-pure cells at level %i", n_rem, lev)
    if (verbose) { message(mess) }
    
    ## How many cells passed the filter
    pure.cells <- colnames(q)[Idents(q)=="Pure"]
    
    if(length(pure.cells)>0){
      # save layer output
      output_by_level[pure.cells, lev] <- "Pure"
      q <- subset(q, idents="Pure")
    } else {
      break  # in case of all cells became filtered, we do not continue with the next layer
    }
  }
  
  #Add 'pure' labels to metadata
  data <- AddMetaData(data,col.name = colname, metadata = rep("Impure",tot.cells))
  meta.cols <- c(meta.cols, colname)
  
  if(length(pure.cells)>0){
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
      meta.cols <- c(meta.cols, combname)
    }
  }
  if (return.as.meta) {
    return(data@meta.data[,meta.cols, drop=F])
  } else {
    return(data)
  }
}


find.nn <- function(q, assay = "RNA", slot="data", signatures=NULL, npca=30,
                    nfeatures=2000, k.param=10,
                    smooth.decay=0.1, smooth.up.only=FALSE,
                    min.cells=30, reduction="calculate", genes.blacklist=NULL) {
  
  DefaultAssay(q) <- assay
  ncells <- length(Cells(q))
  ngenes <- nrow(q)
  
  notfound <- signatures[!signatures %in% colnames(q[[]])]
  signatures <- signatures[signatures %in% colnames(q[[]])]
  
  if (length(notfound)>0) {
    message(paste0("Warning: signatures not found: ", notfound))
  }
  
  if(ncells < min.cells){  #Do not do knn-smoothing
    return(q)
  }  
  
  if (reduction=="calculate") {
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
      VariableFeatures(q) <- rownames(q)
    }
    
    VariableFeatures(q) <- setdiff(VariableFeatures(q), genes.blacklist)
    
    q <- ScaleData(q, verbose=FALSE)
    q <- suppressWarnings(RunPCA(q, features = VariableFeatures(q),
                                 npcs=npca, verbose = FALSE,
                                 reduction.key = "knnPCA_"))
    red.use <- 'pca'
  } else {
    red.use <- reduction
  }
  if (smooth.decay>1) {
    smooth.decay=1
  }
  
  #Smooth scores by kNN neighbors
  q <- SmoothKNN(q, signature.names=signatures, reduction=red.use, k=k.param,
                 decay=smooth.decay, up.only=smooth.up.only, suffix = NULL)
  
  return(q)
  
}

## Filter by mean
filter_bymean <- function(q, positive, negative, pos.thr=0.1, neg.thr=0.2, assay="RNA") {
  
  DefaultAssay(q) <- assay
  ncells <- ncol(q)
  
  positive <- positive[positive %in% colnames(q[[]])]
  negative <- negative[negative %in% colnames(q[[]])]
  
  cols <- c(positive, negative)
  means <- list()
  
  scores <- q[[]][,cols, drop=FALSE]

  if(length(positive)>1){
    pos <- scores[,positive] |> apply(1,max)
  }else{
    pos <- scores[,positive]
  }
  if(length(negative)>1){
    neg <- scores[,negative] |> apply(1,max)
  }else{
    neg <- scores[,negative]
  }
  ispure <- rep("Impure", ncol(q))
  
  if(length(positive)>0 & length(negative)>0) {
    ispure[pos>pos.thr & neg<neg.thr] <- "Pure"
  } else if (length(positive)>0) {
    ispure[pos>pos.thr] <- "Pure"
  } else if (length(negative)>0) {
    ispure[neg<neg.thr] <- "Pure"
  } else {
    stop("No valid signatures were provided.")
  }
      
  q@meta.data[,"is.pure"] <- ispure
  
  return(q)
}

score.computing.for.scGate <- function(data, model, bpp=SerialParam(), assay="RNA", slot="data",
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
  
  data <- AddModuleScore_UCell(data, features = signatures, assay=assay, slot=slot,
                                      BPPARAM = bpp, storeRanks = keep.ranks, maxRank = maxRank)
  
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


load.model.helper <- function(models_path, master.table = "master_table.tsv",  verbose=verbose) {

  df.models.toimpute <- list()
  files.to.impute <- list.files(file.path(models_path),"_scGate_Model.tsv")
  if(length(files.to.impute)==0){
    stop("Please, provide some model table files in your 'model folder' or set models_path = NULL for using the default ones")
  }
  # load models to impute
  for(f in files.to.impute){
    model.name <- strsplit(f,"_scGate_Model.tsv")[[1]][1]
    if(verbose) message(paste0("loading ",model.name))
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

#Serial or parallel lapply
my.lapply <- function(X=NULL, FUN, ncores=1, BPPARAM) {
  if (ncores>1) {
    BiocParallel::bplapply(X=X, FUN=FUN, BPPARAM = BPPARAM)
  } else {
    lapply(X=X, FUN=FUN)
  }
}
