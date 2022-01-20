find.nn <- function(q, assay = "RNA", npca=30, nfeatures=2000, k.param=10, min.cells=30, by.knn = F,
                    genes.blacklist=NULL) {
  
  DefaultAssay(q) <- assay
  ncells <- length(Cells(q))
  
  if(ncells < min.cells){
    q$clusterCT <- 0    #with very few cells, consider them as a single cluster
    return(q)
  }  
  
  q <- NormalizeData(q, verbose = FALSE)
  q <- FindVariableFeatures(q, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
  
  q@assays[[assay]]@var.features <- setdiff(q@assays[[assay]]@var.features, genes.blacklist)
  
  q <- ScaleData(q, verbose=FALSE)
  q <- RunPCA(q, features = q@assays[[assay]]@var.features, npcs=npca, verbose = FALSE)
  q <- suppressMessages(FindNeighbors(q, reduction = "pca", dims = 1:npca, k.param = k.param, verbose=FALSE,
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

score.computing.for.scGate <- function(data, model, ncores=1, assay="RNA", 
                                       add.sign=NULL, keep.ranks=FALSE, maxRank=1500) {
  
  # extract unique signatures
  model.uniq <- model %>%distinct(name, .keep_all = T) 
  
  ## generate list object to be used in computing stage
  signatures <- model.uniq$signature %>% strsplit("[,; ]+") %>% lapply(unlist)   #also allow comma or space
  names(signatures) <- model.uniq$name
  
  if (is.list(add.sign) & length(add.sign)>0) {
     signatures <- append(signatures, add.sign)
  } 
  
  data <- UCell::AddModuleScore_UCell(data, features = signatures, assay=assay, 
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
