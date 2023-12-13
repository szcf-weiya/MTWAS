#' Data preparation
#' @param dat plink bfile
#' @param E original gene expression matrix (list), the length of the list is the number of tissues,  each element is a sample*gene matrix, rownames are sample id, colnames are gene id
#' @param E.info data.frame with columns: chr, start, end, pos
#' @return standard mtwas data format
#' @export
format_twas_dat <- function(dat, E, E.info){
  #### pre-check
  ntissue <- length(E)
  genelist_E <- c()
  sam_E <- sam_use <- c()
  for(k in 1:ntissue){
    genelist_E <- c(genelist_E, colnames(E[[k]]))
    sam_E <- c(sam_E,rownames(E[[k]]))
  }
  gene_use <- unique(intersect(genelist_E, E.info$gene))
  if(length(gene_use)==0){
    stop('There are no overlapped genes between E.info$gene and colnames of E. No data for training')
  }
  sam_use <- unique(intersect(dat$fam$V2, sam_E))
  if(length(sam_use)==0){
    stop('There are no overlapped samples between dat$fam$V2 and rownames of E. No data for training')
  }
  E0 <- E
  E.info <- E.info[match(gene_use,E.info$gene),]
  dat$bed <- dat$bed[match(sam_use,dat$fam$V2),]
  dat$fam <- dat$fam[match(sam_use,dat$fam$V2),]


  #### 0. imputation for each gene
  E.eachgene <- list()
  for(l in 1:length(gene_use)){
    E.eachgene[[l]] <- matrix(NA,nrow=length(sam_use),ncol=ntissue)
    for(k in 1:ntissue){
      E.eachgene[[l]][,k] <- E0[[k]][match(sam_use,rownames(E0[[k]])),
                                     which(colnames(E0[[k]])==gene_use[l])]

    }
  }
  temp <- list(E=E.eachgene)
  E.eachgene.imp <- (impute.E(temp)$imp)$E
  ###
  E <- list()
  for(k in 1:ntissue){
    E[[k]] <-  matrix(NA,nrow=length(sam_use),ncol=length(gene_use))
    for(l in 1:length(gene_use)){
      E[[k]][,l] <-  E.eachgene.imp[[l]][,k]
    }
    colnames(E[[k]]) <- gene_use
    rownames(E[[k]]) <- sam_use
  }
  names(E) <- names(E0)
  #### format
  exp1 <- E.info[,c('chr','start','end','gene')]
  colnames(exp1) <- c('V1','V2','V3','gene')
  exp1$V2 <- pmax(0,exp1$V2-1e6) ## cis SNPs
  exp1$V3 <- exp1$V3+1e6

  ind <- find.index(exp1,dat$bim,type='pos')

  pos <- lapply(1:length(gene_use),function(x){return(ind$df2[which(ind$df1==x)])}) ##length: # genes

  twas_data <- list()
  twas_data$E <- E
  twas_data$E.info <- E.info
  twas_data$dat <- dat
  twas_data$pos <- pos
  return(twas_data)
}


#' match index
#' @param df1 data.frame 1
#' @param df2 data.frame 2
#' @param type type
#' @return matched index
#' @importFrom GenomicRanges GRanges GNCList findOverlaps
#' @importFrom IRanges IRanges
find.index <- function(df1,df2,type='reg'){
  #colnames(df1) <- colnames(df2) <- c('V1','V2','V3')
  df1.gr = GRanges (IRanges(start = df1$V2, end = df1$V3), seqnames=df1$V1)
  if(type=='reg'){
    df2.gr = GRanges(IRanges(start=df2$V2, end = df2$V3), seqnames = df2$V1)
  }
  if(type=='pos'){
    df2.gr = GRanges(IRanges(start=df2$V4, end = df2$V4), seqnames = df2$V1)
  }
  df1.git  = GNCList(df1.gr)
  df2.git  = GNCList(df2.gr)
  overlap_git = findOverlaps(df2.git, df1.git)
  overlap_git
  temp <- as.data.frame(overlap_git)
  colnames(temp) <- c('df2','df1')
  return(temp)
}



#' Impute the list of expression matrix `E`
#' @param x twas data list
#' @return imputed twas data list and the imputation MSE error
#' @importFrom missForest missForest
#' @export
impute.E = function(x) {
  ntissue = length(x$E)
  mse = vector("list", length = ntissue)
  for (j in 1:ntissue) {
    nomiss = which(apply(x$E[[j]], 2, function(z) {any(!is.na(z))}))
    mat = x$E[[j]][, nomiss]
    colnames(mat) = NULL
    res = missForest(mat, variablewise = T)
    x$E[[j]][, nomiss] = res$ximp
    msej = numeric(length(nomiss))
    msej[nomiss] = res$OOBerror
    mse[[j]] = msej
  }
  list(imp = x, mse.error = mse)
}
