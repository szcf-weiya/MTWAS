#' @importFrom stats cor
get.cor <- function(x1,x2){
  nona <- which((!is.na(x1))&(!is.na(x2)))
  if(length(nona)>1){
    return(cor(x1[nona],x2[nona]))
  }else{
    return(NA)
  }
}


#' Create matrix from select terms
#' @param xx raw matrix with colnames with natural indexes
#' @param terms selected terms. If exist interact terms like \code{i*j}, then
#'  the i-th column would multiply j-th column, and result in a new column.
#' @return a matrix
#' @export
create_pmatrix_from_terms = function(xx, terms)
{
  nt = length(terms);
  nr = nrow(xx);
  pmatrix = matrix(0, nr, 0);

  if (nt > 0)
    for(it in 1:nt)
    {
      term = terms[it];
      if (grepl("*",term,fixed=T))
      {
        splits = strsplit(term,"*",fixed=T)[[1]];
        id1 = as.numeric(splits[1]);
        id2 = as.numeric(splits[2]);
        pmatrix = cbind(pmatrix, term=xx[,id1]*xx[,id2]);
      }
      else
      {
        id  = as.numeric(term);
        pmatrix = cbind(pmatrix, term=xx[,id]);
      }
    }
  return(pmatrix);
}

#' Detect the cross-tissue eQTLs
#' @param x TWAS data, which is a list of three elements: \describe{
#'  \item{E}{list of cell types, each element is a matrix of size n_sample x n_gene}
#'  \item{dat}{a list with three named elements. \describe{
#'    \item{bed}{a matrix of size n_sample x n_snp}
#'    \item{bim}{a data frame with n_snp rows, each row is the info for the snp}}
#'  }
#'  \item{pos}{a list of genes, each element is a vector containing the SNP index that
#'  locate in the gene, and it corresponds to each row of \code{E.info}.
#'  It can be derived from \code{E.info} and \code{dat}.}
#' }
#' @param nne number of genes to be considered. Default is -1, which means to use all genes.
#' @param gam the tuning parameter \eqn{\gamma} in the EBIC criterion. Default is 1.
#' @param nF0 the smallest number of main effects. Default is 5.
#' @param ncores if using parallel computing, the number of cores to be used. Default is 1.
#' @param verbose whether to print running message. Default is \code{TRUE}.
#' @param allow_empty whether to allow empty set. Default is \code{TRUE}.
#' @return a list of size n_gene, each element is the shared eQTL for the gene.
#' @importFrom parallel mclapply
#' @export
select.ct.eQTL = function(x, nne = -1, gam = 1, nF0 = 5,
                          verbose = T, allow_empty = T,
                          ncores = 1L) {
  # number of genes
  if (nne == -1)
    nne = length(x$pos)
  if (ncores == 1) {
    res = lapply(1:nne, function(j)
            soda_model_imp_share(j, gam = gam, x$E, x$dat, x$pos,
                                 nF0 = nF0, verbose = verbose, allow_empty = allow_empty))
  } else {
    res = mclapply(1:nne, soda_model_imp_share, gam = gam, nF0 = nF0,
                   E = x$E, dat = x$dat, pos = x$pos,
                   verbose = verbose, allow_empty = allow_empty,
                   mc.cores = ncores, mc.preschedule = F)
  }
  names(res)=x$E.info$gene
  res
}

#' Select cross-tissue via SODA
#' @param j the index of gene
#' @param gam TODO
#' @param E list of cell types, each element is a matrix of size n_sample x n_gene
#' @param dat a list with three named elements. \describe{
#' \item{bed}{a matrix of size n_sample x n_snp}
#' \item{bim}{a data frame with n_snp rows, each row is the info for the snp}}
#' @param pos a list of genes, each element is a vector containing the SNP index that
#' locate in the gene, and it corresponds to each row of \code{E.info}.
#' It can be derived from \code{E.info} and \code{dat}.
#' @param nF0 TODO
#' @param verbose whether to print running messages
#' @param allow_empty whether to allow empty variable set
#' @importFrom stats prcomp
#' @export
soda_model_imp_share <- function(j, gam, E, dat, pos, nF0 = 5,
                                 verbose = TRUE,
                                 allow_empty = TRUE){
  if (length(pos[[j]]) > 0){
    xx <- (as.matrix(dat$bed[,pos[[j]]]))
    #yy <- E[,j]
    ### to remove the SNP with maf<0.05 in tissue-specific dataset
    ind0 <- 1:ncol(xx)
    train.ind <- 1:nrow(xx)
    if(length(ind0)>1){
      # dat <- list(xx=xx,yy=yy)
      # load(paste0(paste0('./dat/chr',chr,'/dat_',j,'.RData')))
      # if(scale){
      #   train.x <- scale(xx[train.ind,])
      #   train.y <- scale(yy[train.ind])
      #   test.x <- scale(xx[test.tissue.ind,])
      #   test.y <- scale(yy[test.tissue.ind])
      # }else{
      #train.x <- (xx[train.ind,])
      train.x <- t(t(xx[train.ind,])-apply(xx[train.ind,],2,mean))
      #train.y <- (yy[train.ind])
      #test.x <- (xx[test.ind,])
      #test.y <- (yy[test.tissue.ind])
      #}
      train.x[is.na(train.x)] <- 0
      train.x[is.nan(train.x)] <- 0

      kct.all <- 1:length(E) # cell type
      if(T){
        for(l in 1:length(kct.all)){
          if(l==1){
            yy0 <- E[[kct.all[l]]][train.ind,j]
          }else{
            yy0 <- cbind(yy0,E[[kct.all[l]]][train.ind,j])
          }
        }
      }
      if(length(which(is.na(yy0[1,])))>0){
        yy0 <- yy0[,-which(is.na(yy0[1,]))]
      }
      temp <- prcomp(yy0)
      yy.pca <- temp$x
      PoV <- temp$sdev^2/sum(temp$sdev^2)
      pca.var <- 0
      i <- 1
      fixset <- numeric(0)

      while(i<6){
        pca.var <- pca.var+PoV[i]
        yy <- yy.pca[,i]
        train.y <- (yy)-mean(yy)
        #test.x[is.na(test.x)] <- 0
        #test.x[is.nan(test.x)] <- 0
        model <- soda_main_allback(train.x,train.y, colnames(train.x),
                                   fixset = fixset, FALSE, verbose, gam, nF0, allow_empty)
        fixset <- unique(c(fixset,model$`Select Set`))
        i <- i+1
      }
      # fn <- paste0('')
      # fn.share <- paste0('/home/songs/gtex/joint_tissue_all/test_hm3/all_tissue/res_soda_linear/chr',
      #                    chr,'/soda_linear_cv',cv,'_gam',gam0,'_nF0',nF0,'_share.RData')
      #
      #fn.share <- paste0(path,'/soda_linear_cv',cv,'_gam',gam0,'_nF0',nF0,'_share.RData')
      return(fixset)
    }
  }
}

#' Select tissue-specific eQTLs
#' @param x data
#' @param ct.eQTL selected cross-tissue eQTLs from \code{select.ct.eQTL}
#' @param id.tissue the index of tissue (cell type) to be considered. Default is -1,
#' which returns the ts.eQTL for all tissues.
#' @param ncores number of cores. If larger than 1, then use parallel computing.
#' @inheritDotParams soda_model_imp_single gam scale nF0
#' @export
select.ts.eQTL = function(x, ct.eQTL, id.tissue = -1, ncores = 1, ...) {
  nne = length(ct.eQTL)
  if (id.tissue == -1) {
    ntissue = length(x$E)
    return(lapply(1:ntissue, function(i) select.ts.eQTL(x, ct.eQTL, i, ncores = ncores, ...)))
  }
  if (ncores == 1) {
    soda.imp = lapply(1:nne, function(j)
      soda_model_imp_single(j, fix = ct.eQTL, x$E[[id.tissue]], x$dat, x$pos, ...))
  } else {
    soda.imp = mclapply(1:nne, soda_model_imp_single,
                        E = x$E[[id.tissue]], dat = x$dat, pos = x$pos,
                        fix = ct.eQTL, ...,
                        mc.cores=ncores, mc.preschedule =F)
  }
  names(soda.imp)=x$E.info$gene
  soda.imp
}

#'
#' Select tissue-specific via SODA
#' @inheritParams soda_model_imp_share
#' @param fix a list of set containing the fixed effect for the gene
#' @param E here E is just for specific cell type, it is a matrix instead of a list
#' @param gam tuning parameter \eqn{\gamma} in EBIC criteria. Default is 1.
#' @param scale whether to scale the expression and the genetic data when regressing the expression on the genetic data. Default is \code{FALSE}.
#' @param nF0 minimum size of the set
#' @importFrom stats coef
#' @export
soda_model_imp_single <- function(j, fix, E, dat, pos, gam = 1, scale = F, nF0 = 5){
  fixset <- fix[[j]]
  if(length(pos[[j]])>0){
    xx <- (as.matrix(dat$bed[,pos[[j]]]))
    rsid00 <- dat$bim$V2[pos[[j]]]
    yy <- E[,j]
    ### to remove the SNP with maf<0.05 in tissue-specific dataset
    ind0 <- 1:ncol(xx)
    train.ind <- 1:nrow(xx)
    ### to remove the SNP with maf<0.05 in tissue-specific dataset
    if(length(ind0)>0){
      # dat <- list(xx=xx,yy=yy)
      # load(paste0(paste0('./dat/chr',chr,'/dat_',j,'.RData')))
      if(scale){
        train.x <- scale(xx[train.ind,])
        train.y <- scale(yy[train.ind])
        # test.x <- scale(xx[test.ind,])
        # test.y <- scale(yy[test.ind])
      }else{
        train.x <- (xx[train.ind,])
        train.y <- (yy[train.ind])
        # test.x <- (xx[test.ind,])
        # test.y <- (yy[test.ind])
      }
      if(all(is.na(train.y))){
        return(NA)
      }else{
        train.x[is.na(train.x)] <- 0
        #test.x[is.na(test.x)] <- 0
        train.x[is.nan(train.x)] <- 0
        #test.x[is.nan(test.x)] <- 0
        if(ncol(xx)==1){
          res <- list()
          res$best_term=1
          res$coef <- coef(stats::lm(train.y~train.x))
          res$single <- setdiff(1,fix[[j]])
          res$single.snp <- rsid00[res$single]
          res$common.snp <- rsid00[fix[[j]]]
          # res$corr <- get.cor(test.y,test.x)
          return(res)
        }else{
          single <- T
          iter <- 0
          allowempty <- T
          nfix <- length(fixset)
          while(length(fixset)==0 || single){
            #while(length(fixset)<=nfix || single){
            iter <- iter+1
            single <- F
            # for(i in 1:1){
            # yy <- yy0
            # dat <- list(xx=xx,yy=yy)
            # load(paste0(paste0('./dat/chr',chr,'/dat_',j,'.RData')))
            #train.x <- t(t(xx)-apply(xx,2,mean))
            #train.y <- (yy)-mean(yy)
            if(iter==5){
              allowempty <- F
            }
            model <- soda_main_allback(train.x,train.y, colnames(train.x),
                                       fixset = fixset, FALSE, FALSE, gam, nF0, allowempty)
            fixset <- unique(c(fixset,model$`Select Set`))
            gam <- gam-0.2
          }

          #model <- soda_main_allback(train.x, train.y,gam=gam,minF0=nF0)
          res <- list()
          res$BIC <- model$EBIC
          #res$Var <- model$Var
          #res$Term <- model$Term
          #res$best_term <- model$Term[[which.min(model$BIC)]]
          res$best_term <- model$`Select Set`
          #model$final_Term
          soda.train <- data.frame(V1=train.y, create_pmatrix_from_terms(train.x,   res$best_term))
          # soda.test <- data.frame(V1=test.y, create_pmatrix_from_terms(test.x,   res$best_term))
          tt <- stats::lm(V1~.,data=soda.train)
          res$coef <- coef(tt)
          res$single <- setdiff(res$best_term,fix[[j]])
          res$single.snp <- rsid00[res$single]
          res$common.snp <- rsid00[fix[[j]]]
          #  pp <- predict(tt,newdata=soda.test)
          #  res$corr <-  get.cor(soda.test$V1,pp)
          return(res)
        }}}}else{
          return(NA)
        }
}
