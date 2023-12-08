agtc <- function(a1,a2,b1,b2){
  sig <- rep(1,length(a1))
  for(i in 1:length(a1)){
    if((is.na(a1[i]))||(is.na(a2[i]))||(is.na(b1[i]))||(is.na(b2[i]))){
      sig[i] <- 0
    }
    else if((a1[i]==b1[i])&(a2[i]==b2[i])){
      sig[i] <- 1
    }else if((a1[i]==b2[i])&(a2[i]==b1[i])){
      sig[i] <- -1
    }
    else{
      if(b1[i]=="A"){temp1 <- "T"
      }else if(b1[i]=="T"){temp1 <- "A"
      }else if(b1[i]=="G"){temp1 <- "C"
      }else if(b1[i]=="C"){temp1 <- "G"
      }else{temp1 <- 'NA'}

      if(b2[i]=="A"){temp2 <- "T"
      }else if(b2[i]=="T"){temp2 <- "A"
      }else if(b2[i]=="G"){temp2 <- "C"
      }else if(b2[i]=="C"){temp2 <- "G"
      }else{temp2 <- 'NA'}

      if((a1[i]==temp1)&(a2[i]==temp2)){
        sig[i] <- 1
      }else if((a1[i]==temp2)&(a2[i]==temp1)){
        sig[i] <- -1
      }else{ sig[i] <- 0}
    }
    # if(i %% 10000==0){
    #   print(i)
    # }
  }
  #cat(length(sig)," variants in total. \n")
  #cat("Skip ", length(which(sig==0)), " variants due to allele code mismatch. \n")
  return(sig)
}

#' Run association test for multiple traits
#' @param list.stats a list of summary statistics, or a list of files to the summary statistics
#' @param x twas data list
#' @param list.ts.eQTL a list of tissue-specific eQTL for each tissue
#' @param ncores number of cores in parallel computing.
#' @return a list of association tests for all traits
#' @export
twas.multi.trait = function(list.stats, x, list.ts.eQTL, ncores = 1) {
  if (ncores == 1) {
    res = lapply(list.stats, function(s) twas.single.trait(s, x, list.ts.eQTL))
  } else {
    res = mclapply(list.stats, twas.single.trait, x = x, list.ts.eQTL = list.ts.eQTL,
             mc.cores = 80, mc.preschedule = F)
  }
  res
}

#' An association test in a tissue-specific manner for a phenotype (trait).
#' @param stats summary statistics for a phenotype. It can be a filename to be read or a data.frame.
#' Suppose it has 5 columns with colnames \code{rsid}, \code{a1}, \code{a2}, \code{chr}, \code{z}.
#' @param x twas data list
#' @param list.ts.eQTL a list of \code{ts.eQTL}, each element corresponds to a tissue
#' @return z-statistics and p-values. Both are of size n_gene x n_tissue.
#' @importFrom utils read.table
#' @export
twas.single.trait = function(stats, x, list.ts.eQTL) {
  if (typeof(stats) == "character") {
    cat("Try to read summary statistics from ", stats, "\n")
    stats = read.table(stats, header = T)
  }
  ntissue = length(list.ts.eQTL)
  ngene = length(x$pos)
  z_all <- p_all <- matrix(0, nrow = ngene, ncol = ntissue)
  for (i in 1:ntissue){
    rr <- twas.res(stats, list.ts.eQTL[[i]], x$dat, x$pos)
    z_all[, i] <- rr$z_g
    p_all[, i] <- rr$p_g
  }
  colnames(z_all) <- colnames(p_all) <- names(x$E)
  rownames(z_all) <- rownames(p_all) <- x$E.info$gene
  list(zstat = z_all, pval = p_all)
}

## stats:: rsid, a1, a2, z
#' MTWAS main function
#' @param stats summary statistics
#' @param soda.all.imp a list of all selected important genes
#' @param dat description
#' @param pos xx
#' @importFrom propagate cor2cov
#' @importFrom stats cov cov2cor pnorm prcomp
twas.res <- function(stats, soda.all.imp,
                     dat, pos # load by fn11,
                     #hm3, # snpinfo
                     #ct.all # cell type
                     ){
  #stats <- stats[stats$chr==chr,]
  #ct <- ct.all[kct]
  #cv <- 1 ## only for index
  # snp.gtex <- dat$bim
  # snp.gtex$snp_hg38 <- paste0(snp.gtex$V1,'_',snp.gtex$V4)
  # hm3$snp_hg38 <-  paste0(hm3$chr,'_',hm3$pos_hg38)
  # rr=substr(hm3$SNP,1,1)
  # hm3 <- hm3[which(rr=='r'),]
  # hm3 <- hm3[!duplicated(hm3$snp_hg38),]
  # dat$bim$rsid <- hm3$SNP[match(snp.gtex$snp_hg38, hm3$snp_hg38)]
  #
  dat$bim$rsid = dat$bim$V2 # temporary use
  bim.order <- 1:nrow(dat$bim) ## change to rsid
  stats <- stats[!duplicated(stats$rsid),]
  dat.bim <- merge(dat$bim, stats, by='rsid',all.x = T)
  dat.bim <- dat.bim[order(bim.order), ]
  sig <- agtc(dat.bim$V5, dat.bim$V6, dat.bim$a1, dat.bim$a2)
  zz <- dat.bim$z*sig
  zz[is.na(zz)] <- 0
  z_g <- c()
  for(l in 1:length(pos)){ ### for each gene
    if(length(soda.all.imp[[l]])>1){

      w <- soda.all.imp[[l]]$coef[-1]
      term <- soda.all.imp[[l]]$best_term
      temp <- dat$bed[,pos[[l]][term]]
      z_l <- zz[pos[[l]][term]]
      temp[is.na(temp)] <- 0
      zz.use<- which(z_l!=0)
      zz.use
      w
      z_l
      #if(length(w)>1){
      if(length(zz.use)>1){  #### main twas functions!!
          w <- w[zz.use]
          temp <- temp[,zz.use]
          z_l <- z_l[zz.use]
          covm <- cov(temp)
          varm <- diag(covm)
          corm <- cov2cor(covm)
          corm1 <- 0.9*corm+0.1*diag(nrow(corm))
          covm <- cor2cov(corm1,var=varm)
          sig_g2 <- t(w)%*%covm%*%w
          z_g[l] <- sum(w*c(sqrt(diag(covm)))*z_l)/sqrt(sig_g2)
      }else if(length(zz.use)==1){
        z_g[l] <-  z_l[zz.use]
      }else{
        z_g[l] <-  0
      }
    }else{
      z_g[l] <-  0
    }
  }
  # avoid directly modify the loaded data
  # E.info$z_g <- z_g
  # E.info$p_g <- pnorm(-abs(z_g))*2
  list(z_g = z_g, p_g = pnorm(-abs(z_g))*2)
}


#' Run association test with downloaded trained data
#' @param stats summary statistics for a phenotype. It can be a filename to be read or a data.frame, should include columns rsid, a1, a2, chr, z
#' @param twas_bim pre-downloaded file 
#' @param twas_eqtl pre-downloaded file 
#' @return mtwas summary statistics
#' @importFrom propagate cor2cov
#' @importFrom stats cov cov2cor pnorm prcomp
#' @export
run_mtwas_easy <- function(stats,
                           twas_bim,
                           #E.info,
                           twas_eqtl){
  twas_bim$rsid = twas_bim$V2 # temporary use
  bim.order <- 1:nrow(twas_bim) ## change to rsid
  stats <- stats[!duplicated(stats$rsid),]
  dat.bim <- merge(twas_bim, stats, by='rsid',all.x = T)
  dat.bim <- dat.bim[order(bim.order), ]
  sig <- agtc(dat.bim$V5, dat.bim$V6, dat.bim$a1, dat.bim$a2)
  zz <- dat.bim$z*sig
  zz[is.na(zz)] <- 0
  z_g <- c()
  for(l in 1:length(twas_eqtl)){ ### for each gene
    if(length(twas_eqtl[[l]])>1){
      
      w <- twas_eqtl[[l]]$coef[-1]
      term <- twas_eqtl[[l]]$pos_index
      covm0 <- twas_eqtl[[l]]$covm
      z_l <- zz[term]
      zz.use<- which(z_l!=0)
      zz.use
      w
      z_l
      #if(length(w)>1){
      if(length(zz.use)>1){  #### main twas functions!!
        w <- w[zz.use]
        #temp <- temp[,zz.use]
        z_l <- z_l[zz.use]
        covm <- covm0[zz.use,zz.use]
        varm <- diag(covm)
        corm <- cov2cor(covm)
        corm1 <- 0.9*corm+0.1*diag(nrow(corm))
        covm <- cor2cov(corm1,var=varm)
        sig_g2 <- t(w)%*%covm%*%w
        z_g[l] <- sum(w*c(sqrt(diag(covm)))*z_l)/sqrt(sig_g2)
      }else if(length(zz.use)==1){
        z_g[l] <-  z_l[zz.use]
      }else{
        z_g[l] <-  0
      }
    }else{
      z_g[l] <-  0
    }
  }
  # avoid directly modify the loaded data
  # E.info$z_g <- z_g
  # E.info$p_g <- pnorm(-abs(z_g))*2
  res <- data.frame(gene=names(twas_eqtl),z_g=z_g,p_g=pnorm(-abs(z_g))*2)
  #list(z_g = z_g, p_g = pnorm(-abs(z_g))*2)
  return(res)
}


                 
