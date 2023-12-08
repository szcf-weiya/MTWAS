#' Run association test with downloaded trained data
#' @param stats summary statistics for a phenotype. It can be a filename to be read or a data.frame, should include columns rsid, a1, a2, chr, z
#' @param twas_bim pre-downloaded file 
#' @param twas_eqtl pre-downloaded file 
#' @return mtwas summary statistics
#' @export

### to use the downloaded gtex dataset and run mtwas directly

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

#aa=run_mtwas_easy(summary_stats,twas_bim,twas_eqtl)