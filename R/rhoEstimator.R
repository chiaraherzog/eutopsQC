#' rhoEstimator
#'
#' Estimates signal to noise ratio
#' @param Mset MSet
#' @param anno EPIC
#' @return signal-to-noise (rho) value
#' @export

rhoEstimator <- function(Mset = Mset,
                          beta_e = 0.3309135,
                          beta_s1 = 0.019570,
                          beta_s2 = 0.970650,
                          anno = "EPIC"){
  
  
  array <- anno
  # load annotation infromation 
  if(anno == "EPIC"){
    data("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
    anno = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  } else if (anno == "mouse"){
    suppressPackageStartupMessages(library(IlluminaMouseMethylationanno.12.v1.mm10))
    suppressPackageStartupMessages(library(IlluminaMouseMethylationmanifest))
    anno <- getAnnotation(IlluminaMouseMethylationanno.12.v1.mm10)
  } else {
    data("IlluminaHumanMethylation450kanno.ilmn12.hg19")
    anno = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  }

  ratioSet <- ratioConvert(Mset, what = "both", keepCN = TRUE)
  beta_raw <- getBeta(ratioSet)
  
  cat(paste("raw beta extracted for samples (n = ", dim(Mset)[2], ")", sep="" ))
  
  # step 1 -------- calculate type I and type II density   
  
  datI <- data.frame(matrix(NA,nrow=512,ncol=ncol(beta_raw)))
  datII <- data.frame(matrix(NA,nrow=512,ncol=ncol(beta_raw)))
  
  colnames(datI) <- colnames(beta_raw)
  colnames(datII) <- colnames(beta_raw)
  
  ind <- match(rownames(beta_raw), anno$Name)
  ind <- ind[!is.na(ind)]
  
  pB <- txtProgressBar(min=1,max=ncol(beta_raw), width =50L, style = 3)
  for(j in 1:ncol(beta_raw)){
    setTxtProgressBar(pB, j)
    
    if(sum(is.na(beta_raw[anno$Type[ind]=='II',j]))<0.10*nrow(beta_raw)){
      dII <- density(beta_raw[anno$Type[ind]=='II',j],
                     na.rm=T,
                     bw=0.01,
                     n=512,
                     from=0,
                     to=1)
      
      dI <- density(beta_raw[anno$Type[ind]=='I',j],
                    na.rm=T,
                    bw=0.01,
                    n=512,
                    from=0,
                    to=1)
      
      datI[,j] <- dI$y  
      datII[,j] <- dII$y  
    }
  }
  
  close(pB)
  
  datI <- cbind(density_x=dI$x, datI)
  datII <- cbind(density_x=dII$x, datII)
  
  print("density distribution for type I and II probes calculated")
  
  #step 2 -------- calculate modes based on the density curve 
  dx <- datI[,1]
  
  datI[,1] <- NULL
  datII[,1] <- NULL
  
  # detect location of both modes
  modes_loc <- data.frame(matrix(NA, ncol=2, nrow=ncol(datII)))
  rownames(modes_loc) <- colnames(datII)
  colnames(modes_loc) <- c('mode0','mode1')
  
  for(i in 1:ncol(datII)){
    dy <- diff(datII[,i])
    d2y <- diff(dy)
    
    ind <- which((diff(sign(dy))!=0))
    ind_modes <- ind[d2y[ind]<0]
    
    if(length(ind_modes>=2)){
      # keep the two biggest peaks
      modes_loc[i,] <- sort(dx[ind_modes[order(datII[,i][ind_modes], decreasing = TRUE)[1:2]]])
    }
  }
  print("modes have been calculated")
  
  #step 3 -------- rho estimation  
  beta_s <- c(beta_s1, beta_s2)
  
  x <- seq(0,100,by=0.01)
  
  b0 <- (x*beta_s[1] + beta_e)/(1+x)
  b1 <- (x*beta_s[2] + beta_e)/(1+x)
  
  rho <- rep(NA,nrow(modes_loc))
  
  for(i in 1:nrow(modes_loc)){
    
    L <- (b0 - modes_loc$mode0[i])^2 + (b1 - modes_loc$mode1[i])^2
    
    if(length(x[which.min(L)]) > 0){
      rho[i] <- x[which.min(L)]
    }
  }
  
  rdat <- data.frame(modes_loc, rho=rho)  
  
  return(rdat)
  
}

