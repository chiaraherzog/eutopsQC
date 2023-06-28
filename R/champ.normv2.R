#' champ.normv2
#'
#' ChAMP normalisation for the Illumina Human Methylation EPIC v2 array. Like the ChAMP champ.norm package, it can use SWAN, funNorm, PBC, or BMIQ normalisation. Depending on normalisation algorithm, different inputs may be needed (beta matrix, rgSet, mset). Please refer to the ChAMP package manual.
#'
#' @param beta Beta Matrix
#' @param rgSet RGset (optional)
#' @param mset Mset (optional)
#' @param cores used cores, default = 4
#' @return Normalised beta matrix
#' @export

champ.normv2 <- function(beta = myLoad$beta, rgSet = myLoad$rgSet, mset = myLoad$mset,
                        resultsDir = "./CHAMP_Normalization/", method = "BMIQ",
                        plotBMIQ = FALSE, arraytype = "mouse", cores = 4)
{

  message("[===========================]")
  message("[>>>>> ChAMP.NORM START <<<<<<]")
  message("-----------------------------")
  if (!file.exists(resultsDir))
    dir.create(resultsDir)
  message("champ.norm Results will be saved in ", resultsDir)
  message("[ SWAN method call for BOTH rgSet and mset input, FunctionalNormalization call for rgset only , while PBC and BMIQ only needs beta value. Please set parameter correctly. ]\n")

  # load probeInfo
  library(IlluminaHumanMethylationEPICv2manifest)
  library(ChAMP)
  library(doParallel)
  library(RPMM)
  probeInfo <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38) |>
    as.data.frame() |>
    dplyr::select(Name, Type) |>
    dplyr::rename(Design = Type) |>
    dplyr::mutate(Design = ifelse(Design == "I", 1, 2))
  design.v <- as.numeric(probeInfo$Design[match(rownames(beta),probeInfo$Name)])


  if (method == "SWAN") {
    beta.p = getBeta(preprocessSWAN(rgSet, mset), "Illumina")
  }

  else if (method == "BMIQ") {
    message("<< Normalizing data with BMIQ Method >>")
    message("Note that,BMIQ function may fail for bad quality samples (Samples did not even show beta distribution).")

    if (min(beta, na.rm = TRUE) == 0) {
      beta[beta == 0] <- 1e-06
      message("Zeros in your dataset have been replaced with 0.000001\n")
    }
    current_cwd <- getwd()
    if (plotBMIQ)
      setwd(resultsDir)
    if (cores > 1) {
      if (cores > detectCores())
        cores <- detectCores()
      message(cores, " cores will be used to do parallel BMIQ computing.")
      cl <- makeCluster(cores)
      registerDoParallel(cl)
      beta.p <- foreach(x = 1:ncol(beta), .combine = cbind) %dopar% ChAMP:::champ.BMIQ(beta[,
                                                                                            x], design.v, sampleID = colnames(beta)[x],
                                                                                       plots = plotBMIQ)$nbeta
      stopCluster(cl)
    }
    else {
      beta.p <- sapply(1:ncol(beta), function(x) ChAMP:::champ.BMIQ(beta[,
                                                                         x], design.v, sampleID = colnames(beta)[x],
                                                                    plots = plotBMIQ)$nbeta)
    }
    setwd(current_cwd)
  }
  else if (method == "PBC") {
    message("<< Normalizing data with PBC Method >>")
    if (min(beta, na.rm = TRUE) == 0) {
      beta[beta == 0] <- 1e-06
      message("Zeros in your dataset have been replaced with 0.000001\n")
    }
    beta.p = ChAMP:::DoPBC(beta, design.v)
  }
  else if (method == "FunctionalNormalization") {
    if (is.null(rgSet))
      stop("rgSet not found, it is required for FunctionalNormalization")
    if (arraytype == "EPIC")
      rgSet@annotation[2] <- "ilm10b4.hg19"
    beta.p <- getBeta(preprocessFunnorm(rgSet))[rownames(beta),
    ]
  }
  else {
    stop("Please Select Normalization Method from: BMIQ,PBC, FunctionalNormalization and SWAN.")
  }
  rownames(beta.p) <- rownames(beta)
  colnames(beta.p) <- colnames(beta)
  message("[>>>>> ChAMP.NORM END <<<<<<]")
  message("[===========================]")
  message("[You may want to process champ.SVD() next.]\n")
  return(beta.p)
}
