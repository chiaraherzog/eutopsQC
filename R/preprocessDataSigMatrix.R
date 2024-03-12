#' preprocessData
#'
#' Pipeline to load and process raw data from the Illumina array
#' Authors: Chiara Herzog
#' Contact: chiara.herzog@uibk.ac.at
#' note to fix? "ChaMP_normalization/" empty folder is left in location of script where preprocessData function is called from
#'
#' @param input Path to the input signal matrix
#' @param output Path to output directory (this is where beta file will be saved)
#' @param report Path to report directory (this is where the report will be saved to)
#' @param array Methylation array type, EPIC by default
#' @param cores cores to be used for ChAMP normalisation. 4 by default
#' @param pheno phenotypic file. NULL by default, should be provided
#' @param overwrite overwrite existing output folder. FALSE by default to prevent any accidental overwriting.
#' @param beta.subset.compatible EPICv2 only, should EPIC v1/V2 compatible subset be provided?
#' @param meth column name (grep) of methylated signal
#' @param unmeth column name (grep) of unmethylated signal
#' @param detPname detP name (grep)
#' @param sep separator for index file, '\t' by default
#' @param run.name NULL by default
#' @return preprocessed beta matrix and QC report
#' @export

preprocessDataSigMatrix <- function(input = "",
                                    output = "",
                                    report = "",
                                    array = "EPIC",
                                    pheno = NULL,
                                    cores = 4,
                                    overwrite = FALSE,
                                    run.name = NULL,
                                    beta.subset.compatible = F,
                                    meth = 'Signal_A',
                                    unmeth = 'Signal_B',
                                    detPname = 'Detection.Pval',
                                    sep = '\t'){

  # Install packages (if missing)
  eutopsQC::installBiocDependencies(eutopsQC::packageList)

  # create Output folder
  if(dir.exists(output)){
    NULL
  } else {
    dir.create(output, recursive = TRUE)
    cat(paste0("Output folder ", output, " created.\n"))
  }

  if("beta_merged.Rdata" %in% list.files(output) & overwrite == FALSE){
    stop("Output folder is not empty, continuing would overwrite existing results.")
  }

  # create Report folder and subdirectories
  if(dir.exists(report)){
    NULL
  } else {
    dir.create(report, recursive = TRUE)
    cat(paste0("Report folder ", report, " created.\n"))
  }

  # append "/" for folder creation
  if(substr(report, nchar(report), nchar(report)) != "/"){
    report <- paste0(report, "/")
  }

  log <- paste0(report, "Log")
  dir.create(log)

  sink(paste0(report, "/Log/log.txt"), split = T)

  # Warning messages
  if(!file.exists(input)) stop('Input file does not exist')
  if(!exists("pheno") | is.null(pheno)) stop('Please provide pheno')
  if(!dir.exists(output)) stop('Output directory does not exist')
  if(!dir.exists(report)) stop('Report directory does not exist')

  if(is.null(run.name)){
      cat("Please provide run name:")
      run.name <- readline()
  }

  # Check if pheno files are present
  bnames <- as.character(stringr::str_split(readLines(input, n = 1), "\t|,", simplify = T))
  pheno$basename <- pheno$sampleid

  if(!all(grep(paste0(pheno$sampleid, collapse = '|'), bnames))){
    warning('Not all pheno ids found in signal matrix.')
  }


  # Begin pipeline
  cat('Beginning signal matrix load and preprocessing pipeline...\n\n')
  cat('Current time:',as.character(Sys.time()),'\n')
  cat('cwd:',getwd(),'\n\n')

  suppressPackageStartupMessages(library(minfi))
  suppressPackageStartupMessages(library(ChAMP))
  suppressPackageStartupMessages(library(impute))
  suppressPackageStartupMessages(library(ewastools))
  suppressPackageStartupMessages(library(stringr))
  suppressPackageStartupMessages(library(kableExtra))
  suppressPackageStartupMessages(library(EpiDISH))
  suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
  suppressPackageStartupMessages(library(IlluminaHumanMethylation450kanno.ilmn12.hg19))

  # Print arguments
  cat('Input arguments:\n')
  cat('Array type = ', array, '\n')
  cat('Path to file =', input, '\n')
  cat('Path to output =',output,'\n')
  cat('Path to log file =',log,'\n')

  # Define global thresholds
  INTENSITY_THRESHOLD <- 9.5     # minimum median intensity required
  DETECTION_P_THRESHOLD <- 0.01  # maximum detection p-value
  FAILED_PROBE_THRESHOLD <- 0.1   # maximum proportion of failed probes per sample

  cat('INTENSITY_THRESHOLD =',INTENSITY_THRESHOLD,'\n')
  cat('DETECTION_P_THRESHOLD =',DETECTION_P_THRESHOLD,'\n')
  cat('FAILED_PROBE_THRESHOLD =',FAILED_PROBE_THRESHOLD,'\n\n')

  # Set up detailed array anno
  if(array == '450k') {
    array.det = 'IlluminaHumanMethylation450k'
    anno = 'ilmn12.hg19'
  } else if(array == 'EPIC'){
    array.det = 'IlluminaHumanMethylationEPIC'
    anno = 'ilm10b4.hg19'
  } else if(array == 'EPICv2'){
    array.det = 'IlluminaHumanMethylationEPICv2'
    anno = '20a1.hg38'
  }

  # Set up log
  log_data <- list(plate_name = run.name, # name of plate
                   n_samples=NA,          # no of samples on plate
                   n_probes=NA,           # no of probes in matrix
                   n_samples_removed=NA,  # no of samples that fail QC
                   rm_sample_list=NA      # list of any removed samples
  )

  # Unmeth = A, Meth = B (based on order of other signal intensity files)
  Mset <- minfi::readGEORawFile(filename = input,
                         sep = sep,
                         Uname = Uname,
                         Mname = Mname,
                         array = array.det,
                         annotation = anno,
                         showProgress = F,
                         row.names = 2)

  qc <- getQC(Mset)
  cat(' done\n')

  # Filter any samples with median (un)methylated intensity less than threshold
  low_intensity_samples <- rownames(qc)[qc$mMed<INTENSITY_THRESHOLD | qc$uMed<INTENSITY_THRESHOLD]

  # Read in detP
  select <- sort(grep(detPname, bnames))
  detP <- data.table::fread(input,
                            sep = sep,
                            select = select)

  # Filter samples with too many failed probes
  failed_samples <- colnames(detP)[colSums(detP>DETECTION_P_THRESHOLD) > (nrow(detP) * FAILED_PROBE_THRESHOLD)]
  samples_to_remove <- unique(c(low_intensity_samples, failed_samples))
  log_data$n_samples_removed <- length(samples_to_remove)

  rm_ind <- match(samples_to_remove, colnames(Mset))

  if(length(rm_ind)>0){
    Mset_filtered <- Mset[,-rm_ind]
    cat('\nRemoved',length(rm_ind),'failed sample(s):\n')
    for(r in 1:length(samples_to_remove)){
      cat(samples_to_remove[r],'\n')
    }
    log_data$rm_sample_list <- samples_to_remove #log samples removed
  } else {
    Mset_filtered <- Mset
    cat('\nNo samples removed... \n')
  }
  cat("\n")

  # Convert to RatioSet and then beta
  RSet <- minfi::ratioConvert(Mset, what = "both", keepCN = TRUE)
  rm(Mset);gc()
  beta <- minfi::getBeta(RSet)
  beta <- na.omit(beta)
  rm(RSet);gc()

  log_data$n_samples <- ncol(beta)
  log_data$n_probes <- nrow(beta)

  # check that no bad samples are present:
  beta <- beta[,!colnames(beta) %in% samples_to_remove]

  # Probe bias correction using BMIQ
  if (grepl("v2", array, ignore.case = T)){
    cat('Begin BMIQ normalisation using', array, 'version...\n')
    beta_ssNOOB_filtered_norm <- invisible(eutopsQC:::champ.normv2(beta=beta,
                                                                   arraytype=array,
                                                                   cores=cores))
    cat('done\n')
  } else {
    cat('Begin BMIQ normalisation using', array, 'version...\n')
    beta_ssNOOB_filtered_norm <- invisible(champ.norm(beta=beta,
                                                      arraytype=array,
                                                      cores=cores))
    cat('done\n')
  }

  DETECTION_P_THRESHOLD <- 0.01   # maximum detection p-value
  FAILED_SAMPLE_THRESHOLD <- 0.1   # maximum proportion of failed samples per probe
  cat('DETECTION_P_THRESHOLD =',DETECTION_P_THRESHOLD,'\n')
  cat('FAILED_SAMPLE_THRESHOLD =',FAILED_SAMPLE_THRESHOLD,'\n\n')

  # Remove any failed probes
  failed_probes <- rownames(detP)[rowSums(detP>DETECTION_P_THRESHOLD)>(ncol(detP)*FAILED_SAMPLE_THRESHOLD)]

  # combine all CpGs to remove
  if (grepl("v2", array, ignore.case = T)){
    rm_names <- unique(c(chrY_0_M_names_v2,non_CpG_names_v2,snp_names_v2,zhou_list_v2,failed_probes))
    cat('Removing', length(chrY_0_M_names_v2),'chrY, chr0, chrM probes\n')
    cat('Removing', length(non_CpG_names_v2),'non-CpG probes\n')
    cat('Removing', length(snp_names_v2),'SNP probes flagged previously on v1\n')
    cat('Removing', length(failed_probes),'failed probes (detP)\n')
    cat('Removing', length(zhou_list_v2),'remaining Zhou SNP probes on v2\n\n')
  } else {
    rm_names <- unique(c(chrY_names,non_CpG_names,snp_names,zhou_list,failed_probes))
    cat('Removing', length(chrY_names),'chrY probes\n')
    cat('Removing', length(non_CpG_names),'non-CpG probes\n')
    cat('Removing', length(snp_names),'SNP probes\n')
    cat('Removing', length(failed_probes),'failed probes (detP)\n')
    cat('Removing', length(zhou_list),'Zhou SNP probes\n\n')
  }


  # remove probes and any bad samples
  ind.row <- match(rownames(beta_ssNOOB_filtered_norm),rm_names)
  beta_filtered <- beta_ssNOOB_filtered_norm[is.na(ind.row),]

  # Imputation
  r <- rowSums(is.na(beta_filtered))
  if(sum(r > 0.8*ncol(beta_filtered))>0){
    beta_filtered <- beta_filtered[(r < 0.8*ncol(beta_filtered))>0, ]
    cat('!! Probe failure in more than 80% of samples for',
        sum(r > 0.8*ncol(beta_filtered)),
        'probes on plate',p,'(removed from beta matrix)\n')
  }

  out <- capture.output(beta_filtered_imputed <- impute::impute.knn(beta_filtered,k=10,rowmax=0.8)$data)
  rm(beta_plate_filtered);invisible(gc())

  beta <- beta_filtered_imputed

  cat('\nBeta matrix has',
      nrow(beta),
      'CpGs and',
      ncol(beta),
      'samples\n\n')

  # Save beta density as plotly plot
  cat('Begin plot beta densities...')
  d <- density(beta[,1], na.rm=TRUE, bw=0.02,
                 from = -0.05, to = 1.05)
  p <- plotly::plot_ly(x=d$x,y=d$y,
                 mode = 'line',
                 type='scatter',
                 text=colnames(beta)[1])

  for(j in 2:ncol(beta)){
      d <- density(beta[,j], na.rm=TRUE, bw=0.02,
                   from = -0.05, to = 1.05)
      p <- p |> plotly::add_trace(x=d$x,y=d$y,
                           mode = 'line',
                           type='scatter',
                           text=colnames(beta)[j])
  }

  p <- p |> plotly::layout(title = run.name,
                        yaxis = list(title = 'Density'),
                        xaxis = list(title = 'Beta'),
                        showlegend=FALSE) |>
      plotly::config(displayModeBar = FALSE)
  cat(' done\n')

  # Save outputs
  save(qc, file=paste(log,'/',
                        run.name,'_qc.Rdata',sep=''))

  save(detP, file=paste(log,'/', run.name, '_detP.Rdata',sep=''))

  save(p, file=paste(log,'/', run.name, '_beta_density_plot.Rdata',sep=''))

  save(log_data, file=paste(log,'/', run.name, '_log_data.Rdata',sep=''))

  save(beta, file = paste0(output, '/beta.Rdata', sep = ''))

  if(array == 'EPICv2' & beta.subset.compatible == T){
    # create a EPIC v1 v2 compatible subset of beta_merged
    beta_merged_compatible <- subset_versionshared_CpGs(beta_merged, array)
  }

  cat(' done\n\n')

  # Create a report
  cat('Creating RMarkdown Report...\n')
  file.copy(from = paste0(system.file("rmd", "_site_SigMatrix.yml",
                                      package = "eutopsQC")),
            to = paste0(report, "_site.yml"))

  # file.copy(from = "~/Documents/Work/Code/preprocessing/eutopsQC/inst/rmd/_site_SigMatrix.yml",
  #           to = paste0(report, "_site.yml"))

  rmarkdown::render(input = paste0(system.file("rmd", "index_SigMatrix.Rmd",
                                               package = "eutopsQC")),
                    output_file = paste0(report, "index.html", sep = ""))

  rmarkdown::render(input = paste0(syste.file("rmd", "1-plate-summary_SigMatrix.Rmd",
                                              package = "eutopsQC")),
                    output_file = paste0(report, "1-plate-summary.html"))

  rmarkdown::render(input = paste0(system.file("rmd", "2-qc_SigMatrix.Rmd",
                                               package = "eutopsQC")),
                    output_file = paste0(report, "2-qc.html", sep = ""))

  rmarkdown::render(input = paste0(system.file("rmd", "3-beta-distributions_SigMatrix.Rmd",
                                               package = "eutopsQC")),
                    output_file = paste0(report, "3-beta-distributions.html", sep = ""))


  if(!identical(pheno$basename, colnames(beta))){

    if(all(grep(paste0(pheno$basename, collapse = "|"), colnames(beta)))){
    ind <-grep(paste0(pheno$basename, collapse = "|"), colnames(beta))
    beta <- beta[,ind]
    colnames(beta)<- pheno$basename
    }
  }

  if (!grepl("v2", array, ignore.case = T)){
    # note age, ic and smk indices calculated from beta

    if(exists("pheno") & array != "mouse"){
      out <- epidish(beta.m = beta,
                     ref.m = centEpiFibIC.m,
                     method = "RPC")$estF
      ind <- match(pheno$basename, rownames(out))
      pheno$ic <- out[ind,3]
    }

    rmarkdown::render(input = paste0(system.file("rmd", "4-age-ic-qc_SigMatrix.Rmd",
                                                 package = "eutopsQC")),
                      output_file = paste0(report, "4-age-ic-smk.html", sep = ""))

  } else{
    # EPIC v2
    # age, ic and smk indices calculated from beta_merged_compatible (based on v1 IDs)
    # suggest to add retrained, v1.v2 compatible age index to WID clocks packages

    if(exists("pheno")){
      out <- epidish(beta.m = beta_merged_compatible,
                     ref.m = centEpiFibIC.m,
                     method = "RPC")$estF
      ind <- match(pheno$basename, rownames(out))
      pheno$ic <- out[ind,3]
    }

    local({
      beta <- beta_merged_compatible
      rmarkdown::render(input = paste0(system.file("rmd", "4-age-ic-qc_SigMatrix.Rmd",
                                                   package = "eutopsQC")),
                        output_file = paste0(report, "4-age-ic-smk.html", sep = ""))
    })

  }

  rmarkdown::render(input = paste0(system.file("rmd", "5-pca.Rmd",
                                               package = "eutopsQC")),
                      output_file = paste0(report, "5-pca.html", sep = ""))

  cat('Session info:\n\n')
  print(sessionInfo())
  sink()

}
