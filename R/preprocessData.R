#' preprocessData
#'
#' Pipeline to load and process raw data from the Illumina array
#' Authors: James E. Barrett, Chiara Herzog
#' Contact: james.barrett@ucl.ac.uk, chiara.herzog@uibk.ac.at
#' note to fix? "ChaMP_normalization/" empty folder is left in location of script where preprocessData function is called from
#' 
#' @param input Path to the input directory (raw IDAT folder)
#' @param output Path to output directory (this is where beta file will be saved)
#' @param report Path to report directory (this is where the report will be saved to)
#' @param array Methylation array type, EPIC by default
#' @param cores cores to be used for ChAMP normalisation. 4 by default
#' @param pheno phenotypic file. NULL by default. If pheno file is provided, it must have a column named basename with identical names to IDAT files. File can be provided as Rdata, RDS, .txt or .csv. Each column will be included in the batch effect correction. Ensure that no column exists with entirely identical variables
#' @param by.dir process by directory instead. FALSE by default. The result will not be different, but by.dir = T can be slower (for smaller projects), yet it is recommended for large projects.
#' @param overwrite overwrite existing output folder. FALSE by default to prevent any accidental overwriting.
#' @param save.rs save SNP (rs) probe values. FALSE by default.
#' @return preprocessed beta matrix and QC report
#' @export

preprocessData <- function(input = "",
                           output = "",
                           report = "",
                           array = "EPIC",
                           pheno = NULL,
                           cores = 4,
                           by.dir = FALSE,
                           path_to_bad_sample_list = "",
                           overwrite = FALSE,
                           save.rs = FALSE,
                           create.shiny = F){
  
  # Install packages

  if (!require("devtools")){
    install.packages("devtools")
  }
  
  if (!require("ewastools")){
      devtools::install_github("hhhh5/ewastools")
  }
  
  # create Output folder
  if(dir.exists(output)){
    NULL
  } else {
    dir.create(output, recursive = TRUE)
    cat(paste0("Output folder ", output, " created.\n"))
  }
  if(grepl("beta_merged.Rdata", list.files(output)) & overwrite == FALSE){
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
  if(!dir.exists(input)) stop('Raw file directory does not exist')
  if(!dir.exists(output)) stop('Output directory does not exist')
  if(!dir.exists(report)) stop('Report directory does not exist')
  
  
  
  # Load in pheno file  if present. If pheno file present, check that basenames are in the files
  if(!is.null(pheno)){
    if(grepl(".csv", pheno)){
      pheno <- read.table(pheno, sep = ",")
    } else if(grepl(".txt", pheno)){
      pheno <- read.table(pheno, sep = "\t")
    } else if (grepl(".xlsx", pheno)){
      pheno <- readxl::read_xlsx(pheno, sheet = 1)
    } else if (grepl(".Rdata", pheno)){
      assign('pheno', get(load(pheno)))
    } else if (grepl(".Rds", pheno)){
      pheno <- readRDS(pheno)
    } else {
      stop("Pheno not in a compatible file (.csv, .txt, .xlsx, or .Rdata). Please amend the input file.")
    }
  }
  
  if(exists("pheno") && any(grepl(paste0(pheno$basename, collapse = "|"), list.files(input, pattern = ".idat", recursive = T, include.dirs = F)) == FALSE)){
    stop("Pheno names not entirely overlappying with basenames")
  }
  
  # Begin pipeline
  cat('Beginning idat load and preprocessing pipeline...\n\n')
  cat('Current time:',as.character(Sys.time()),'\n')
  cat('cwd:',getwd(),'\n\n')
  
  suppressPackageStartupMessages(library(minfi))
  suppressPackageStartupMessages(library(ChAMP))
  suppressPackageStartupMessages(library(impute))
  suppressPackageStartupMessages(library(plotly))
  suppressPackageStartupMessages(library(ewastools))
  suppressPackageStartupMessages(library(stringr))
  suppressPackageStartupMessages(library(kableExtra))
  suppressPackageStartupMessages(library(EpiDISH))   
  suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))   
  suppressPackageStartupMessages(library(IlluminaHumanMethylation450kanno.ilmn12.hg19))   
  
  # Print arguments
  cat('Input arguments:\n')
  cat('Array type = ', array, '\n')
  cat('Path to idat =', input, '\n')
  cat('Path to bad sample list =', path_to_bad_sample_list, '\n')
  cat('Path to output =',output,'\n')
  cat('Path to log file =',log,'\n')
  cat('Pheno present =', ifelse(is.null(pheno), "FALSE", "TRUE"), "\n")
  cat('By directory = ', ifelse(by.dir == FALSE, "FALSE", "TRUE"), "\n\n")
  
  # Define global thresholds
  INTENSITY_THRESHOLD <- 9.5     # minimum median intensity required
  DETECTION_P_THRESHOLD <- 0.01  # maximum detection p-value
  FAILED_PROBE_THRESHOLD <- 0.1   # maximum proportion of failed probes per sample
  
  cat('INTENSITY_THRESHOLD =',INTENSITY_THRESHOLD,'\n')
  cat('DETECTION_P_THRESHOLD =',DETECTION_P_THRESHOLD,'\n')
  cat('FAILED_PROBE_THRESHOLD =',FAILED_PROBE_THRESHOLD,'\n\n')
  
  # Initalise a list to log various parameters
  plates <- list.dirs(input,
                      full.names = FALSE, recursive = F)
  
  if(by.dir == F){ # keep entire folder as one "set"
    tmp <- stringr::str_split(input, "/", simplify = TRUE)
    plates <- tmp[length(tmp)-1]
    rm(tmp)
  }
  
  # by plate (or directory)
  for (i in 1:length(plates)){
    
    log_data <- list(plate_name = plates[i], # name of plate
                     n_samples=NA,                # no of samples on plate
                     n_probes=NA,                 # no of probes in matrix
                     n_samples_removed=NA,        # no of samples that fail QC
                     snp_outlier_metric=NA,       # from ewastools::snp_outliers
                     rm_sample_list=NA          # list of any removed samples
    )
    
    #----------------------------------------------------#
    # Load and extract data
    
    # Reading in sample sheet and targets
    cat('Begin load idats...')
    if (length(plates) == 1){
      RGset <- read.metharray.exp(base = input,
                                  verbose = T,
                                  force = TRUE,
                                  recursive = TRUE)
    } else {
      RGset <- read.metharray.exp(base = paste(input, plates[i], sep = ""),
                                  verbose = T,
                                  force = TRUE,
                                  recursive = TRUE)
    }
    
    # Add annotation for mouse array if required
    if(grepl("mouse", array, ignore.case = T)){
      suppressPackageStartupMessages(library(IlluminaMouseMethylationanno.12.v1.mm10))
      suppressPackageStartupMessages(library(IlluminaMouseMethylationmanifest))
      RGset@annotation <- c(array = "IlluminaMouseMethylation", annotation = "12.v1.mm10")
    }
    
    # Save shinyMethyl page
    if(create.shiny == T){
      suppressPackageStartupMessages(library(shinyMethyl))
      invisible(summary <- shinySummarize(RGset))
      
      save(summary, file = paste0(output, "/", plates[i], "_shinyMethyl.Rdata")) # This file can be opened later to run the shiny app
    }
    
    # Save RS (only if not mouse, otherwise will throw error)
    if(save.rs == T & !grepl("mouse", array, ignore.case = T)){
      rs <- getSnpBeta(RGset)
      save(rs, file = paste0(log, "/", plates[i], "_rs.Rdata"))
    }
    
    # Save the controls for later
    if (is(RGset, "rgDataSet")) {
      ctrls <- ENmix::getCGinfo(RGset, type = "ctrl")
    } else if (is(RGset, "RGChannelSet")) {
      ctrls <- minfi::getProbeInfo(RGset, type = "Control")
    }
    ctrls <- ctrls[ctrls$Address %in% rownames(RGset) & !is.na(ctrls$Color), ]
    ctrl_r <- assays(RGset)$Red[ctrls$Address, ]
    ctrl_g <- assays(RGset)$Green[ctrls$Address, ]
    # Extract detection p-values from RGset
    detP <- minfi::detectionP(RGset, type = "m+u")
    
    # Extract quality control information
    cat('Begin QC (median intensity) extraction...')
    Mset <- preprocessRaw(RGset)
    qc <- getQC(Mset)
    cat(' done\n')
    
    # Extract Rho
    cat('Calculating rho...\n')
    rho <- data.frame(matrix(NA, ncol=3))
    colnames(rho) <- c("mode0", "mode1", "rho")
    rho <- eutopsQC::rhoEstimator(Mset, anno = array)
    anno <- array
    cat('done\n')
    
    # Extract SNP outlier metric using ewastools
    cat('Begin QC (SNP outlier metric) extraction...')
    beta_snp <- getSnpBeta(RGset)
    genotypes <- call_genotypes(beta_snp, learn = FALSE, maxiter = 50)
    log_data$snp_outlier_metric <- snp_outliers(genotypes)
    cat(' done\n')
    
    # Filter any samples with median (un)methylated intensity less than threshold
    low_intensity_samples <- rownames(qc)[qc$mMed<INTENSITY_THRESHOLD | qc$uMed<INTENSITY_THRESHOLD]
    
    # Load any bad samples (e.g. flagged by lab team or ICH)
    if(path_to_bad_sample_list!=''){
      bad_sample_list <- tryCatch({
        dat <- read.csv(path_to_bad_sample_list, header = FALSE) #load csv containing any bad samples
        ind <- match(dat$V1,colnames(RGset)) # see if any are present in beta matrix
        colnames(RGset)[ind[!is.na(ind)]]    # return sample names of any that are present
      }, warning = function(warning_condition) {
        return(NULL)
      }, error = function(error_condition) {
        return(NULL)
      })
      if(is.null(bad_sample_list)){
        warning('Failed to read bad sample list')
        cat('\n!! Warning: Failed to read bad sample list\n')
      }
    } else {
      bad_sample_list <- NULL
    }
    
    # Filter samples with too many failed probes
    failed_samples <- colnames(detP)[colSums(detP>DETECTION_P_THRESHOLD) > (nrow(detP) * FAILED_PROBE_THRESHOLD)]
    
    samples_to_remove <- unique(c(low_intensity_samples,
                                  bad_sample_list,
                                  failed_samples))
    
    log_data$n_samples_removed <- length(samples_to_remove)
    
    rm_ind <- match(samples_to_remove, colnames(RGset))
    if(length(rm_ind)>0){
      RGset_filtered <- RGset[,-rm_ind]
      cat('\nRemoved',length(rm_ind),'failed sample(s):\n')
      for(r in 1:length(samples_to_remove)){
        cat(samples_to_remove[r],'\n')
      }
      log_data$rm_sample_list <- samples_to_remove #log samples removed
    } else {
      RGset_filtered <- RGset
      cat('\nRemoved 0 failed samples... \n')
    }
    cat('\n')     
    
    
    # Bias correction
    # Background intensity correction and dye bias correction
    cat('Begin ssNOOB preprocessing...')
    ssNOOB_filtered <- preprocessNoob(RGset_filtered, dyeCorr=TRUE, verbose=TRUE, dyeMethod='single')
    cat('done\n')
    
    # Extract corresponding beta values
    cat('Begin beta extraction...')
    beta_ssNOOB_filtered <- getBeta(ssNOOB_filtered)
    cat('done\n')
    log_data$n_samples <- ncol(beta_ssNOOB_filtered)
    log_data$n_probes <- nrow(beta_ssNOOB_filtered)
    
    # check that no bad samples are present:
    beta_ssNOOB_filtered <- beta_ssNOOB_filtered[,!colnames(beta_ssNOOB_filtered) %in% samples_to_remove]
    
    # Probe bias correction using BMIQ
    if(grepl("mouse", array, ignore.case = T)){
      cat('Begin BMIQ normalisation using', array, 'version...\n')
      beta_ssNOOB_filtered_norm <- invisible(IlluminaMouseMethylationmanifest::champ.normm(beta=beta_ssNOOB_filtered,arraytype=array,cores=cores))
      cat('done\n')
    } else {
      cat('Begin BMIQ normalisation using', array, 'version...\n')
      beta_ssNOOB_filtered_norm <- invisible(champ.norm(beta=beta_ssNOOB_filtered,arraytype=array,cores=cores))
      cat('done\n')
    }
    
    # Save beta density as plotly plot
    cat('Begin plot beta densities...')
    d <- density(beta_ssNOOB_filtered_norm[,1], na.rm=TRUE, bw=0.02,
                 from = -0.05, to = 1.05)
    p <- plot_ly(x=d$x,y=d$y,
                 mode = 'line',
                 type='scatter',
                 text=colnames(beta_ssNOOB_filtered_norm)[1])
    
    for(j in 2:ncol(beta_ssNOOB_filtered_norm)){
      d <- density(beta_ssNOOB_filtered_norm[,j], na.rm=TRUE, bw=0.02,
                   from = -0.05, to = 1.05)
      p <- p %>% add_trace(x=d$x,y=d$y,
                           mode = 'line',
                           type='scatter',
                           text=colnames(beta_ssNOOB_filtered_norm)[j])
    }
    
    p <- p %>%   layout(title = plates[i],
                        yaxis = list(title = 'Density'),
                        xaxis = list(title = 'Beta'),
                        showlegend=FALSE)%>%
      config(displayModeBar = FALSE)
    cat(' done\n')
    
    # Save outputs
    cat('Begin save outputs...')
    save(qc, file=paste(log,'/',
                        plates[i],'_qc.Rdata',sep=''))
    
    save(ctrl_g, ctrl_r, ctrls, file=paste(log,'/',
                                           plates[i],'_ctrl.Rdata',sep=''))
    
    if(length(plates) != 1){
      save(beta_ssNOOB_filtered_norm, file=paste(output,'/',
                                                 plates[i],'_beta_ssNOOB_filtered_norm.Rdata',sep=''))
    } 
    
    save(detP, file=paste(log,'/',
                          plates[i],'_detP.Rdata',sep=''))
    
    save(p, file=paste(log,'/',
                       plates[i],'_beta_density_plot.Rdata',sep=''))
    
    save(rho, file=paste(log,'/',
                         plates[i],'_rho.Rdata',sep=''))
    
    save(log_data, file=paste(log,'/',
                              plates[i],'_log_data.Rdata',sep=''))
    
    cat(' done\n\n')
  }
  
  rm(RGset, p, detP, log_data, qc, ssNOOB_filtered, beta_ssNOOB_filtered, 
     beta_snp, d, genotypes, Mset, probeInfoALL.lv, RGset_filtered, DETECTION_P_THRESHOLD,
     FAILED_PROBE_THRESHOLD, failed_samples, INTENSITY_THRESHOLD, i, j, rm_ind, samples_to_remove);invisible(gc())
  
  cat('Beginning plate combination pipeline...\n\n')
  DETECTION_P_THRESHOLD <- 0.01   # maximum detection p-value
  FAILED_SAMPLE_THRESHOLD <- 0.1   # maximum proportion of failed samples per probe
  cat('DETECTION_P_THRESHOLD =',DETECTION_P_THRESHOLD,'\n')
  cat('FAILED_SAMPLE_THRESHOLD =',FAILED_SAMPLE_THRESHOLD,'\n\n')
  
  for(p in plates){
    if(!file.exists(paste(log,'/',p,'_detP.Rdata',sep=''))){
      stop(paste('detP matrix for plate ',p,' is missing',sep=''))
    }
  }
  
  load(paste(log,'/',plates[1],'_detP.Rdata',sep=''))
  detP_merged <- detP
  
  cat(paste('Merged plate ',plates[1],'\n',sep=''))
  
  if(length(plates)>1){
    for(p in plates[2:length(plates)]){
      load(paste(log,'/',p,'_detP.Rdata',sep=''))
      
      detP_merged <- merge(detP_merged, detP,
                           by.x='row.names', by.y='row.names', sort=FALSE)
      
      # reassign the row names
      row.names(detP_merged) <- detP_merged$Row.names
      
      # remove this column created by the merge function
      detP_merged$Row.names <- NULL
      
      cat(paste('Merged plate ',p,'\n',sep=''))
    }
    cat('\n')
  }
  
  # Remove any failed probes
  failed_probes <- rownames(detP_merged)[rowSums(detP_merged>DETECTION_P_THRESHOLD)>(ncol(detP_merged)*FAILED_SAMPLE_THRESHOLD)]
  
  # Combine and save output
  DETECTION_P_THRESHOLD <- 0.01   # maximum detection p-value
  cat('DETECTION_P_THRESHOLD =',DETECTION_P_THRESHOLD,'\n\n')
  
  # combine all CpGs to remove
  rm_names <- unique(c(chrY_names,non_CpG_names,snp_names,zhou_list,failed_probes))
  cat('Removing', length(chrY_names),'chrY probes\n') 
  cat('Removing', length(non_CpG_names),'non-CpG probes\n') 
  cat('Removing', length(snp_names),'SNP probes\n') 
  cat('Removing', length(failed_probes),'failed probes (detP)\n') 
  cat('Removing', length(zhou_list),'Zhou SNP probes\n\n') 
  
  input_dir_list <- output
  
  #check beta matrices and detP matrices are present for each plate
  if(length(plates)>1){
    for(p in plates){
      if(!file.exists(paste(output,'/',p,'_beta_ssNOOB_filtered_norm.Rdata',sep=''))){
        stop(paste('beta matrix for plate ',p,' is missing',sep=''))
      }
      if(!file.exists(paste(log,'/',p,'_detP.Rdata',sep=''))){
        stop(paste('detP matrix for plate ',p,' is missing',sep=''))
      }
    }
  }
  
  # Load and combine the beta and detP matrices
  na_count <- 0
  beta_merged <- NULL
  
  for(p in plates){
    
    if(length(plates) != 1){ # load in extra file; if plates == 1, beta is still in environment
      load(paste(output,'/',p,'_beta_ssNOOB_filtered_norm.Rdata',sep=''))
    }
    
    load(paste(log,'/',p,'_detP.Rdata',sep=''))
    
    # remove probes and any bad samples   
    ind.row <- match(rownames(beta_ssNOOB_filtered_norm),rm_names)
    beta_plate_filtered <- beta_ssNOOB_filtered_norm[is.na(ind.row),]
    rm(beta_ssNOOB_filtered_norm);invisible(gc())
    
    # match detection p-values to beta matrix
    ind.row <- match(rownames(beta_plate_filtered),rownames(detP))
    ind.col <- match(colnames(beta_plate_filtered),colnames(detP))
    detP_plate <- detP[ind.row, ind.col]
    
    # replace failed data points with NA
    beta_plate_filtered[detP_plate > DETECTION_P_THRESHOLD] <- NA
    na_count <- na_count + sum(is.na(beta_plate_filtered))
    rm(detP);rm(detP_plate);invisible(gc())
    
    # Imputation
    r <- rowSums(is.na(beta_plate_filtered))
    if(sum(r > 0.8*ncol(beta_plate_filtered))>0){
      beta_plate_filtered <- beta_plate_filtered[(r < 0.8*ncol(beta_plate_filtered))>0, ]
      cat('!! Probe failure in more than 80% of samples for',
          sum(r > 0.8*ncol(beta_plate_filtered)),
          'probes on plate',p,'(removed from beta matrix)\n')
    }
    
    out <- capture.output(beta_plate_filtered_imputed <- impute.knn(beta_plate_filtered,k=10,rowmax=0.8)$data)
    rm(beta_plate_filtered);invisible(gc())
    
    
    if(is.null(beta_merged)){
      beta_merged <- beta_plate_filtered_imputed
      rm(beta_plate_filtered_imputed);invisible(gc())
      
    } else {
      beta_merged <- merge(beta_merged, beta_plate_filtered_imputed,
                           by.x='row.names', by.y='row.names', sort=FALSE)
      rm(beta_plate_filtered_imputed);invisible(gc())
      
      # reasign the row names
      row.names(beta_merged) <- beta_merged$Row.names
      
      # remove this column created by the merge function
      beta_merged$Row.names <- NULL
    }
    
    cat(paste('Merged plate ',p,'\n',sep=''))
  }
  cat('\n')
  
  
  cat(paste(na_count,
            ' (',
            round(na_count/(ncol(beta_merged) * nrow(beta_merged)),digits=4),
            '%) ',
            'data points failed detection p-value test\n',sep=''))
  
  
  cat('\nMerged beta matrix has',
      nrow(beta_merged),
      'CpGs and',
      ncol(beta_merged),
      'samples\n\n')
  
  # Save output
  cat('Begin save outputs...')
  save(beta_merged, file=paste(output,'/beta_merged.Rdata',sep=''))
  cat(' done\n\n')
  
  # Remove intermediate files
  cat('Delete intermediate files... ')
  files <- list.files(output, full.names = TRUE)
  ind <- grepl("filtered_norm", files)
  file.remove(files[ind])
  cat(' done\n\n')
  
  
  # Load rho
  if(length(plates) == 1){
    load(paste(log,'/',plates,'_rho.Rdata',sep=''))
  } else {
    for(c in 1:length(plates)){
      load(paste(log,'/',plates[c],'_rho.Rdata',sep=''))
      
      if(c == 1){
        rho_tmp <- rho
      } else {
        rho_tmp <- rbind(rho_tmp, rho)
      }
    }
    rho <- rho_tmp
  }
  
  if(!exists("pheno")){
    pheno <- data.frame(matrix(nrow = ncol(beta_merged),
                               ncol = 3))
    colnames(pheno) <- c("basename", "sentrix_id", "sentrix_pos")
    pheno$basename <- colnames(beta_merged)
    pheno$sentrix_id <- stringr::str_split(pheno$basename, "_", simplify = T)[,1]
    pheno$sentrix_pos <- stringr::str_split(pheno$basename, "_", simplify = T)[,2]
  }
  
  pheno$rho <- numeric(length = nrow(pheno))
  pheno$rho <- rho$rho[match(pheno$basename, rownames(rho))]
  
  # Load controls
  if(length(plates) == 1){
    load(paste(log,'/',plates,'_ctrl.Rdata',sep=''))
  } else {
    for(c in 1:length(plates)){
      load(paste(log,'/',plates[c],'_ctrl.Rdata',sep=''))
      
      if(c ==1){
        ctrl_tmp_g <- ctrl_g
        ctrl_tmp_r <- ctrl_r
      } else {
        ctrl_tmp_g <- cbind(ctrl_tmp_g, ctrl_g)
        ctrl_tmp_r <- cbind(ctrl_tmp_r, ctrl_r)
      }
    }
    
    ctrl_g <- ctrl_tmp_g
    ctrl_r <- ctrl_tmp_r
  }
  
  
  # Load log
  if(length(plates) == 1){
    load(paste(log,'/',plates,'_log_data.Rdata',sep=''))
  } else {
    for(p in 1:length(plates)){
      load(paste(log,'/',plates[p],'_log_data.Rdata',sep=''))
      
      if(p == 1){
        log_data_tmp <- log_data
      } else {
        
        for(x in 1:length(log_data_tmp)){
          log_data_tmp[[x]] <- c(log_data_tmp[[x]], log_data[[x]])
        }
      }
      
      log_data <- log_data_tmp
    }
  }
  
  if(length(plates) != 1){
  save(log_data, file = paste0(log, "/merged_log_data.Rdata"))
  save(rho, file = paste0(log, "/merged_rho.Rdata"))
  }
  
  # Create a report
  cat('Creating RMarkdown Report...\n')
  file.copy(from = paste0(system.file("rmd", "_site.yml", package = "eutopsQC")),
            to = report)
  rmarkdown::render(input = paste0(system.file("rmd", "index.Rmd",
                                               package = "eutopsQC")), 
                    output_file = paste0(report, "index.html", sep = ""))
  rmarkdown::render(input = paste0(system.file("rmd", "1-plate-summary.Rmd",
                                               package = "eutopsQC")), 
                    output_file = paste0(report, "1-plate-summary.html", sep = ""))
  rmarkdown::render(input = paste0(system.file("rmd", "2-control-plots.Rmd", package = "eutopsQC")), 
                    output_file = paste0(report, "2-control-plots.html", sep = ""))
  
  rmarkdown::render(input = paste0(system.file("rmd", "3-qc.Rmd", package = "eutopsQC")),
                    output_file = paste0(report, "3-qc.html", sep = ""))
  rmarkdown::render(input = paste0(system.file("rmd", "4-beta-distributions.Rmd", package = "eutopsQC")),
                    output_file = paste0(report, "4-beta-distributions.html", sep = ""))
  rmarkdown::render(input = paste0(system.file("rmd", "5-snr.Rmd", package = "eutopsQC")),
                    output_file = paste0(report, "5-snr.html", sep = ""))
  rmarkdown::render(input = paste0(system.file("rmd", "6-age-ic-smk.Rmd", package = "eutopsQC")),
                    output_file = paste0(report, "6-age-ic-smk.html", sep = ""))
  if(exists("pheno")){
    rmarkdown::render(input = paste0(system.file("rmd", "7-dimensred.Rmd", package = "eutopsQC")),
                      output_file = paste0(report, "7-dimensred.html", sep = ""))
  }
  
  
  
  # Save snp
  if(!grepl("mouse", array, ignore.case = T) && save.rs==T){
    if(length(plates) == 1){
      load(paste(log,'/',plates,'_rs.Rdata',sep=''))
    } else {
      for(p in 1:length(plates)){
        load(paste(log,'/',plates[p],'_rs.Rdata',sep=''))
        
        if(p == 1){
          rs_tmp <- rs
        } else {
          rs_tmp <- cbind(rs_tmp, rs)
        }
        
        rs <- rs_tmp
      }
    }
  }
  
  # Save merged files
  if(length(plates) != 1){
    if(!grepl("mouse", array, ignore.case = T) && save.rs == T){
      save(rs, file = paste0(log, "/merged_rs.Rdata"))
    }
  }
  
  if(exists("pheno")){
    save(pheno, file = paste0(log, "/pheno_qc.Rdata"))
  }
  
  
  cat('Session info:\n\n')
  print(sessionInfo())
  sink()
  
}
