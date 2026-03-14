#' generate HTML Report
#'
#' @param report report directory
#' @param array array version (EPIC, mouse, v2, or other)
#' @param pheno pheno file (if exists)
#' @param beta_merged beta file
#'
#' @returns html report
#' @export
generateReport <- function(report = reportDir,
                           array = arrayInput,
                           log = logDir,
                           pheno = NULL,
                           beta_merged = betaInput){

  require(EpiDISH)

  # Normalize report path for rmarkdown render, otherwise fails
  report_dir <- normalizePath(report, mustWork = FALSE)

  # Helper function
  render_report <- function(rmd, output_file) {
    rmarkdown::render(
      input         = system.file("rmd", rmd, package = "eutopsQC"),
      output_file   = output_file,
      output_dir    = report_dir,
      knit_root_dir = getwd()
    )
  }

  # Loading log files
  f <- list.files(log, full.names = T, pattern = "*.Rdata")
  for (i in f){
    load(i)
  }

  # Get plate name
  ctrlnames <- f[grepl("_ctrl", f)]
  plates <- unique(gsub("_ctrl.Rdata", "", basename(ctrlnames)))

  cat('Creating RMarkdown Report...\n')
  file.copy(from = system.file("rmd", "_site.yml", package = "eutopsQC"),
            to = report_dir)
  render_report("index.Rmd",                 "index.html")
  render_report("1-plate-summary.Rmd",       "1-plate-summary.html")
  render_report("2-control-plots.Rmd",       "2-control-plots.html")
  render_report("3-qc.Rmd",                  "3-qc.html")
  render_report("4-beta-distributions.Rmd",  "4-beta-distributions.html")
  render_report("5-snr.Rmd",                 "5-snr.html")

  if (!grepl("v2", array, ignore.case = T)){
    # note age, ic and smk indices calculated from beta_merged

    if(exists("pheno") & array != "mouse"){
      out <- epidish(beta.m = beta_merged,
                     ref.m = centEpiFibIC.m,
                     method = "RPC")$estF
      ind <- match(pheno$basename, rownames(out))
      pheno$ic <- out[ind,3]
    }

    render_report("6-age-ic-smk.Rmd", "6-age-ic-smk.html")

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
      beta_merged <- beta_merged_compatible
      render_report("6-age-ic-smk.Rmd", "6-age-ic-smk.html")
    })

  }

  if(exists("pheno")){
    render_report("7-dimensred.Rmd", "7-dimensred.html")
    save(pheno, file = paste0(report_dir, "/pheno_postQC.Rdata"))
  }

  }
