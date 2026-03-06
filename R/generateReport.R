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
                           pheno = NULL,
                           beta_merged = betaInput){
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

  if (!grepl("v2", array, ignore.case = T)){
    # note age, ic and smk indices calculated from beta_merged

    if(exists("pheno") & array != "mouse"){
      out <- epidish(beta.m = beta_merged,
                     ref.m = centEpiFibIC.m,
                     method = "RPC")$estF
      ind <- match(pheno$basename, rownames(out))
      pheno$ic <- out[ind,3]
    }

    rmarkdown::render(input = paste0(system.file("rmd", "6-age-ic-smk.Rmd", package = "eutopsQC")),
                      output_file = paste0(report, "6-age-ic-smk.html", sep = ""))

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
      rmarkdown::render(input = paste0(system.file("rmd", "6-age-ic-smk.Rmd", package = "eutopsQC")),
                        output_file = paste0(report, "6-age-ic-smk.html", sep = ""))
    })

  }

  if(exists("pheno")){
    rmarkdown::render(input = paste0(system.file("rmd", "7-dimensred.Rmd", package = "eutopsQC")),
                      output_file = paste0(report, "7-dimensred.html", sep = ""))

    save(pheno, file = paste0(report, file = 'pheno_postQC.Rdata'))
  }

  }
