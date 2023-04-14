#' installBiocDependencies
#'
#' Helper function to install packages when required. List of packages is contained in packageList data
#'
#' @param packageList packages to be installed
#' @return Installs missing packages
#' @export

installBiocDependencies <- function(packageList){

  packageList <- packageList
  if(!all(packageList %in% rownames(installed.packages()))){

    to.install = packageList[!packageList %in% rownames(installed.packages())]
    cat("Note: You have missing Bioconductor packages. These will be installed.\n")
    for (i in to.install){
      suppressMessages({
          cat("Installing package", i)
          BiocManager::install(i, quiet = TRUE, update = F, ask = F)
          cat("... done\n")
        })

    }}

    if (!requireNamespace("WIDclocks", quietly = TRUE)){
      cat("Installing package WIDclocks from github")
      devtools::install_github("chiaraherzog/WIDclocks", quiet = TRUE, upgrade = "never")
      cat("... done\n")
    }

  if (!requireNamespace("ewastools", quietly = TRUE)){
    cat("Installing package ewastools from github")
    devtools::install_github("hhhh5/ewastools", quiet = TRUE, upgrade = "never")
    cat("... done\n")
  }

}

