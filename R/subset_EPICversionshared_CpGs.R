
#' subset_versionshared_CpGs
#' 
#' Subset a beta matrix for v1.v2 compatible CpG-probes
#'
#' @param matrix beta_merged
#' @param character "EPIC or "EPICv2"
#' @return beta matrix subset
#' @export
#'
#' @examples
#' subset_versionshared_CpGs(beta_merged, "EPICv2")

subset_versionshared_CpGs <- function(beta, array){
  
  cat("Making EPIC v1.v2 compatible beta-subset, rownames will be EPICv1 IDs")
  
  #load packages if needed
  packages <- c("dplyr","tidyverse", "janitor")
  for (package in packages){
    if (!package %in% .packages()) {
      # Load it
      library(
        package,
        character.only = TRUE
      )
    }
  }
  
  if(grepl("v2", array, ignore.case = T)){

      # keep only v1 compatible probes, and rename to v1 ID
      beta_subset = beta |> 
        as.data.frame() |>
        rownames_to_column(var="IlmnID") |>
        filter(IlmnID %in% c(mapping_compatible_v1v2$IlmnID)) |>
        left_join(mapping_compatible_v1v2) |>
        select(!IlmnID)
      
      # Compute mean across replicate probes with exact sequence match 
      # (non-sequence match already removed from the map)
      
      # Find reps
      dupes <- beta_subset |> janitor::get_dupes(EPICv1_Loci) |> pull(EPICv1_Loci) |> unique()
      
      # Prepare df non reps
      beta_subset_nondupe = beta_subset |> 
        filter(!EPICv1_Loci %in% dupes)
       
      # calculate mean beta and prepare df for reps 
      beta_subset_dupe = beta_subset |> 
        filter(EPICv1_Loci %in% dupes) |>
        group_by(EPICv1_Loci) |>
        dplyr::reframe(across(everything(), ~ mean(.x))) |>
        ungroup() |> distinct() 
      
      # final matrix, with EPICv1 IDs as rowname
      beta_subset = bind_rows(beta_subset_dupe, beta_subset_nondupe) |>
        tibble::column_to_rownames("EPICv1_Loci") |>
        as.matrix()
      
    } else{

      # keep only v2 compatible probes, and rename to v1 ID
      beta_subset = beta |> 
        as.data.frame() |>
        rownames_to_column(var="EPICv1_Loci") |>
        filter(EPICv1_Loci %in% unique(mapping_compatible_v1v2$EPICv1_Loci)) |>
        tibble::column_to_rownames("EPICv1_Loci") |>
        as.matrix()
    }
  
  return(beta_subset)
}

