# ---
# title: National Models 5.0 - covariate interpretation
# author: Mannfred Boehm
# created: May 22, 2024
# ---


#NOTES################################

# This script extracts covariate contributions to model predictions. 

# The outputs will be analysed for identifying covariates of importance across several metrics (e.g. BCR, ecology, etc.)

# The first part of this script generates the exported data for the R package function. 

# The second part is the function that summarises the v.5.0 models. 






#WRITE FUNCTION#######################
#'
#'@param covariate_data An exported `data.frame` (see: `data(bam_covariate_importance)` where rows are covariates and columns denote the relative influence of for a given bootstrap replicate by species by BCR permutation. 
#'
#'@param group_by `<tidy-select>` An unquoted expression specifying the grouping variable(s) by which to summarise the relative influence of model covariates. E.g. `bcr` summarises the relative influence of model covariates by Bird Conservation Region, while `species` summarises covariates by taxon. Expressions must be columns that exist in `bam_covariate_importance` or one of the trait databases (see: `traits` argument).  
#'
#'@param species The species to be summarised. Default is `"all"` but can also be a `character` with names (common, scientific, or FLBCs) denoting the species of interest. We follow the taxonomy of the AOS 64th Supplement (July 21, 2023). 
#' Perhaps we can work a list of available species into an exported data.frame, e.g. via `usethis::use_data(unique(sample_id$sci_name))`. 
#' 
#'@param bcr The Bird Conservation Regions to be summarised. Works in the same manner as `species`.
#' 
#'@param traits One of `avonet` (Tobias et al 2022) or `acad`. Can also be a `data.frame` with species as rows (see `species` argument above) and traits as columns (can we handle continuous traits?). 
#' 
#'@param plot If `TRUE`, creates a stacked bar plot with relative influence.
#' 
#'@param colours A `character` of hex codes specifying the colours if `plot = TRUE`. 
#'
#'@param export If `TRUE`, exports to an object the dataframe underlying any plots created by this function.
#' 
#'@return ...tbd
#'
#'@examples ...tbd


bam_explore <- function(data = covs_all, ..., species=c("all", ...), traits = NULL, plot = FALSE){

  
  # Filter the dataset based on the species if specified
  if (species != "all") {
    covs_all <- covs_all |> dplyr::filter(species == !!species)
  }
  
  if (rlang::enquos(...) %in% colnames(.data)) 
  
  
  # sum relative influence by user-specified variable and variable class
  rel_inf_sum <- 
    covs_all |> 
    group_by(..., var_class) |> 
    summarise(sum_influence = sum(rel.inf))
  
}
  
  # sum rel. influence of variable classes 
  # var_sum <- 
  #   covs_all |> 
  #   group_by(var_class)  |>  
  #   summarise(sum_var_class = sum(rel.inf))
  
  
  # sum rel. influence of summarising_var
  var_sum <-
    covs_all |>  
    group_by(summarising_var) |>  
    summarise(sum_bcr = sum(rel.inf))
  
  
  proportion_inf <- 
    rel_inf_sum |> 
    left_join(x = _, bcr_sum, by="bcr") |> 
    mutate(prop = sum_influence/sum_bcr ) 
  
  
  # colour blind palette from: https://colorbrewer2.org/#type=qualitative&scheme=Set3&n=12
  
  cbPalette <- c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5","#d9d9d9","#bc80bd","#ccebc5","#ffed6f")
  
  ggplot() +
    geom_bar(aes(x=summarising_var, y=prop, fill=var_class), data=proportion_inf, stat = "identity") + 
    scale_fill_manual(values = cbPalette) +
    theme_classic()
    # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.title =element_blank()) +
    # theme_bw() 
  
  print(paste("plotting proportion of variable influence per"), summarising_var, sep=" ")
  
}
                        


# root <- "G:/Shared drives/BAM_NationalModels/NationalModels5.0/output/bootstraps"
# 
# prediction_files <- list.files(root)
# 
# loop <- 
#   prediction_files %>% 
#   stringr::str_split_fixed(pattern="_", n=3) %>% 
#   dplyr::as_tibble() %>% 
#   magrittr::set_colnames(c("spp", "bcr", "boot")) %>% 
#   dplyr::arrange(spp, bcr)
# 
# bcr.i <- loop$bcr[i]
# spp.i <- loop$spp[i]
# boot.i <- loop$boot[i]




#4. create nested dataframe ----

# Main Dataframe
# main_df <- data.frame(
#   species = unique(covs_all$sci_name),
#   common_name = unique(covs_all$common_name)
# )
# 
# # Species Dataframes with nested BCR and Trait Dataframes
# species1_df <- data.frame(
#   Habitat = "Habitat 1",
#   Population = 100,
#   BCRData = I(list(data.frame(BCRValue = 0.8))),
#   TraitData = I(list(data.frame(Trait1 = "Trait1 Value 1", Trait2 = "Trait2 Value 1")))
# )

