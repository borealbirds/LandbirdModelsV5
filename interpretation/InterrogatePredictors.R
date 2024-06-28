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

#PREAMBLE############################

#1. Attach packages----
library(gbm)
library(tidyverse)



#2. Create a list of dataframes containing relative influence per covariate----
#   Every list element represents a bootstrap sample 

# connect to BAM Drive and find bootstrap files 
root <- "G:/Shared drives/BAM_NationalModels5"

gbm_objs <- list.files(file.path(root, "output", "bootstraps"))[1:3]


# import extraction lookup table to obtain covariate classes (`var_class`)
# lookup table is missing "Year" and "Method", so manually adding here
varclass_lookup <- readxl::read_xlsx(file.path(root, "NationalModels_V5_VariableList.xlsx"), sheet="ExtractionLookup") |>
  dplyr::select(Category, Label) |> 
  dplyr::rename(var_class = Category, var = Label) |> 
  tibble::add_row(var_class = "Method", var = "method") |> 
  tibble::add_row(var_class = "Year", var = "year")


# import Institute for Bird Populations database for appending scientific names to FLBCs
ibp <- read_csv(file.path(root, "data", "Extras", "sandbox_data", "trait_data_for_summarising_covariates", "institute_for_bird_populations_species_codes.csv")) |> 
  dplyr::select(SPEC, SCINAME, COMMONNAME) |> 
  dplyr::rename(spp = SPEC, sci_name = SCINAME, common_name = COMMONNAME)


# create an index containing the species (FLBC), BCR, and bootstrap replicate
# append binomial names to FLBC
sample_id <- 
  gbm_objs |> 
  stringr::str_split_fixed(pattern="_", n=3) |> 
  gsub("\\.R", "", x = _) |>
  tibble::as_tibble() |> 
  magrittr::set_colnames(c("spp", "bcr", "boot")) |> 
  dplyr::arrange(spp, bcr) |> 
  dplyr::left_join(ibp)




#3. import trait databases----


avonet <- readxl::read_excel(file.path(root, "data", "Extras", "sandbox_data", "trait_data_for_summarising_covariates", "avonet_database_eBird_taxonomy.xlsx"), sheet="AVONET2_eBird")



# check that all BAM species are represented in trait databases 
all(sample_id$sci_name %in% avonet$Species2)




# create a list of dataframes containing relative influence per covariate
# every element of this list comes from a single bootstrap
covs <- list()
for(i in 1:length(gbm_objs)){
  
   # loads a `gbm` object named `b.i`
  load(file=file.path(root, "output", "bootstraps", gbm_objs[i]))
  
  # `summary(b.i)` is a `data.frame` with columns `var` and `rel.inf`
  # `cross_join()` attaches the spp/bcr/boot info to every covariate of a given iteration i 
  covs[[i]] <- b.i |>
    gbm::summary.gbm(plotit = FALSE) |>
    dplyr::left_join(varclass_lookup, by="var") |>
    as_tibble() |> 
    dplyr::cross_join(x = _, sample_id[i,]) 
   
  # print progress
  cat(paste("\riteration", i))
    Sys.sleep(0.05)
}

# flatten list of dataframes
covs_all <- purrr::reduce(covs, full_join)




#4. create nested dataframe ----

# Main Dataframe
main_df <- data.frame(
  species = unique(covs_all$sci_name),
  common_name = unique(covs_all$common_name)
)

# Species Dataframes with nested BCR and Trait Dataframes
species1_df <- data.frame(
  Habitat = "Habitat 1",
  Population = 100,
  BCRData = I(list(data.frame(BCRValue = 0.8))),
  TraitData = I(list(data.frame(Trait1 = "Trait1 Value 1", Trait2 = "Trait2 Value 1")))
)




#WRITE FUNCTION#######################
#'
#'@param .data An exported `data.frame` (see: `data(xyz)` where rows are covariates and columns denote the relative influence of for a given bootstrap replicate by species by BCR permutation. 
#'
#'@param ... `<tidy-select>` An unquoted expression specifying the grouping variable by which to summarise the relative influence of model covariates. E.g. `bcr` summarises the relative influence of model covariates by Bird Conservation Region, while `species` summarises covariates by taxon.  
#'
#' @param species The species to be summarised. Default is `"all"` but can also be a `character` with the FLBC or scientific name denoting the species of interest. We follow the taxonomy of the AOS 64th Supplement (July 21, 2023). 
#' Perhaps we can work a list of available species into an exported data.frame, e.g. via `usethis::use_data(unique(sample_id$sci_name))`. 
#' 
#' @param traits One of `avonet` (Tobias et al 2022), `x`, or `y`. Can also be a `data.frame` with species as rows (see `species` argument above) and traits as columns (can we handle continuous traits?). 
#' 
#' @param plot If `TRUE`, creates a stacked bar plot with relative influence.
#' 
#' @param colours A `character` of hex codes specifying the colours if `plot = TRUE`. 
#'
#' @return a `data.frame` To be passed to a plotting function...tbd
#'
#' @examples ...tbd

# what a user might want:
# rel.inf by species
# rel.inf by bcr
# to be able to choose what covariates they want to compare (i.e. exclude some covariates)? Or does that misrepresent the models somehow?


bam_relative_influence <- function(.data = covs_all, ..., species=c("all", ...), traits = NULL, plot = FALSE){

  
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




