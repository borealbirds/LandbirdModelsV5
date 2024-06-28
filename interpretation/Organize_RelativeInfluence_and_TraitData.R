# ---
# title: National Models 5.0 - data organization for analysing relative influence of covariates and species traits
# author: Mannfred Boehm
# created: June 28, 2024
# ---

#NOTES################################

# This script extracts covariate contributions to model predictions, and also synthesizes various trait databases.

# The outputs will be analysed for identifying covariates of importance across several metrics (e.g. BCR, ecology, etc.).

# This script generates the exported data for the R package function. 



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






#3. import trait databases----


avonet <- readxl::read_excel(file.path(root, "data", "Extras", "sandbox_data", "trait_data_for_summarising_covariates", "avonet_database_eBird_taxonomy.xlsx"), sheet="AVONET2_eBird")



# check that all BAM species are represented in trait databases 
all(sample_id$sci_name %in% avonet$Species2)




