# ---
# title: National Models 5.0 - data organization for analysing relative influence of covariates and species traits
# author: Mannfred Boehm
# created: June 28, 2024
# ---

#NOTES################################

# This script extracts covariate contributions to model predictions, and also synthesizes various trait databases.

# The outputs will be analysed for identifying covariates of importance across several metrics (e.g. BCR, ecology, etc.).

# This script generates the exported data for R package functions that summarise covariate performance. 



#PREAMBLE############################

#1. Attach packages----
library(gbm)
library(tidyverse)



#2. Create a list of dataframes containing relative influence per covariate----
#   Every list element represents a bootstrap sample 

# connect to BAM Drive and find bootstrap files 
root <- "G:/Shared drives/BAM_NationalModels5"

gbm_objs <- list.files(file.path(root, "output", "bootstraps"))[sample(1:500, 250)]


# import extraction lookup table to obtain covariate classes (`var_class`)
# lookup table is missing "Year" and "Method", so manually adding here
varclass_lookup <- readxl::read_xlsx(file.path(root, "NationalModels_V5_VariableList.xlsx"), sheet="ExtractionLookup") |>
  dplyr::select(Category, Label) |> 
  dplyr::rename(var_class = Category, var = Label) |> 
  tibble::add_row(var_class = "Method", var = "method") |> 
  tibble::add_row(var_class = "Year", var = "year")


# import IBP nomenclature for appending scientific names to FLBCs
ibp <- read_csv(file.path(root, "data", "Extras", "sandbox_data", "trait_data_for_summarising_covariates", "institute_for_bird_populations_species_codes.csv")) 


# create an index containing the species (FLBC), BCR, and bootstrap replicate
# append binomial names to FLBC
sample_id <- 
  gbm_objs |> 
  stringr::str_split_fixed(pattern="_", n=3) |> 
  gsub("\\.R", "", x = _) |>
  tibble::as_tibble() |> 
  magrittr::set_colnames(c("spp", "bcr", "boot")) |> 
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
    dplyr::cross_join(x = _, sample_id[i,]) |> 
    dplyr::mutate(file_name = gbm_objs[i])
  
  # print progress
  cat(paste("\riteration", i))
  Sys.sleep(0.001)
}

# mergelist of dataframes
bam_covariate_importance <- suppressMessages(purrr::reduce(covs, full_join))





#3. import trait databases----

# AVONET: https://doi.org/10.1111/ele.13898
avonet <- readxl::read_excel(file.path(root, "data", "Extras", "sandbox_data", "trait_data_for_summarising_covariates", "avonet_database_eBird_taxonomy.xlsx"), sheet="AVONET2_eBird")
pif


# ACAD Scores
# Population Size (`PS-g`) 
# Breeding and Non-breeding Distributions (`BD-g` and `ND-g`)
# Threats to Breeding (`TB-c`) and Non-breeding (`TN-c`)
# Population Trend (`PT-c`)
acad <- read_csv(file.path(root, "data", "Extras", "sandbox_data", "trait_data_for_summarising_covariates", "acad_global_2024.csv"))


# check that all BAM species are represented in trait databases 
all(sample_id$sci_name %in% avonet$Species2)
all(sample_id$sci_name %in% acad$sci_name)




#4. Create list of bcr x species x covariate permutations (with bootstraps)----
# for partial dependence plotting  


# this loop creates a list of lists of `data.frame`s. 
# each top level element is a bootstrap replicate, each second-level element represents a covariate for that bootstrap.
# each covariate has a corresponding `data.frame`, with column 1 being the covariate domain and column 2 being the predicted range.
boot_pts <- list()
for (q in 1:length(gbm_objs)){
  
  # load a bootstrap replicate 
  load(file.path(root, "output", "bootstraps", gbm_objs[q]))
  
  # find evaluation points (`data.frame`) for every covariate (indexed by `i`)  
  pts <- list()
  for (i in 1:length(b.i$var.names)){
    pts[[i]] <- plot.gbm(x=b.i, return.grid = TRUE, i.var = i)
  }
  
  boot_pts[[q]] <- pts
  
  # print progress
  cat(paste("\riteration", q))
  Sys.sleep(0.001)
}



# group `b.i` objects by bcr x species x var (or var_class?)
# group_keys() finds the names of the group permutations
# 0000 permutations of  bcr x common_name x var 
boot_group_keys <- 
  bam_covariate_importance |> 
  dplyr::group_by(bcr, common_name, var) |> 
  dplyr::group_keys()



# for every zth species x bcr x covariate permutation:
# gather the relevant bootstrap predictions from `boot_pts`
# and enter each wth dataframe as an element of [[z]]
# output: a list of lists where the top level elements are species x bcr x var permutations, 
# and the second-level elements are dataframes of the associated bootstrapped model predictions

boot_pts_sorted <- list()
for (z in 1:nrow(boot_group_keys)){
  
  # zth covariate
  key_z <- boot_group_keys[z,]
  
  # `sample_id` and `boot_pts` have the same length and order
  # so we can annotate elements in `boot_pts` using info from `sample_id`
  boot_pts_index <- which(sample_id$bcr == key_z$bcr & sample_id$common_name == key_z$common_name)
  
  # a list of bootstrap predictions for zth species x bcr permutation 
  spp_bcr_list <- boot_pts[boot_pts_index]
  
  # search for all bootstrap predictions for the zth covariate
  spp_bcr_var <- list()
  for (w in 1:length(spp_bcr_list)) {
    
    # get covariate names for the current spp x bcr permutation 
    var_names <- 
      lapply(spp_bcr_list[[w]], colnames) |> 
      lapply(X=_, `[[`, 1) |> #don't need the name of the y variable
      purrr::flatten_chr()
    
    # find which elements match the current covariate of interest
    spp_bcr_var[w] <- spp_bcr_list[[w]][which(var_names %in% key_z)] 
    names(spp_bcr_var)[w] <- paste("bootstrap replicate", w, sep="_")
  }
  
  boot_pts_sorted[[z]] <- spp_bcr_var
  names(boot_pts_sorted)[z] <- paste(boot_group_keys[z, "bcr"], boot_group_keys[z, "common_name"], boot_group_keys[z, "var"], sep = "_")
  
  # print progress
  cat(paste("\riteration", z))
  Sys.sleep(0.001)
}

