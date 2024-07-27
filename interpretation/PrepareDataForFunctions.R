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
library(RcppAlgos)
library(tidyverse)



#2. Create a list of dataframes containing relative influence per covariate----
#   Every list element represents a bootstrap sample 

# connect to BAM Drive and find bootstrap files 
root <- "G:/Shared drives/BAM_NationalModels5"

# gbm objects stored on the BAM drive
# 129.2 GB for 10 bootstraps
gbm_objs <- list.files(file.path(root, "output", "bootstraps"))[sample(1:300, 100)]


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





#3. Create list of bcr x species x covariate permutations (with bootstraps)----
# for partial dependence plotting  


# this loop creates a nested list of `data.frame`s. 
# each top level element a bcr x species x bootstrap tuple, with second-level elements as the rth covariate for that bootstrap.
boot_pts <- list()
for (q in 1:length(gbm_objs)){
  
  # load a bootstrap replicate 
  load(file.path(root, "output", "bootstraps", gbm_objs[q]))
  
  # find evaluation points (`data.frame`) for every covariate (indexed by `r`)
  # this can sped up by reducing `continuous.resolution=100` in `plot.gbm()`
  pts <- list()
  for (r in 1:length(b.i$var.names)){
    pts[[r]] <- plot.gbm(x=b.i, return.grid = TRUE, i.var = r, type="response")
  }
  
  boot_pts[[q]] <- pts
  
  # print progress
  cat(paste("\riteration", q))
  Sys.sleep(0.001)
}



# from the family of sets (BCRs, species, covariates) 
# group_keys() finds the names of the unique tuples
# There are 00000 bcr x common_name x var tuples 
boot_group_keys <- 
  bam_covariate_importance |> 
  dplyr::group_by(bcr, common_name, var) |> 
  dplyr::group_keys()



# for every zth bcr x species x var tuple (rows in `boot_group_keys`):
# search for the relevant bootstrap predictions in `boot_pts` by matching species and bcr at the top level of `boot_pts` then matching `var` at the second level
# then, enter bootstrap dataframes for a given bcr x species x var tuple as a sub-element of [[z]]
# output: a list of lists where the top level elements are species x bcr x var tuples, 
# and the second-level elements are dataframes of that tuple's bootstrapped model predictions
# you can query every zth tuple for its bcr x species x var identity by: names(boot_pts_sorted)[z] 
boot_pts_sorted <- list()
for (z in 1:nrow(boot_group_keys)){
  
  # zth bcr x species x var tuple
  key_z <- boot_group_keys[z,]
  
  # `sample_id` and `boot_pts` have the same length and order (they are derived from `gbm_objs`)
  # so we can identify elements in the latter using info from the former
  # search for all bootstraps for a given bcr x sample tuple
  boot_pts_index <- which(sample_id$bcr == key_z$bcr & sample_id$common_name == key_z$common_name)
  
  # a list of bootstrap predictions for species x bcr tuple
  spp_bcr_list <- boot_pts[boot_pts_index]
  
  # for a given bcr x species tuple, search within every bootstrap model for
  spp_bcr_var <- list()
  for (w in 1:length(spp_bcr_list)) {
    
    # get all possible covariate names for the current bcr x species bootstrap
    var_names <- 
      lapply(spp_bcr_list[[w]], colnames) |> 
      lapply(X=_, `[[`, 1) |> #don't need the name of the y variable
      purrr::flatten_chr()
    
    # search the wth bootstrap to find the current covariate of interest
    spp_bcr_var[w] <- spp_bcr_list[[w]][which(var_names %in% key_z$var)]
    names(spp_bcr_var)[w] <- paste("bootstrap replicate", w, sep="_")
  }
  
  boot_pts_sorted[[z]] <- spp_bcr_var
  names(boot_pts_sorted)[z] <- paste(boot_group_keys[z, "bcr"], boot_group_keys[z, "common_name"], boot_group_keys[z, "var"], sep = "_")
  
  # print progress
  cat(paste("\riteration", z))
  Sys.sleep(0.001)
}



#4. Create a list of two-way covariate interactions by bcr x species x bootstrap replicate----

# this loop creates a list of lists of `data.frame`s. 
# each top level element is a bootstrap replicate, each second-level element represents a 2-way covariate interaction for that bootstrap.
# each covariate interaction has a corresponding `data.frame`, with columns 1 and 2 being the covariate domains and column 3 being the response (density)
# NOTE: In computing interactions involving discrete variables (e.g. `MODISLCC_1km`) the resultant plot will differ compared to interactions between two continuous variables (heatmap).
# `plot.gbm()` produces faceted line plots instead of a heatmap. Each line plot corresponds to one of the 17 levels of `MODISLCC_1km`, showing how `MODISLCC_1km` affects the response variable at each level.


boot_pts_i2 <- list()
for (q in 1:length(gbm_objs)){
  
  # load a bootstrap replicate 
  load(file.path(root, "output", "bootstraps", gbm_objs[q]))
  
  # find evaluation points (`data.frame`) for every covariate permutation of degree 2 (indexed by i,j) 
  pts <- list()
  n <- length(b.i$var.names) # get the number of variables
  interaction_index <- 1 # starts at 1 and increases for every covariate interaction computed. Resets at one when moving to the next bootstrap model.
  
  # end at n-1 to avoid finding the interaction of variable n x variable n
  for (i in 1:(n-1)) {  
    
    # start at i+1 to avoid finding the interaction of variable 1 x variable 1
    for (j in (i+1):n) { 
      
      # `continuous.resolution` is defaulted at 100: with two covariates this produces a dataframe of length=100*100 (10000 rows)
      # by lowering to 25, we get 625 rows, which should still be enough resolution to find local maximums
      pts[[interaction_index]] <- plot.gbm(x = b.i, return.grid = TRUE, i.var = c(i, j), continuous.resolution=25, type="response")
      names(pts)[interaction_index] <- paste(b.i$var.names[i], b.i$var.names[j], sep = "_") # label the interaction
      interaction_index <- interaction_index + 1
    }
  }
  
  boot_pts_i2[[q]] <- pts
  
  # print progress
  cat(paste("\riteration", q))
  Sys.sleep(0.001)
}


data <- plot.gbm(x=b.i, i.var=c("SCANFIBlackSpruce_1km", "year"), type="response", continuous.resolution = 12, return.grid = TRUE)
data2 <- plot.gbm(x=b.i, i.var=c("WetSeason_1km", "year"), type="response", continuous.resolution = 12, return.grid = TRUE)
data3 <- plot.gbm(x=b.i, i.var=c("MODISLCC_1km", "year"), type="response", continuous.resolution = 12, return.grid = TRUE)

# take the mean and SD of the response (y) across the 2-way interaction space. 


# create an index of every bcr x common_name x 2-way interactions 
# there are 0000 permutations of  bcr x common_name x 2-way interactions 
unique_bcr_spp <- dplyr::distinct(bam_covariate_importance, bcr, common_name)

boot_group_keys_i2 <- 
  comboGrid(unique(bam_covariate_importance$var), 
            unique(bam_covariate_importance$var), repetition =  FALSE) |> 
  as_tibble() |> 
  tidyr::crossing(unique_bcr_spp) |> 
  dplyr::rename(var1=Var1, var2=Var2)
  
  

# NEED TO BE CAREFUL THAT the horizontal order of var1 and var2 in boot_group_keys_i2 
# does not affect the loop's ability to find interactions (e.g. the loop searches for var2, var1 and fails when var1, var2 exists)

# for every zth species x bcr x 2-way interaction tuple (rows in `boot_group_keys_i2`):
# gather the relevant bootstrap predictions from `boot_pts_i2`
# and enter each wth dataframe as a sub-element of [[z]]
# output: a list of lists where the top level elements are species x bcr x 2-way interaction permutations, 
# and the second-level elements are dataframes of the associated bootstrapped model predictions

boot_pts_sorted_i2 <- list()
for (z in 1:nrow(boot_group_keys_i2)){
  
  # zth bcr x species x 2-way interaction permutation
  key_z <- boot_group_keys_i2[z,]
  
  # `sample_id` and `boot_pts_i2` have the same length and order 
  # so we can annotate elements in `boot_pts_i2` using info from `sample_id`
  # identify a unique species x bcr combo, then gather covariate interactions into it
  boot_pts_index <- which(sample_id$bcr == key_z$bcr & sample_id$common_name == key_z$common_name)
  
  # a list of bootstrap predictions for zth species x bcr permutation 
  spp_bcr_list <- boot_pts_i2[boot_pts_index]
  
  # search for all bootstrap predictions for the zth species x bcr permutation
  spp_bcr_var <- list()
  for (w in 1:length(spp_bcr_list)) {
    
    # get covariate names for the current spp x bcr permutation 
    i2_names <- 
      lapply(spp_bcr_list[[w]], colnames) |> 
      lapply(X=_, `[[`, 1) |> # don't need the name of the y variable
      purrr::flatten_chr()
    
    # find which elements match the current covariate combination of interest
    matching_indices <- which((i2_names == key_z$var1) | (i2_names == key_z$var2))
    spp_bcr_var[[w]] <- spp_bcr_list[[w]][matching_indices]
    names(spp_bcr_var)[w] <- paste("bootstrap replicate", w, sep="_")
  }
  
  interaction_name <- paste(key_z$bcr, key_z$common_name, paste(key_z$var1, key_z$var2, sep = "_"), sep = "_")
  if (!interaction_name %in% names(boot_pts_sorted_i2)) {
    boot_pts_sorted_i2[[interaction_name]] <- list()
  }
  
  boot_pts_sorted_i2[[interaction_name]][[length(boot_pts_sorted_i2[[interaction_name]]) + 1]] <- spp_bcr_var
  
  # print progress
  cat(paste("\riteration", z))
  Sys.sleep(0.001)
}

# list with 37487 elements
# 1 species x 3 bcrs x 100 bootstraps x 66! 2-way interactions
# 12 var_classes = 12! = 479 million











#0000. import and structure trait datasets----

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