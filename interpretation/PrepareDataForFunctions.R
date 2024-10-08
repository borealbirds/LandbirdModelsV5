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
gbm_objs <- list.files(file.path(root, "output", "bootstraps"))[1:5]


# import extraction lookup table to obtain covariate classes (`var_class`)
# lookup table is missing "Year" and "Method", so manually adding here
varclass_lookup <- readxl::read_xlsx(file.path(root, "NationalModels_V5_VariableList.xlsx"), sheet="ExtractionLookup") |>
  dplyr::select(Category, Label) |> 
  dplyr::rename(var_class = Category, var = Label) |> 
  tibble::add_row(var_class = "Method", var = "method") |> 
  tibble::add_row(var_class = "Year", var = "year")


# import IBP nomenclature for appending scientific names to FLBCs
ibp <- read_csv(file.path(root, "data", "Extras", "sandbox_data", "trait_data_for_summarising_covariates", "institute_for_bird_populations_species_codes.csv")) 


# create an index from `gbm_objs` containing the species (FLBC), BCR, and bootstrap replicate
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
saveRDS(gbm_objs,  file="C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/gbm_objs.rds")


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
      # we discard the grid (too much data) and keep the mean and std. dev. of the response for covariates i,j
      grid_ij <- plot.gbm(x = b.i, return.grid = TRUE, i.var = c(i, j), continuous.resolution=25, type="response")  
      
      pts[[interaction_index]] <- matrix(data=c(mean(grid_ij$y), sd(grid_ij$y)), ncol=2, nrow=1)
      colnames(pts[[interaction_index]]) <- c("y_mean", "y_sd")
      
      names(pts)[interaction_index] <- paste(b.i$var.names[i], b.i$var.names[j], sep = ".") # label the interaction
       
      interaction_index <- interaction_index + 1
      
    }
  }
  
  boot_pts_i2[[q]] <- pts
  
  # print progress
  cat(paste("\riteration", q))
  Sys.sleep(0.001)
}

saveRDS(boot_pts_i2, file="C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/boot_pts_i2.rds")

# create an index of every bcr x common_name x 2-way interactions 
# `RcppAlgo::comboGrid` is like `expand.grid` and `tidyr::crossing` but avoids duplicates 
# e.g. for our purposes (var1=x, var2=y) is a duplicate of (var1=y, var2=x) 
# so the Cartesian Product of comboGrid is smaller
boot_group_keys_i2 <- 
  RcppAlgos::comboGrid(unique(bam_covariate_importance$var), 
            unique(bam_covariate_importance$var), 
            bam_covariate_importance$bcr, 
            bam_covariate_importance$common_name, repetition =  FALSE) |> 
  tibble::as_tibble() |> 
  dplyr::rename(var_1=Var1, var_2=Var2, bcr=Var3, common_name=Var4)



# for every zth species x bcr x 2-way interaction tuple (rows in `boot_group_keys_i2`):
# gather the relevant bootstrap predictions from `boot_pts_i2`
# (each element of `boot_pts_i2` has many `y_means`; one for every 2-way interaction)
# and enter each wth `y_mean` (that matches the current 2-way interaction of interest) 
# as a sub-element of [[z]]
# output: a list of lists where the top level elements are species x bcr x 2-way interaction permutations, 
# and the second-level elements are matrices of the associated average and sd of the bootstrapped prediction spaces

boot_pts_sorted_i2 <- list()
for (z in 1:nrow(boot_group_keys_i2)){
  
  # zth bcr x species x 2-way interaction permutation
  key_z <- boot_group_keys_i2[z,]
  
  # problem: which elements in boot_pts_i2 match the species x bcr tuple defined in key_z?
  # approach: `sample_id` and `boot_pts_i2` have the same length and order (they both are derived from gbm_objs)
  # so we can identify elements in `boot_pts_i2` using info from `sample_id`
  # first, identify a unique species x bcr tuple, then gather covariate interactions into it
  boot_pts_index <- which(sample_id$bcr == key_z$bcr & sample_id$common_name == key_z$common_name)
  
  # a list of bootstrap predictions for zth species x bcr permutation 
  spp_bcr_list <- boot_pts_i2[boot_pts_index]
  
  # gather all relevant 2-way interaction bootstrap averages for the zth species x bcr permutation
  spp_bcr_var_i2 <- list()
  for (w in 1:length(spp_bcr_list)) {
    
    # get covariate names for the current spp x bcr permutation 
    # NOTE: this is admittedly wasteful because we're only interested in the interactions within `key_z`
    i2_names <- names(spp_bcr_list[[w]])
    
    # find which sub-element of spp_bcr_list[[w]] match the current 2-way interaction of interest
    # account for the possibility of being indexed as x*y or y*x
    if (paste(key_z$var_1, key_z$var_2, sep=".") %in% i2_names |
        paste(key_z$var_2, key_z$var_1, sep=".") %in% i2_names){
      
      matching_index <- 
        which(i2_names == paste(key_z$var_1, key_z$var_2, sep=".") | 
              i2_names == paste(key_z$var_2, key_z$var_1, sep="."))
    
      # assign current bcr x species x 2-way interaction x bootstrap tuple as a top-level element of `spp_bcr_var_i2[[w]]`
      # these bootstraps will eventually be gathered under a single bcr x species x 2-way interaction element in `boot_pts_sorted_i2`
      spp_bcr_var_i2[w] <- spp_bcr_list[[w]][matching_index]
      names(spp_bcr_var_i2)[w] <- paste("bootstrap replicate", w, sep="_")
    
    } #close if()
  
  } #close nested loop
  
  # take mean across bootstraps of the zth bcr x spp x 2-way interaction tuple
  if (purrr::is_empty(spp_bcr_var_i2) == FALSE){
    mean_z <- 
      spp_bcr_var_i2 |> 
      purrr::list_c() |> 
      colMeans() |> 
      matrix(data=_, ncol=2) |> 
      data.frame() |> 
      setNames(c("mean_y_mean", "mean_y_sd"))
    
    if (mean_z$mean_y_mean >= 0.10){
      boot_pts_sorted_i2[[z]] <- mean_z
    } else {
      boot_pts_sorted_i2[[z]] <- list()
    } #close nested else
    
  } # close top level if()
  else {
    boot_pts_sorted_i2[[z]] <- list()
  } # close top level else
  
  names(boot_pts_sorted_i2)[z] <- paste(key_z$bcr, key_z$common_name, key_z$var_1, key_z$var_2, sep=".")
  
  # print progress
  cat(paste("\riteration", z, "of", nrow(boot_group_keys_i2)))
  Sys.sleep(0.000001)
  
} # close top level loop


# remove empty elements (output reduced from 8.9MB to 905.2kB;with threshold set and non-threshold elements removed)
boot_pts_reduced_i2 <- purrr::discard(boot_pts_sorted_i2, purrr::is_empty)
saveRDS(boot_pts_reduced_i2, file="C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/boot_pts_reduced_i2.rds")








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