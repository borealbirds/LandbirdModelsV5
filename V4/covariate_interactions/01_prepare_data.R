# ---
# title: National Models 4.0 - estimating and ranking 2-way covariate interactions for CONWA and CAWA
# author: Mannfred Boehm
# created: September 13, 2024
# ---

#NOTES################################

# This script prepares data alongside e.g. conw_bcr60.R-type scripts 
# that are needed for `summarise_interactions()`




#1. attach packages----
print("* attaching packages on master *")
library(gbm)
library(RcppAlgos)
library(tidyverse)


# connect to BAM Drive and find bootstrap files 
root <- "G:/Shared drives/BAM_NationalModels4/NationalModels4.0/Feb2020/out/boot"
dirs <- c("CAWA/BCR_60", "CAWA/BCR_12", "CONW/BCR_60", "CONW/BCR_81")

# gbm objects stored on the BAM drive
# filter for CONW (BCR60, BCR81) and CAWA (BCR60, BCR12)
gbm_objs <- 
  lapply(file.path(root, dirs), list.files, full.names = TRUE) |> 
  unlist()

# create an index from `gbm_objs` containing the species (FLBC), BCR, and bootstrap replicate
sample_id <- 
  gbm_objs |> 
  sub("^.*gnmboot-", "", x = _) |> 
  sub("\\.RData", "", x = _) |>
  stringr::str_split_fixed(pattern="-", n=3) |> 
  tibble::as_tibble() |> 
  magrittr::set_colnames(c("spp", "bcr", "boot"))

# create a list of dataframes containing relative influence per covariate
# every element of this list comes from a single bootstrap
covs <- list()
for(i in 1:length(gbm_objs)){
  
  # loads a `gbm` object named `b.i`
  load(gbm_objs[i])
  
  # `summary(b.i)` is a `data.frame` with columns `var` and `rel.inf`
  # `cross_join()` attaches the spp/bcr/boot info to every covariate of a given iteration i 
  covs[[i]] <- out |>
    gbm::summary.gbm(plotit = FALSE) |>
    as_tibble() |> 
    dplyr::cross_join(x = _, sample_id[i,]) |> 
    dplyr::mutate(file_name = gbm_objs[i])
  
  # print progress
  cat(paste("\riteration", i))
  Sys.sleep(0.00001)
}

covariate_importance <- suppressMessages(purrr::reduce(covs, full_join))
#saveRDS(covariate_importance,  file="C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/boreal_avian_modelling_project/NationalModelsv5/v4/covariate_interactions/02_covariate_importance.rds")




# create an index of every bcr x common_name x 2-way interactions 
# `RcppAlgo::comboGrid` is like `expand.grid` and `tidyr::crossing` but avoids duplicates 
# e.g. for our purposes (var1=x, var2=y) is a duplicate of (var1=y, var2=x) 
# so the Cartesian Product of comboGrid is smaller
covariate_importance <- readRDS(file="C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/boreal_avian_modelling_project/NationalModelsv5/v4/covariate_interactions/02_covariate_importance.rds")



# this loop creates a nested list of `data.frame`s. 
# each top level element a bcr x species x bootstrap tuple, 
# with second-level elements as the rth covariate for that bootstrap.
boot_pts <- list()
for (q in 1:length(gbm_objs)){
  
  # load a bootstrap replicate 
  load(file.path(gbm_objs[q]))
  
  # find evaluation points (`data.frame`) for every covariate (indexed by `r`)
  # this can sped up by reducing `continuous.resolution=100` in `plot.gbm()`
  pts <- list()
  for (r in 1:length(out$var.names)){
    pts[[r]] <- plot.gbm(x=out, return.grid = TRUE, i.var = r, type="response")
  }
  
  boot_pts[[q]] <- pts
  
  # print progress
  cat(paste("\riteration", q))
  Sys.sleep(0.001)
}

saveRDS(boot_pts, file="C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/boreal_avian_modelling_project/NationalModelsv5/v4/covariate_interactions/02_boot_pts.rds")


# from the family of sets (BCRs, species, covariates) 
# group_keys() finds the names of the unique tuples
# There are 373 bcr x common_name x var tuples 
# for CONW BCR60, CAWA BCR60, CONW BCR81, CAWA BCR12
boot_group_keys <- 
  covariate_importance |> 
  dplyr::group_by(bcr, spp, var) |> 
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
  boot_pts_index <- which(sample_id$bcr == key_z$bcr & sample_id$spp == key_z$spp)
  
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
  names(boot_pts_sorted)[z] <- paste(boot_group_keys[z, "bcr"], boot_group_keys[z, "spp"], boot_group_keys[z, "var"], sep = "_")
  
  # print progress
  cat(paste("\riteration", z))
  Sys.sleep(0.001)
}

saveRDS(boot_pts_sorted, file="C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/boreal_avian_modelling_project/NationalModelsv5/v4/covariate_interactions/02_boot_pts_sorted.rds")


# create an index of every bcr x common_name x 2-way interactions 
# `RcppAlgo::comboGrid` is like `expand.grid` and `tidyr::crossing` but avoids duplicates 
# e.g. for our purposes (var1=x, var2=y) is a duplicate of (var1=y, var2=x) 
# so the Cartesian Product of comboGrid is smaller
boot_group_keys_i2 <- 
  RcppAlgos::comboGrid(unique(covariate_importance$var), 
                       unique(covariate_importance$var), 
                       covariate_importance$bcr, 
                       covariate_importance$spp, repetition =  FALSE) |> 
  tibble::as_tibble() |> 
  dplyr::rename(var_1=Var1, var_2=Var2, bcr=Var3, spp=Var4)

# saveRDS(boot_group_keys_i2, file = "C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/boreal_avian_modelling_project/NationalModelsv5/v4/covariate_interactions/boot_group_keys_i2.rds")
