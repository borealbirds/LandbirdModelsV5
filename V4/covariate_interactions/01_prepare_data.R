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
dirs <- c("CAWA/BCR_60")

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
#saveRDS(covariate_importance,  file="C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/boreal_avian_modelling_project/NationalModelsv5/v4/covariate_interactions/covariate_importance.rds")




# create an index of every bcr x common_name x 2-way interactions 
# `RcppAlgo::comboGrid` is like `expand.grid` and `tidyr::crossing` but avoids duplicates 
# e.g. for our purposes (var1=x, var2=y) is a duplicate of (var1=y, var2=x) 
# so the Cartesian Product of comboGrid is smaller
covariate_importance <- readRDS(file="C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/boreal_avian_modelling_project/NationalModelsv5/v4/covariate_interactions/covariate_importance.rds")

boot_group_keys_i2 <- 
  RcppAlgos::comboGrid(unique(covariate_importance$var), 
                       unique(covariate_importance$var), 
                       covariate_importance$bcr, 
                       covariate_importance$spp, repetition =  FALSE) |> 
  tibble::as_tibble() |> 
  dplyr::rename(var_1=Var1, var_2=Var2, bcr=Var3, spp=Var4)

# saveRDS(boot_group_keys_i2, file = "C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/boreal_avian_modelling_project/NationalModelsv5/v4/covariate_interactions/boot_group_keys_i2.rds")
