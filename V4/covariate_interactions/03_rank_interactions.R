# ---
# title: National Models 4.0 - estimating and ranking 2-way covariate interactions for CONWA and CAWA
# author: Mannfred Boehm
# created: September 25, 2024
# ---


#1. attach packages----
print("* attaching packages on master *")
library(gbm)
library(tidyverse)



#2. set root path----
print("* setting root file path *")
root <- "C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/boreal_avian_modelling_project/NationalModelsv5/v4/covariate_interactions/"


boot_group_keys_i2 <- 
  readRDS(file = file.path(root, "rds_files", "boot_group_keys_i2.rds")) |> 
  dplyr::filter(bcr == "BCR_12" & spp == "CAWA")

boot_pts_i2 <- readRDS(file = file.path(root, "results", "boot_pts_i2_cawa12.rds"))


#7. sort and filter results for CONW and CAWA----

# for every zth species x bcr x 2-way interaction tuple (rows in `boot_group_keys_i2`):
# gather the relevant bootstrap predictions from `boot_pts_i2`
# (each element of `boot_pts_i2` has many `y_means`; one for every 2-way interaction)
# and enter each wth `y_mean` (that matches the current 2-way interaction of interest) 
# as a sub-element of [[z]]
# output: a list of lists where the top level elements are species x bcr x 2-way interaction permutations, 
# and the second-level elements are matrices of the associated average and sd of the bootstrapped prediction spaces

# NOTE: `sample_id` must be created in `01_prepare_data.R`
boot_pts_sorted_i2 <- list()
for (z in 1:nrow(boot_group_keys_i2)){
  
  # zth bcr x species x 2-way interaction permutation
  key_z <- boot_group_keys_i2[z,]
  
  # problem: which elements in boot_pts_i2 match the species x bcr tuple defined in key_z?
  # approach: `sample_id` and `boot_pts_i2` have the same length and order (they both are derived from gbm_objs)
  # so we can identify elements in `boot_pts_i2` using info from `sample_id`
  # first, identify a unique species x bcr tuple, then gather covariate interactions into it
  boot_pts_index <- which(sample_id$bcr == key_z$bcr & sample_id$spp == key_z$spp)
  
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
      unlist() |> 
      tibble(value = _) |> 
      summarise(mean = mean(value), sd = sd(value)) 
      
    
    if (!is.na(mean_z$mean) && mean_z$mean >= 0.00010) {
      boot_pts_sorted_i2[[z]] <- mean_z
    } else {
      boot_pts_sorted_i2[[z]] <- list()
    } #close nested else
    
  } # close top level if()
  else {
    boot_pts_sorted_i2[[z]] <- list()
  } # close top level else
  
  names(boot_pts_sorted_i2)[z] <- paste(key_z$bcr, key_z$spp, key_z$var_1, key_z$var_2, sep=".")
  
  # print progress
  cat(paste("\riteration", z, "of", nrow(boot_group_keys_i2)))
  Sys.sleep(0.000001)
  
} # close top level loop


# remove empty elements (output reduced from 8.9MB to 905.2kB;with threshold set and non-threshold elements removed)
nice_var_names <- read_csv("C:/Users/mannf/Proton Drive/mannfredboehm/My files/Drive/boreal_avian_modelling_project/NationalModelsv5/v4/covariate_interactions/results/nice_var_names.csv")


boot_pts_reduced_i2 <- 
  boot_pts_sorted_i2 |> 
  purrr::discard(purrr::is_empty) |>    # Remove empty elements
  purrr::imap_dfr(~ tibble(name = str_remove(.y, "^BCR_12\\.CAWA\\."), mean = .x$mean,  sd = .x$sd)) |> 
  tidyr::extract(name, into = c("covariate_1", "covariate_2"), regex = "(.+)\\.(.+)$") |>
  dplyr::arrange(desc(mean)) |> 
  dplyr::left_join(nice_var_names, by = c("covariate_1" = "var")) |> 
  dplyr::left_join(nice_var_names, by = c("covariate_2" = "var")) |> 
  dplyr::select(var_nice.x, var_nice.y, mean, sd) |> 
  dplyr::slice_max(mean, n=20, with_ties = FALSE) |> 
  

readr::write_csv(boot_pts_reduced_i2, file=file.path(root, "results", "cawa12_interactions.csv"))


