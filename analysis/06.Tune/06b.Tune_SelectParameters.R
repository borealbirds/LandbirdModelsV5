# ---
# title: National Models 5.0 - tune: pick learning rate
# author: Elly Knight
# created: January 2, 2024
# ---

#NOTES################################

#PREAMBLE############################

#1. Load packages----
library(tidyverse) #basic data wrangling

#2. Set root path for data on google drive----
root <- "G:/Shared drives/BAM_NationalModels/NationalModels5.0"

#3. Load the dataset----
load(file.path(root, "Data", "04_NM5.0_data_stratify.R"))

#4. Turn off scientific notation---
options(scipen=99999)

#A. PICK BEST LR###################

#1. Get BCR & bird combo list----
bcr.spp <- birdlist %>% 
  pivot_longer(ABDU:YTVI, names_to="spp", values_to="use") %>% 
  dplyr::filter(use==TRUE)

#2. Get list of results files----
files <- list.files(file.path(root, "Output", "Tuning"))

#3. Read them in----
tune <- purrr::map(.x = files, .f = ~ read.csv())

#4. Pick the learning rate with 1000 - 9999 trees AND highest deviance
ratelist <- tune %>% 
  dplyr::filter(trees < 10000,
                trees > 1000) %>% 
  group_by(spp, bcr) %>% 
  dplyr::filter(deviance.mean = max(deviance.mean))

#5. Save-----
save(visit, cov, bird, offsets, covlist, birdlist, bcrlist, gridlist, ratelist, file=file.path(root, "Data", "06_NM5.0_data_tune.R"))