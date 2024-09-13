# ---
# title: National Models 4.0 - estimating and ranking 2-way covariate interactions for CONWA and CAWA
# author: Mannfred Boehm
# created: September 13, 2024
# ---

#NOTES################################

# This script extracts covariate contributions to model predictions, and also synthesizes various trait databases.

# The outputs will be analysed for identifying covariates of importance across BCRs



#1. Attach packages----
library(gbm)
library(RcppAlgos)
library(tidyverse)



#2. Create a list of dataframes containing relative influence per covariate----
#   Every list element represents a bootstrap sample 

# connect to BAM Drive and find bootstrap files 
root <- "G:/Shared drives/BAM_NationalModels4/NationalModels4.0"

# gbm objects for CAWA and CONWA stored on the BAM drive
gbm_objs <- list.files(file.path(root, "output", "bootstraps"))[sample(1:60, 3)]

# 


