# ---
# title: National Models 4.1 - extract covariates
# author: Elly Knight, Melina Houle
# created: December 22, 2022
# ---

#NOTES################################

#the coordinates for ABMI Ecosystem health data are buffered within 5.5 km of the actual coordinate and should not be used for covariate extraction.

library(tidyverse)
library(rgee)
library(terra)
library(sf)

#A. EXTRACT COVARIATES FROM GEE####

#B. EXTRACT COVARIATES FROM GD####

#C. GET BCR####

#D. REPACKAGE#####

save(visit, bird, "")