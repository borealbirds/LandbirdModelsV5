# ---
# title: National Models 4.1 - extract covariates
# author: Elly Knight, Melina Houle
# created: December 22, 2022
# ---

#NOTES################################

#PREAMBLE############################

#1. Load packages----

library(tidyverse) #basic data wrangling
library(rgee) #to extract data from google earth engine
library(terra) #basic raster handling
library(sf) #basic shapefile handling

#2. Set root path for data on google drive----
root <- "G:/.shortcut-targets-by-id/0B1zm_qsix-gPbkpkNGxvaXV0RmM/BAM.SharedDrive/RshProjs/PopnStatus/NationalModelsV4.1/PointCount/"

#A. DATA PREP####

#1. Load data----
load(file.path(root, "01_NM4.1_data_offsets.R"))

#B. EXTRACT COVARIATES FROM GEE####

#C. EXTRACT COVARIATES FROM GD####

#D. GET BCR####

#E. SAVE#####

save(visit, bird, "")