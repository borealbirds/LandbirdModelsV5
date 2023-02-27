# ---
# title: National Models 4.1 - extract covariates
# author: Elly Knight, Melina Houle
# created: December 22, 2022
# ---

#NOTES################################

#This script uses CONUS Albers equal-area conic (EPSG: 5072)

#PREAMBLE############################

#1. Load packages----

library(tidyverse) #basic data wrangling
library(rgee) #to extract data from google earth engine
library(terra) #basic raster handling
library(sf) #basic shapefile handling
library(exactextractr) #fast & efficient raster extraction

#2. Set root path for data on google drive----
root <- "G:/Shared drives/BAM_NationalModels/NationalModels4.1"

#A. DATA PREP####

#1. Load data----
load(file.path(root, "BirdData", "02_NM4.1_data_offsets.R"))

#2. Create sf object of just location*year for annual layers----
#transform to 5072
loc.yr <- visit %>% 
  dplyr::select(project, location, lat, lon, year) %>% 
  unique() %>% 
  st_as_sf(coords=c("lon", "lat"), crs=4326, remove=FALSE) %>% 
  st_transform(crs=5072)

#3. Create sf object of just location for static layers----
loc <- loc.yr %>% 
  dplyr::select(-year) %>% 
  unique()

#EXTRA. Test sample dataset----
set.seed(1234)
loc.n <- loc %>% 
  sample_n(1000)

#4. Buffer location objects----
loc.buff <- st_buffer(loc.n, 200)

#B. EXTRACT COVARIATES FROM GOOGLE DRIVE####

#1. Climate Normals----

#1a. Get list of layers
cn.files <- data.frame(filepath = list.files(file.path(root, "Covariates", "ClimateNormals"), full.names = TRUE)) %>% 
  separate(filepath, into=c("f1", "f2", "f3", "f4", "f5", "f6", "file"), remove=FALSE, sep="/") %>% 
  separate(file, into=c("n1", "beginning", "end", "variable", "tif"), remove=FALSE) %>% 
  mutate(variable = ifelse(tif!="tif", paste0(variable, tif), variable)) %>% 
  dplyr::filter(tif!="18")

#1b. Read in layers
cn <- rast(cn.files$filepath)
names(cn) <- cn.files$variable

#1c. Extract values
loc.cn <- exact_extract(cn, loc.buff, 'mean')




#C. EXTRACT COVARIATES FROM GEE####


#D. GET BCR####

#E. STRATIFY####

#F> 

#F. SAVE#####

save(visit, bird, "")