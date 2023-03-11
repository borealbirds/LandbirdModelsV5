# ---
# title: National Models 4.1 - extract covariates
# author: Elly Knight, Melina Houle
# created: December 22, 2022
# ---

#NOTES################################

#This script uses CONUS Albers equal-area conic (EPSG: 5072)

#This script extracts the list of model covariates agreed upon by the BAM National Model team. The full list of those covariates can be found at https://docs.google.com/spreadsheets/d/1XATuq8BOYC2KkbJObaturaf4NFD2ufvn/edit?usp=sharing&ouid=104837701987164094932&rtpof=true&sd=true

#This script runs through each category of covariates separately to avoid mistakes. Future versions should consider aligning file names with a lookup table that details extraction method and then looping through all files at once.

#PREAMBLE############################

#1. Load packages----

library(tidyverse) #basic data wrangling
library(rgee) #to extract data from google earth engine
library(terra) #basic raster handling
library(sf) #basic shapefile handling
library(exactextractr) #fast & efficient raster extraction
library(rgee) #access & extract layers stored on Google Earth Engine

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

#5. Plain dataframe for joining to output----
loc.df <- data.frame(loc.n) %>% 
  dplyr::select(-geometry)

#B. EXTRACT COVARIATES FROM GOOGLE DRIVE####

#1. Climate Normals - static - point value----

#1a. Get list of layers
cn.files <- data.frame(filepath = list.files(file.path(root, "Covariates", "ClimateNormals"), full.names = TRUE)) %>% 
  separate(filepath, into=c("f1", "f2", "f3", "f4", "f5", "f6", "file"), remove=FALSE, sep="/") %>% 
  separate(file, into=c("n1", "beginning", "end", "variable", "tif"), remove=FALSE) %>% 
  mutate(variable = ifelse(tif!="tif", paste0(variable, tif), variable)) %>% 
  dplyr::filter(tif!="18")

#1b. Read in layers
cn <- rast(cn.files$filepath)
names(cn) <- cn.files$variable

#1c. Extract values - use point value
loc.cn <- cbind(loc.df,
                loc.n %>% 
                  st_transform(crs(cn)) %>% 
                  terra::extract(x=cn, ID=FALSE))

#1d. Save output
write.csv(loc.cn, file.path(root, "CovariateExtraction", "ClimateNormals_GD.csv"), row.names = FALSE)

#1e. Cleanup
#rm(cn, cn.files, loc.cn)

#2. Topograhic - static - point or buffer value----

#NOTE: TRI & roughness not working - "cannot read values" - see "Issues" column in tracking sheet for details
#NOTE: Need to fix column name of exact_extract code

#2a. Get list of layers
tp.files <- data.frame(filepath = list.files(file.path(root, "Covariates", "Topography"), full.names = TRUE, pattern="*.tif")) %>% 
  separate(filepath, into=c("f1", "f2", "f3", "f4", "f5", "f6", "file"), remove=FALSE, sep="/") %>% 
  mutate(variable = str_sub(file, -100, -5),
         variable = ifelse(variable=="mTPI_Canada", "mTPI", variable),
         extract = ifelse(variable %in% c("TRI", "roughness"), "radius", "point"))

#2b. Set up loop
#extract iteratively because extents and projections do not match and mix of extraction
loc.tp <- loc.df
for(i in 1:nrow(tp.files)){
  
  #2c. Read in layer
  tp <- rast(tp.files$filepath[i])
  names(tp) <- tp.files$variable[i]
  
  #2d. Extract values - use point value
  if(tp.files$extract[i]=="point"){
    loc.tp <- cbind(loc.tp,
                    loc.n %>% 
                      st_transform(crs(tp)) %>% 
                      terra::extract(x=tp, ID=FALSE))
  }
  
  #2e. Extract values - use buffer mean
  if(tp.files$extract[i]=="radius"){
    loc.tp <- cbind(loc.tp,
                    loc.buff %>% 
                      st_transform(crs(tp)) %>% 
                      exact_extract(x=tp, "mean", colname_fun="values"))
  }
  
  print(paste0("Finished layer ", i, " of ", nrow(tp.files)))
}

#2f. Save output
write.csv(loc.tp, file.path(root, "CovariateExtraction", "Topography_GD.csv"), row.names = FALSE)

#2g. Cleanup
#rm(tp, tp.files, loc.tp)

#3. Wetland - static - buffer value----

#3a. Read in layers
pt <- rast(file.path(root, "Covariates", "Wetlands", "Peat-ML_global_peatland_extent.tif"))
names(pt) <- "peat"

#3b. Extract values - use point value
loc.pt <- cbind(loc.df,
                loc.buff %>% 
                  st_transform(crs(pt)) %>% 
                  exact_extract(x=pt, "mean"))

#3c. Fix column names
colnames(loc.pt) <- c(colnames(loc.df), names(pt))

#3d. Save output
write.csv(loc.pt, file.path(root, "CovariateExtraction", "Wetland_GD.csv"), row.names = FALSE)

#3e. Cleanup
#rm(pt, pt.files, loc.pt)

#4. Disturbance-----

#4a. Get list of layers
dt.files <- data.frame(filepath = list.files(file.path(root, "Covariates", "Disturbance"), full.names = TRUE, pattern="*.tif")) %>% 
  separate(filepath, into=c("f1", "f2", "f3", "f4", "f5", "f6", "file"), remove=FALSE, sep="/") %>% 
  mutate(variable = str_sub(file, -100, -5),
         variable = case_when(variable=="cum_threat2020.02.18" ~ "CAN_HF_1km",
                              variable=="Global_HF_CONUS" ~ "Global_HF_1km",
                              !is.na(variable) ~ variable))

#4b. Set up loop
#extract iteratively because extents and projections do not match
loc.dt <- loc.df
for(i in 1:nrow(dt.files)){
  
  #2c. Read in layer
  dt <- rast(dt.files$filepath[i])
  names(dt) <- dt.files$variable[i]
  
  #2d. Extract values - use point value
  if(dt.files$extract[i]=="point"){
    loc.dt <- cbind(loc.dt,
                    loc.n %>% 
                      st_transform(crs(dt)) %>% 
                      terra::extract(x=dt, ID=FALSE))
  }
  
  #2e. Extract values - use buffer mean
  if(dt.files$extract[i]=="radius"){
    loc.dt <- cbind(loc.dt,
                    loc.buff %>% 
                      st_transform(crs(dt)) %>% 
                      exact_extract(x=dt, "mean", colname_fun="values"))
  }
  
  print(paste0("Finished layer ", i, " of ", nrow(dt.files)))
}

#2f. Save oudtut
write.csv(loc.dt, file.path(root, "CovariateExtraction", "Topography_GD.csv"), row.names = FALSE)

#2g. Cleanup
#rm(dt, dt.files, loc.dt)


#5. Disturbance----

#6. Biomass----


#C. EXTRACT COVARIATES FROM GEE####

#1. Initialize GEE----
ee_Initialize()
ee_check()

#2. Create GCS bucket----
#This only needs to be done once!
# project_id <- ee_get_earthengine_path() %>%
#   list.files(., "\\.json$", full.names = TRUE) %>%
#   jsonlite::read_json() %>%
#   '$'(project_id) # Get the Project ID
# 
# googleCloudStorageR::gcs_create_bucket("national_models", projectId = project_id)

#3. Send polygons to GEE---
poly <- sf_as_ee(loc.buff)
point <- sf_as_ee(loc.n)

#4. Get layers----


#D. GET BCR####

#E. STRATIFY####

#F. 

#G. SAVE#####

save(visit, bird, "")