# ---
# title: National Models 4.1 - extract covariates
# author: Elly Knight, Melina Houle
# created: December 22, 2022
# ---

#NOTES################################

#This script uses CONUS Albers equal-area conic (EPSG: 5072)

#This script extracts the list of model covariates agreed upon by the BAM National Model team. The full list of those covariates can be found at https://docs.google.com/spreadsheets/d/1XATuq8BOYC2KkbJObaturaf4NFD2ufvn/edit?usp=sharing&ouid=104837701987164094932&rtpof=true&sd=true

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
loc.cn <- cbind(data.frame(loc.n),
                loc.n %>% 
                  st_transform(crs(cn)) %>% 
                  terra::extract(x=cn, ID=FALSE)) %>% 
  dplyr::select(-geometry)

#1d. Save output
write.csv(loc.cn, file.path(root, "CovariateExtraction", "ClimateNormals.csv"), row.names = FALSE)

#2. Topograhic - static - point value----

#2a. Get list of layers
tp.files <- data.frame(filepath = list.files(file.path(root, "Covariates", "Topography"), full.names = TRUE, pattern="*.tif")) %>% 
  separate(filepath, into=c("f1", "f2", "f3", "f4", "f5", "f6", "file"), remove=FALSE, sep="/") %>% 
  mutate(variable = str_sub(file, -100, -5),
         variable = ifelse(variable=="mTPI_Canada", "mTPI", variable))

#2b. Set up loop
#extract iteratively because extents and projections do not match
for(i in 1:nrow(tp.files)){
  
  #2c. Read in layer
  tp <- rast(tp.files$filepath[i])
  names(tp) <- tp.files$variable[i]
  
  #2d. Extract values - use point value
  if(i==1){
    loc.tp <- cbind(data.frame(loc.n),
                    loc.n %>% 
                      st_transform(crs(tp)) %>% 
                      terra::extract(x=tp, ID=FALSE)) %>% 
      dplyr::select(-geometry)
  } else {
    loc.tp <- cbind(loc.tp,
                    loc.n %>% 
                      st_transform(crs(tp)) %>% 
                      terra::extract(x=tp, ID=FALSE))
    
    loc.tp <- terra::extract(tp, loc.n)
  }
  print(paste0("Finished layer ", i, " of ", nrow(tp.files)))
}

#1d. Save output
write.csv(loc.tp, file.path(root, "CovariateExtraction", "Topography.csv"), row.names = FALSE)

#3. Peatland depth - static - point value----

#4. Road-----

#5. Human footprint----


#C. EXTRACT COVARIATES FROM GEE####


#D. GET BCR####

#E. STRATIFY####

#F. 

#G. SAVE#####

save(visit, bird, "")