# ---
# title: National Models 4.1 - stratify by BCR*country
# author: Elly Knight
# created: September 29, 2022
# ---

#NOTES################################



#PREAMBLE############################

#1. Load packages----

library(tidyverse) #basic data wrangling
library(terra) #basic raster handling
library(sf) #basic shapefile handling
library(fasterize) #fast rasterization of shapefiles
library(exactextractr) #fast & efficient raster extraction

#2. Set root path for data on google drive----
root <- "G:/Shared drives/BAM_NationalModels/NationalModels4.1"

#3. Load data packages with offsets and covariates----
load(file.path(root, "03_NM4.1_data_covariates.R"))

#BCR ATTRIBUTION#####################

#1. Make visit a vector object----
visit.v <- visit %>% 
  dplyr::select(id, lat, lon) %>% 
  st_as_sf(coords=c("lon", "lat"), crs=4326) %>% 
  st_transform(5072) %>% 
  vect()

#2. Get list of buffered BCRs----
bcrs <- data.frame(path=list.files(file.path(root, "Regions", "Buffered_BCRs"), pattern="*.shp", full.names = TRUE),
                   filename=list.files(file.path(root, "Regions", "Buffered_BCRs"), pattern="*.shp")) %>% 
  mutate(bcr=str_sub(filename, -100, -5))

#3. Set up loop----
bcr.list <- list(bcr=visit$id)
for(i in 1:nrow(bcrs)){
  
  #4. Read in buffered bcr shapefile----
  bcr.shp <- try(read_sf(bcrs$path[i]) %>% 
    st_transform(crs=5072) %>% 
    mutate(bcr=1))
  
  if(class(bcr.shp)[1]=="sf"){
    
    #5. Convert to raster for fast extraction----
    r <- rast(ext(bcr.shp), resolution=100, crs=crs(bcr.shp))
    bcr.r <- rasterize(x=bcr.shp, y=r, field="bcr")
    
    #6. Extract raster value----
    bcr.list[[i+1]] <- extract(x=bcr.r, y=visit.v)[,2]
  }
  

  
  print(paste0("Finished bcr ", i, " of ", nrow(bcrs)))
  
}

#7. Collapse to dataframe and name----
bcr <- data.frame(do.call(cbind, bcr.list))
colnames(bcr) <- c("id", bcrs$bcr)

#UPDATE BIRD LIST BY BCR##############




#THIN FOR EACH BOOTSTRAP##############