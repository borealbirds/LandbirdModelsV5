# ---
# title: National Models 4.1 - stratify by BCR*country
# author: Elly Knight
# created: September 29, 2022
# ---

#NOTES################################

#In this script, we stratify & prepare the data for modelling. Steps include:
#1. BCR Attribution: Diving the data into BCRs for separate models. Each BCR is buffered by 100km. We do this so that we can feather predictions from adjacent regions together. The exception is the international boundaries between US and Canada, which we don't buffer because spatial layers for the covariates are different on either side of the border. In this case, we use a shapefiles of country boundaries to intersect with the buffered regions.




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

#2. Read in country shapefiles----
can <- read_sf(file.path(root, "Regions", "CAN_adm", "CAN_adm0.shp")) %>% 
  st_transform(crs=5072) 
usa <- read_sf(file.path(root, "Regions", "USA_adm", "USA_adm0.shp")) %>% 
  st_transform(crs=5072) 

#3. Read in BCR shapefile----
bcr <- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel.shp")) %>% 
  mutate(bcr=paste0("bcr", subUnit)) %>% 
  st_transform(crs=5072)

ggplot() +
  geom_sf(data=bcr, aes(fill=bcr))
  
#3. Set up loop for Canada BCR----
bcr.ca.list <- list(id=visit$id)
for(i in 1:nrow(bcr)){
  
  #4. Filter, buffer, & crop bcr shapefile to Canada----
  bcr.ca <- bcr %>% 
    dplyr::filter(row_number()==i) %>% 
    st_buffer(100000) %>% 
    mutate(bcr=1) %>% 
    st_intersection(can)
  
  if(nrow(bcr.ca) > 0){
    #5. Convert to raster for fast extraction----
    r <- rast(ext(bcr.ca), resolution=1000, crs=crs(bcr.ca))
    bcr.r <- rasterize(x=bcr.ca, y=r, field="bcr")
    
    #6. Extract raster value----
    bcr.ca.list[[i+1]] <- extract(x=bcr.r, y=visit.v)[,2]
    names(bcr.ca.list)[i+1] <- paste0("ca-", bcr$bcr[i])
  }
  
  print(paste0("Finished bcr ", i, " of ", nrow(bcr)))
  
}

#7. Set up loop for USA BCR----
bcr.usa.list <- list(id=visit$id)
for(i in 1:nrow(bcr)){
  
  #8. Filter, buffer, & crop bcr shapefile to Canada----
  bcr.usa <- bcr %>% 
    dplyr::filter(row_number()==i) %>% 
    st_buffer(100000) %>% 
    mutate(bcr=1) %>% 
    st_intersection(usa)
  
  if(nrow(bcr.usa) > 0){
    #9. Convert to raster for fast extraction----
    r <- rast(ext(bcr.usa), resolution=1000, crs=crs(bcr.usa))
    bcr.r <- rasterize(x=bcr.usa, y=r, field="bcr")
    
    #10. Extract raster value----
    bcr.usa.list[[i+1]] <- extract(x=bcr.r, y=visit.v)[,2]
    names(bcr.usa.list)[i+1] <- paste0("usa-", bcr$bcr[i])
  }

  print(paste0("Finished bcr ", i, " of ", nrow(bcr)))
  
}

#11. Collapse to dataframe and name----
bcr.df <- data.frame(do.call(cbind, bcr.ca.list)) %>% 
  left_join(data.frame(do.call(cbind, bcr.usa.list)))

#12. Remove columns (combos of BCR & country) that don't exist---
bcr.df <- bcr.df[colSums(!is.na(bcr.df))>0]

#13. Remove points outside of study area----
bcr.df <- bcr.df[rowSums(!is.na(bcr.df[,c(2:ncol(bcr.df))]))>0,]
visit.bcr <- visit %>% 
  dplyr::filter(id %in% bcr.df$id)




#THIN FOR EACH BOOTSTRAP##############




#UPDATE BIRD LIST BY BCR##############

#SAVE#####

save(visit.bcr, file=file.path(root, "03_NM4.1_data_stratify.Rdata"))
