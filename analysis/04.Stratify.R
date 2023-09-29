# ---
# title: National Models 4.1 - stratify by BCR*country
# author: Elly Knight
# created: September 29, 2022
# ---

#NOTES################################

#TODO: FIND GOOD LAYER TO INTERSECT BCRS WITH
#TODO: THINK ABOUT BORDER TRANSITION & BUFFERING VS LAYER CONSISTENCY




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

#2. Read in us/canada shapefile----
uscan <- read_sf(file.path(root, "Regions", "USA_Canada", "USA_Canada_ShapefileMErge.shp")) %>% 
  mutate(country = ifelse(StateName %in% c("NUNAVUT", "NORTHWEST TERRITORIES", "YUKON TERRITORY", "BRITISH COLUMBIA", "QUEBEC", "NEWFOUNDLAND AND LABRADOR", "ALBERTA", "SASKATCHEWAN", "MANITOBA", "ONTARIO", "QUEBEC", "NOVA SCOTIA", "NEW BRUNSWICK", "PRINCE EDWARD ISLAND"), "CA", "USA"))

can <- uscan %>% 
  dplyr::filter(country=="CA") %>% 
  st_union()

#2. Read in BCR shapefile----
bcrs <- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel.shp")) %>% 
  mutate(bcr=paste0("bcr", subUnit))

#3. Set up loop----
bcr.list <- list(bcr=visit$id)
for(i in 1:nrow(bcrs)){
  
  #4. Filter & buffer bcr shapefile----
  bcr.shp <- bcrs %>% 
    dplyr::filter(row_number()==i) %>% 
    st_transform(crs=5072) %>% 
    st_buffer(100000) %>% 
    mutate(bcr=1)
  
  if(class(bcr.shp)[1]=="sf"){
    
    #5. Convert to raster for fast extraction----
    r <- rast(ext(bcr.shp), resolution=1000, crs=crs(bcr.shp))
    bcr.r <- rasterize(x=bcr.shp, y=r, field="bcr")
    
    #6. Extract raster value----
    bcr.list[[i+1]] <- extract(x=bcr.r, y=visit.v)[,2]
  }
  

  
  print(paste0("Finished bcr ", i, " of ", nrow(bcrs)))
  
}

#7. Collapse to dataframe and name----
bcr <- data.frame(do.call(cbind, bcr.list))
colnames(bcr) <- c("id", bcrs$bcr)

#8. Remove points outside of study area----
visit.bcr <- 




#THIN FOR EACH BOOTSTRAP##############




#UPDATE BIRD LIST BY BCR##############