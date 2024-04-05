#-----------------------------------
#  Produce weighting rasters for
#  Prediction Mosaicking function
#  Author: Anna Drake
#  Created: April 2024
#-----------------------------------

# This code produces 3 outputs: 
# 1,2. A scored raster of BCR overlap for US and Canada + a shapefile of overlap zones 
# 3. Distance-weighting rasters for every BCR (decreasing weight toward edge)
# These outputs are generic and are used within species-prediction mosaic code 
# Distance weighting uses approach by Dominic Roye (dominicroye.github.io)
# All outputs go into "MosaicWeighting" folder

#Load packages ----

library(sf)
library(tidyverse)
library(terra)

#1. NAD83(NSRS2007)/Conus Albers projection (epsg:5072) ----
crs<-"+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"

#2. Set file paths ----
root<-"G:/Shared drives/BAM_NationalModels/NationalModels5.0" #PC link
setwd(file.path(root,"MosaicWeighting"))

#3.Get Canada and USA BCR boundaries ----

usa <- read_sf(file.path(root, "Regions", "USA_adm", "USA_adm0.shp")) %>% 
  st_transform(crs=crs) # US boundary

can <- read_sf(file.path(root, "Regions", "CAN_adm", "CAN_adm0.shp")) %>% 
  st_transform(crs=crs) # Canadian boundary

bcr.us <- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel.shp")) %>% 
  dplyr::filter(subUnit!=1) %>% # remove BCR 1
  st_transform(crs=crs) %>%
  st_intersection(usa) %>% 
  mutate(country="us") # Restrict to USA BCRs

bcr.ca <- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel.shp")) %>% 
  dplyr::filter(subUnit!=1) %>% # remove BCR 1
  st_transform(crs=crs) %>%
  st_intersection(can) %>% 
  mutate(country="ca") # Restrict to Canadian BCRs

#4. Determine overlap zones for Canada and US  -----

Baselayer<-terra::rast(file.path(root, "PredictionRasters/ClimateNormal/AHM_1km.tif")) #any prediction raster
Baselayer[!is.na(Baselayer)] <- 1
e<-ext(Baselayer)

OverlapCan<-c()

for (i in (1:nrow(bcr.ca))) {
  bcr.i <-  bcr.ca %>% 
    dplyr::filter(row_number()==i) %>% 
    st_buffer(100000)%>%  # Filter & buffer shapefile
    st_intersection(., can) #crop at USA border
  
  OverlapCan[[i]]<-Baselayer %>% 
    project(.,crs)  %>%
    crop(.,e) %>%
    mask(.,bcr.i) 
}

c.mosaic <- sprc(OverlapCan) %>% mosaic(., fun="sum")  # == amount of prediction layer overlap
writeRaster(c.mosaic,paste(root,"MosaicWeighting/","ModelOverlap_Can.tif"))

#reduce to shapefile of overlap areas
c.mosaic[c.mosaic==1]<-0
c.mosaic[c.mosaic>1]<-1 
OverlapCan <- as.polygons(c.mosaic) %>% .[2]
writeVector(OverlapCan,paste(root,"MosaicWeighting/","BCR_Overlap_Can",sep=""))

#Now USA
OverlapUS<-c()

for (i in (1:nrow(bcr.us))) {
  bcr.i <-  bcr.us %>% 
    dplyr::filter(row_number()==i) %>% 
    st_buffer(100000)%>%  # Filter & buffer shapefile
    st_intersection(., usa) #crop at USA border
  
  OverlapUS[[i]]<-Baselayer %>% 
    project(.,crs)  %>%
    crop(.,e) %>%
    mask(.,bcr.i) 
}

u.mosaic <- sprc(OverlapUS)%>%mosaic(., fun="sum")
writeRaster(u.mosaic,paste(root,"MosaicWeighting/","ModelOverlap_US.tif", sep=""))

#reduce to a shapefile of overlap areas
u.mosaic[u.mosaic==1]<-0
u.mosaic[u.mosaic>1]<-1 

OverlapUS <- as.polygons(u.mosaic) %>% .[2]
writeVector(OverlapUS,paste(root,"MosaicWeighting/","BCR_Overlap_US", sep=""))

#5.Produce distance-to-the-edge raster for every BCR (for model weighting) ---------------
# ***Note*** this extracts Canada, change selection where noted to produce US BCR outputs

for (i in (1:nrow(bcr.ca))) {  #***for US BCRs: bcr.ca==bcr.us/can==usa
  bcr.i <-  bcr.ca %>% 
    dplyr::filter(row_number()==i) %>% 
    st_buffer(100000)%>%  # Filter & buffer shapefile
    st_intersection(., can) #crop at CAN border
  
# Make fishnet (10km2 to speed processing - could make this finer)
  grid <- st_make_grid(bcr.i, cellsize = 10000, what = "centers")
  grid <- st_intersection(grid, bcr.i)  
  
# Polygon -> line
  bcr.i <- st_cast(bcr.i, "MULTILINESTRING")
  
# Get distance from the edge for each point
  dist <- st_distance(bcr.i, grid)
  
# Dataframe as proportion of max distance
  df <- data.frame(dist = as.vector(dist)/as.vector(max(dist)),
                   st_coordinates(grid))
  
# Rasterize + match resolution
  ext <- ext(bcr.i)
  r <- terra::rast(resolution = 10000, ext = ext, crs = crs)
  dist_sf <- st_as_sf(df, coords = c("X", "Y")) %>% st_set_crs(crs)
  dist_raster <- terra::rasterize(dist_sf, r, "dist", fun = mean) %>% terra::disagg(.,fact=10, method="bilinear")
  writeRaster(dist_raster,paste(root,"MosaicWeighting/CAN_BCR_Weighting/","EdgeWeighting_BCR_", bcr.i$subUnit, ".tif", sep=""))
}  

#*** For US BCRs change save to: writeRaster(dist_raster,paste(root,"MosaicWeighting/USA_BCR_Weighting/","EdgeWeighting_USA_BCR_", bcr.i$subUnit, ".tif", sep=""))

###################### END OF CODE #####################################
