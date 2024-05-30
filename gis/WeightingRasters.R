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

#PREAMBLE##########

#1. Load packages ----
library(sf)
library(tidyverse)
library(terra)

#2. NAD83(NSRS2007)/Conus Albers projection (epsg:5072) ----
crs <- "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"

#3. Set file paths ----
root <- "G:/Shared drives/BAM_NationalModels/NationalModels5.0" #PC link
setwd(file.path(root,"MosaicWeighting"))

#4. Read in country shapefiles----
can <- read_sf(file.path(root, "Regions", "CAN_adm", "CAN_adm0.shp")) |> 
  st_transform(crs=5072) 
usa <- read_sf(file.path(root, "Regions", "USA_adm", "USA_adm0.shp")) |> 
  st_transform(crs=5072) 

#5. Read in BCR shapefile----
#Remove subunit 1
bcr <- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel.shp")) |> 
  dplyr::filter(subUnit!=1) |> 
  st_transform(crs=5072)

#6. Identify BCRs for each country----
bcr.ca <- bcr |> 
  st_intersection(can) |> 
  mutate(country="can") |> 
  dplyr::select(subUnit, country)

bcr.usa <- bcr |> 
  st_intersection(usa) |> 
  mutate(country="usa") |> 
  dplyr::select(subUnit, country)

#OVERLAP#######

#1. Determine overlap zones for Canada and US  -----

Baselayer<-terra::rast(file.path(root, "PredictionRasters/ClimateNormal/AHM_1km.tif")) #any prediction raster
Baselayer[!is.na(Baselayer)] <- 1
e<-ext(Baselayer)

#2. Canada----
OverlapCan <- c()

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

#2, Mosaic----
c.mosaic <- sprc(OverlapCan) %>% mosaic(., fun="sum")  # == amount of prediction layer overlap
writeRaster(c.mosaic,paste(root,"MosaicWeighting/","ModelOverlap_Can.tif"))

#3. Reduce to shapefile of overlap areas----
c.mosaic[c.mosaic==1]<-0
c.mosaic[c.mosaic>1]<-1 
OverlapCan <- as.polygons(c.mosaic) %>% .[2]
writeVector(OverlapCan,paste(root,"MosaicWeighting/","BCR_Overlap_Can",sep=""))

#4. Now USA----
OverlapUS<-c()

for (i in (1:nrow(bcr.usa))) {
  bcr.i <-  bcr.usa %>% 
    dplyr::filter(row_number()==i) %>% 
    st_buffer(100000)%>%  # Filter & buffer shapefile
    st_intersection(., usa) #crop at USA border
  
  OverlapUS[[i]]<-Baselayer %>% 
    project(.,crs)  %>%
    crop(.,e) %>%
    mask(.,bcr.i) 
}

#5. Mosaic----
u.mosaic <- sprc(OverlapUS)%>%mosaic(., fun="sum")
writeRaster(u.mosaic,paste(root,"MosaicWeighting/","ModelOverlap_US.tif", sep=""))

#6. Reduce to a shapefile of overlap areas----
u.mosaic[u.mosaic==1]<-0
u.mosaic[u.mosaic>1]<-1 

OverlapUS <- as.polygons(u.mosaic) %>% .[2]
writeVector(OverlapUS,paste(root,"MosaicWeighting/","BCR_Overlap_US", sep=""))

#7. Combine the two----
MosaicOverlap1<-terra::rast(file.path(root,"MosaicWeighting","ModelOverlap_Can.tif"))
MosaicOverlap2<-terra::rast(file.path(root,"MosaicWeighting","ModelOverlap_US.tif"))
MosaicOverlap <- mosaic(MosaicOverlap1, MosaicOverlap2)
writeRaster(MosaicOverlap, file.path(root, "gis","ModelOverlap.tif", sep=""), overwrite=TRUE)

#DISTANCE TO EDGE######

#For model weighting when mosaicing

#1. Make merged subunits----

#See "gis/BufferedSubunits.R" for rationale

bcr.can4142 <- bcr.ca |> 
  dplyr::filter(subUnit %in% c(41, 42)) |> 
  st_union() |> 
  nngeo::st_remove_holes() |> 
  st_sf() |> 
  mutate(country="can",
         subUnit = 4142)
st_geometry(bcr.can4142) <- "geometry"

bcr.usa41423 <- bcr.ca |> 
  dplyr::filter(subUnit %in% c(41, 42, 3)) |> 
  st_union() |> 
  nngeo::st_remove_holes() |> 
  st_sf() |> 
  mutate(country="usa",
         subUnit = 41423)
st_geometry(bcr.usa41423) <- "geometry"

bcr.usa414232 <- bcr.usa |> 
  dplyr::filter(subUnit %in% c(41, 42, 3, 2)) |> 
  st_union() |> 
  nngeo::st_remove_holes() |> 
  st_sf() |> 
  mutate(country="usa",
         subUnit = 414232)
st_geometry(bcr.usa414232) <- "geometry"

bcr.can8182 <- bcr.ca |> 
  dplyr::filter(subUnit %in% c(81, 82)) |> 
  st_union() |> 
  nngeo::st_remove_holes() |> 
  st_sf() |> 
  mutate(country="can",
         subUnit = 8182)
st_geometry(bcr.can8182) <- "geometry"

#polygons to remove
bcr.remove <- data.frame(country=c("can", "usa", "usa", "can", "usa"),
                         subUnit=c(41, 42, 3, 42, 41))

#2. Put together----
bcr.country <- rbind(bcr.ca, bcr.usa, bcr.can4142, bcr.usa41423, bcr.usa414232, bcr.can8182) |> 
  anti_join(bcr.remove) |> 
  mutate(bcr = paste0(country, subUnit))

#3. Set up loop for BCR buffering----
for(i in 1:nrow(bcr.country)){
  
  #4. Filter & buffer shapefile----
  bcr.i <- bcr.country[i,] |> 
    st_buffer(100000)
  
  #5. Make fishnet----
  #(10km2 to speed processing - could make this finer)
  grid <- st_make_grid(bcr.i, cellsize = 10000, what = "centers") |> 
    st_intersection(bcr.i)
  
  #6. Polygon -> line----
  bcr.line <- st_cast(bcr.i, "MULTILINESTRING")
  
  #7. Get distance from the edge for each point---
  dist <- st_distance(bcr.line, grid)
  
  #8. Dataframe as proportion of max distance----
  df <- data.frame(dist = as.vector(dist)/as.vector(max(dist)),
                   st_coordinates(grid))
  
  #9. Rasterize + match resolution----
  r <- terra::rast(resolution = 10000, ext = ext(bcr.line), crs = crs)
  
  dist_sf <- st_as_sf(df, coords = c("X", "Y")) |> 
    st_set_crs(crs)
  
  dist_raster <- rasterize(dist_sf, r, "dist", fun = mean) |> 
    terra::disagg(fact=10, method="bilinear")
  
  #10. Make bcr raster----
  bcr.r <- bcr.country[i,] |> 
    rasterize(r) |> 
    terra::disagg(fact=10, method="bilinear")
  
  #11. Invert----
  bcr.inv <- bcr.r
  bcr.inv[bcr.inv==1] <- 0 
  bcr.inv[is.na(bcr.inv)] <- 1
  
  #11. Zero out the center----
  dist_zero <- dist_raster*bcr.inv
  
  #12. Rescale----
  dist_scale <- dist_zero/dist_zero@cpp[["range_max"]]
  
  #13. Fill in the center----
  dist_fill <- list(dist_scale, bcr.r) |> 
    sprc() |> 
    mosaic(fun="sum") |> 
    crop(st_transform(bcr.i, crs(MosaicOverlap)))
  
  #14. Reduce other areas of overlap----
  overlap.i <- crop(MosaicOverlap, st_transform(bcr.i, crs(MosaicOverlap)), mask=TRUE) |>
    round() |> 
    project(crs(dist_fill))
  overlap.i[overlap.i < 3] <- 1
  overlap.i[overlap.i==3] <- 2/3
  overlap.i[overlap.i==4] <- 1/2
  
  dist_out <- overlap.i*resample(dist_fill, overlap.i)
  dist_out[dist_out > 1] <- 1
  plot(dist_out)
  
  #15. Crop to international boundary----
  if(bcr.i$country=="can"){ bcr.out <- mask(dist_out, can)}
  if(bcr.i$country=="usa"){ bcr.out <- mask(dist_out, usa)}
  
  #16. Save----
  writeRaster(bcr.out, paste0(root, "/gis/edgeweights/", bcr.i$bcr, ".tif"), overwrite=TRUE)
  
  print(paste0("Finished bcr ", i, " of ", nrow(bcr.country)))
  
}

#17. Check it----
check <- list.files(file.path(root, "gis", "edgeweights"), full.names = TRUE)[-c(17,31)] |> 
  lapply(rast) |> 
  sprc() |> 
  mosaic(fun="sum")
plot(check)

#Nope nope nopeity nope


###################### END OF CODE #####################################
