###############################################################
# Water proximity raster - using CEC NA Lakes and Rivers 2023
# https://www.cec.org/north-american-environmental-atlas/lakes-and-rivers-2023/
# March 20, 2025
# Anna Drake
###############################################################

# This code uses a shapefile of North American Lakes and Rivers to produce a binary
# "proximity to water" layer at 100m resolution. A pixel is proximate to water (1) if it falls within a 200m buffer of a lake or river
# This layer is coarsened to 1km predictive layer using nearest-neighbour sampling

require(terra)
require(data.table)
require(dplyr)
require(sf)
require(tidyverse)

# Directory ----
root<-"G:/Shared drives/BAM_NationalModels5"
out<-"C:/Users/andrake/Desktop/WaterExtraction"
setwd(out)

# BCRs ----
bcr<- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel.shp")) |>
  st_transform(crs=crs("EPSG:5072"))

# Template -----
temp<-rast(file.path(root,"PredictionRasters/Wetland/1km/WetOccur_1km.tif"))

# Buffer lakes by 200m ---
Lakes<-read_sf("C:/Users/andrake/Desktop/rivers_and_lakes_shapefile/rivers_and_lakes_shapefile/NA_Lakes_and_Rivers/data/lakes_p/northamerica_lakes_cec_2023.shp") |> 
  st_transform(crs=crs("EPSG:5072"))
Lakes<-Lakes|>st_crop(bcr)|> st_buffer(200)|>st_union()  # 200 m buffer

# Rasterize lakes ---
template = rast(vect(Lakes),res=100)
LakesRast <- rasterize(vect(Lakes), template)
LakesRast<- subst(LakesRast, NA, 0)

# Buffer rivers by 200 m ---
River<-read_sf("C:/Users/andrake/Desktop/rivers_and_lakes_shapefile/rivers_and_lakes_shapefile/NA_Lakes_and_Rivers/data/rivers_l/northamerica_rivers_cec_2023.shp")|>
  st_transform(crs=crs("EPSG:5072"))
River<-River|>st_crop(bcr)|>st_buffer(200)|>st_union()  # 200 m buffer

# Rasterize rivers ----
template = rast(vect(River),res=100)
RiverRast <- rasterize(vect(River), template)
RiverRast<- subst(RiverRast, NA, 0)

# Match layer extents -----
RiverRast_e<-RiverRast|>extend(LakesRast)|>crop(LakesRast)
LakesRast_e<-LakesRast|>project(RiverRast_e, method="near")|>
             extend(RiverRast_e)|>crop(RiverRast_e)

# Merge ----
WaterRast<-RiverRast_e+LakesRast_e
WaterRast<- subst(WaterRast, 2, 1)

# Write out 100m res ---
writeRaster(WaterRast,"Buffered200_Water_100mres.tif", overwrite=T)

# Coarsen to 1km using template ----
WaterRast_1km<-WaterRast|>project(temp,method="near")|>extend(temp)|>crop(temp)

# Write out 1km res ---
writeRaster(WaterRast_1km,"Buffered200_Water_1kmres.tif", overwrite=T)

################ End of Processing ###################
