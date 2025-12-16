###############################################
# title: "Smoothing data range extent layer"
# author: Elly Knight
# date: August 8, 2025
##################################################

#This code takes the NA tif output from "Probable_Range_Extent.R" and smooths it for masking packaged outputs to the area of available survey data.

#PREAMBLE############################

#1. Load packages----

library(tidyverse) #basic data wrangling
library(sf) #read in & handle study area
library(sp) #for adehabitatHR compatibility
library(terra) #raster management
library(adehabitatHR) #area calculation
library(ks) #KDE

#2. Set root path for data on google drive----
root <- "G:/Shared drives/BAM_NationalModels5"

#3. Load dataset ----
load(file.path(root, "data", "04_NM5.0_data_stratify.Rdata"))

#4. Load BCRs ----
bcr <- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel.shp")) |> 
  st_transform(crs=4326)

#CALCULATE EXTENT#################

#1. Get bootstrap 1 ----
visit_1 <- visit |> 
  dplyr::filter(id %in% bootlist$b1)

#2. Make it spatial ----
pts <- visit_1 |> 
  st_as_sf(coords=c("lon", "lat"), crs=4326) 

#3. Try convex hull ----
mcp <- pts |> 
  st_union() |> 
  st_concave_hull(ratio = 0.2) 

#4. Visualize ----
ggplot() +
  geom_sf(data=bcr, aes(fill=factor(subUnit))) +
  geom_sf(data=st_as_sf(mcp), fill=NA, linewidth = 1) +
  geom_point(data=visit_1, aes(x=lon, y=lat))

#5. Make it a raster ----
r <- bcr |> 
  st_transform(crs=5072) |> 
  rast(resolution = 1000)

mcp_v <- mcp |> 
  st_as_sf() |> 
  st_transform(crs=5072) |> 
  vect()

mcp_r <- rasterize(mcp_v, r, field=1, background=NA)
plot(mcp_r)

#compare
na1 <- rast(file.path(root, "MosaicWeighting", "Range", "ABS_NA_Layer450.tiff"))
plot(na1)

#6. Smooth it ----
w <- matrix(1, 299, 299)
mcp_smooth <- focal(mcp_r, w=w, fun=mean, na.policy="omit", expand=TRUE)
plot(mcp_smooth)

#7. Save ----
writeRaster(mcp_smooth, file.path(root, "gis", "DataExtentMask.tif"))


