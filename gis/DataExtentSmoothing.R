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
library(terra) #raster management

#2. Set root path for data on google drive----
root <- "G:/Shared drives/BAM_NationalModels5"

#3. Get tif----
na1 <- rast(file.path(root, "MosaicWeighting", "Range", "ABS_NA_Layer450.tiff"))

#4. Get Canada extent----
can <- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel_Unbuffered_Canada.shp")) |> 
  st_union() |> 
  st_transform(crs(na1)) |> 
  vect()

#5. Crop ----
na2 <- na1 |> crop(can, mask=TRUE)

#6. Smooth it ----
w <- matrix(1, 299, 299)
na3 <- focal(na2, w=w, fun=mean, na.policy="omit", expand=TRUE)

#7. Rebinarize ----
na4 <- na3 > 0.5

#8. Fix boundary ----
na5 <- resample(na4, na2)
na5[is.na(na5)] <- 0
na6 <- sum(na5, na2)
na6[values(na6) > 0] <- 1

#8. Plot ----
plot(na6)

#9. Save ----
writeRaster(na6, file.path(root, "gis", "DataLimitationsMask.tif"))

