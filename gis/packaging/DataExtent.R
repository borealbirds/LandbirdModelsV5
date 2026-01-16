###############################################
# title: "Data extent layer"
# author: Elly Knight, Anna Drake
# date: August 8, 2025
##################################################

#This code creates a data availability mask for the final packaged models. We use a threshold of at least 450 surveys within a 250 km radius as the threshold based on occurrence probability analyses.
# false positives at > 450 records are rare (~2% if true occurrence rate=0.1% or 1/1000 (our "absent" threshold)); 
# false negatives at < 450 records are substantial (~19% if true occurrence rate =1% or 1/100 (a rare "present" species))

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
  st_transform(crs=5072) |>
  st_union() |> 
  vect() |> 
  aggregate() |> 
  fillHoles()

#5. Blank raster ----
rast<-rast(file.path(root,"PredictionRasters/AnnualClim/ERAMAP_1km_2019.tif"))
rast[]<-0

#SURVEY DENSITY #############

#1. Make it spatial ----
visit_sp <- st_as_sf(visit, coords=c("lon", "lat"), crs=4326) |> 
  st_transform(crs=5072)

#2. Buffer by 250 km ----
All <- st_buffer(visit_sp, 250000)

#3. Stack and sum buffers - 'rasterize' is memory intensive, so split up the file ----
All2<-list()
rows<-c(seq(1,nrow(All), by=round(nrow(All)/10,0)))
rows[11]<-nrow(All) # include the last few rows in group 10

for (b in c(1:10)){
  All.b<-All[c(rows[b]:rows[b+1]),] # subset
  All2[b]<- rasterize(vect(All.b), rast, fun="sum") # summarize surveys per 250km
}

#4. Put together  -----
Total<-rast(All2)
Total <- app(Total, fun = "sum", na.rm = TRUE)

#5. Filter to areas with at least 450 surveys ----
Total[Total < 450] <- NA
plot(Total)

#POLISHING ##########

#1. Downsampling (for speed) ----
Down <- aggregate(Total, 10, "mean")

#2. Binarize ----
Down[Down > 1] <- 1
plot(Down)

#3. Make it a shapefile and fill holes ----
Shp <- as.polygons(Down, dissolve = TRUE, na.rm=TRUE) |> 
  aggregate() |> 
  fillHoles()

#4. Smooth ----
Smooth <- Shp |> 
  st_as_sf() |> 
  st_buffer(500000) |> 
  st_buffer(-500000)

#5. Crop by bcr extent ----
Out <- Smooth |> 
  vect() |> 
  crop(bcr)
  
#6. Plot ----
plot(bcr)
plot(Out, add=TRUE, col="red")

#7. Save ----
terra::writeVector(Out, file.path(root, "gis", "DataLimitationsMask.shp"))