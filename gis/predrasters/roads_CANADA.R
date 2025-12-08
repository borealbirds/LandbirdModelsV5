
################################################3
#  BAM NAM 5.0 - Producing Prediction Rasters for Canada roads
#  The script download, reproject and rescale to 1km and uses a 5x5 moving window
#  Temporal raster are saved on disk and transfered to Google SharedDrive
#################################################
##### road CANADA
library(terra)
library(utils)
library(sf)
library(dplyr)


setwd("E:/MelinaStuff/BAM/NationalModelv5.0")

# Set extent 
rast1k <- rast(nrows=4527, ncols=7300, xmin=-4100000, xmax=3200000, ymin=1673000, ymax=6200000, crs = "EPSG:5072")

## Access urls and extract yearlist
webData <- read.csv("./Scripts/webDataAccess.csv", fileEncoding="UTF-8-BOM")
canroad <- subset(webData, Dataset == "Canada road")
yrlist <- sort(unique(canroad$year))

dwd.folder <- "./CovariateRasters/road/canada"
if (!file.exists(dwd.folder)) {
  dir.create(dwd.folder, showWarnings = FALSE, recursive = TRUE)
}

# Step 1 - Download
for (yr in yearlist) {
  url <- canroad$url[canroad$year==yr]
  download.file(url, file.path(dwd.folder, basename(url)), method = "auto", mode = "wb")
  unzip.folder <- file.path(dwd.folder, yr) 
  if (!file.exists(unzip.folder)) {
    dir.create(unzip.folder, showWarnings = FALSE)
  }
  unzip(file.path(dwd.folder, basename(url)), exdir = unzip.folder)
  file.remove(file.path(dwd.folder, basename(url)))
}

# Step 2 - Process
out.folder <- "./PredictionRasters/road"
if (!file.exists(out.folder)) {
  dir.create(out.folder, showWarnings = FALSE)
}
  
for (yr in yrlist) {
  road_yr <- list.files(file.path(dwd.folder, yr), pattern = "*.shp", full.names = TRUE)
  if(length(road_yr)>0){ # road are shapefile
    roads <- st_read(road_yr)
  }else{  # road are gdb
    road_yr <- list.files(unzipfolder, pattern = "*.gdb", full.names = TRUE)
    layername <- sub('\\..*$', '', list.files(file.path(unzip.folder, yr), pattern = "*.gdb"))
    roads <- sf::st_read(road_yr, layer = layername)
  }
  
  roadcast <- roads %>%
    st_transform(st_crs(rast1k)) %>%
    st_cast("MULTILINESTRING")
  
  roadvec <- vect(roadcast)
  rm(roadcast)
  road.1kden <- rasterizeGeom(roadvec, rast1k, fun="length", unit="km")
  road.5kden<-focal(road.1kden,w=matrix(1,5,5),fun=mean,na.rm=TRUE) #get 5k average at 1k resolution

  out.rast1k <- file.path(out.folder, paste0("canroad_1kden_", yr, ".tif"))
  out.rast5k <- file.path(out.folder, paste0("canroad_5kden_", yr, ".tif"))
  writeRaster(road.1kden, out.rast1k)
  writeRaster(road.5kden, out.rast5k)
}




## When considering the rank
road <- road[,c("isroad", "RANK")]
for(rk in unique(road$RANK)){
    roadcast <- road %>%
      filter(RANK == rk) %>%
      st_transform(st_crs(rast1k)) %>%
      st_cast("MULTILINESTRING") 
    
    roadrast <-rasterize(roadcast, rast1k, fun='min', touches = TRUE, background=NA)
    out.rast100 <- file.path(out.folder, paste0("canroad_1km", as.character(yr), "_rk", rk, ".tif"))
    writeRaster(roadrast, out.rast100)
    rm(roadrast)
}
roadcast <- road %>%
    st_transform(st_crs(rast1k)) %>%
    st_cast("MULTILINESTRING")
  
roadvec <- vect(roadcast)
road.1kden <- rasterizeGeom(roadvec, rast1k, fun="length", unit="km")
rm(roadcast)
out.rast1k <- file.path(out.folder, paste0("canroad_1kmden_", as.character(yr), ".tif"))
writeRaster(road.1kden, out.rast1k)
rm(roadvec)

                                                                                                                                                                                                                                