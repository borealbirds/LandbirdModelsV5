
library(terra)
library(googledrive)

# Set extent 
rast1k <- rast(nrows=4527, ncols=7300, xmin=-4100000, xmax=3200000, ymin=1673000, ymax=6200000, crs = "EPSG:5072")

## Access urls and extract yearlist
setwd("E:/MelinaStuff/BAM/NationalModelv4.1")
out.pred <- "./PredictionRasters/wetland"
if (!file.exists(out.pred)) {
  dir.create(out.pred, showWarnings = FALSE)
}

## Access urls and extract yrlist
webData <- read.csv("./Scripts/webDataAccess.csv", fileEncoding="UTF-8-BOM")
peat_lnd <- subset(webData, Dataset == "Peatland")
peat_parentid <- peat_lnd$url

dr <- drive_ls(as_id(peat_parentid), pattern = "Peat")
Peat_rast <- file.path("./CovariateRasters/wetland", dr$name)
drive_download(as_id(dr$id), path = Peat_rast)

## PEATLAND VARIABLES 5x5
Peat<-rast(Peat_rast)
Peat_pj <- terra::project(Peat, rast1k, crop = TRUE)
Peat5x5<-focal(Peat_pj,w=matrix(1,5,5),fun=mean,na.rm=TRUE)#get 5k average at 1k resolution
out.rast5k <- file.path(out.pred, paste0("Peatland5k.tif"))
writeRaster(Peat5x5, out.rast5k, overwrite=TRUE)
