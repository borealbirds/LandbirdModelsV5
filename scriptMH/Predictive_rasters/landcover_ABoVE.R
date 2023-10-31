################################################3
#  BAM NAM 4.1 - Producing Prediction Rasters for Landcover ABoVE
#  The script download from GSD, reproject and rescale to 1km and uses a 5x5 moving window
#  Temporal raster are saved on disk and transfered to Google SharedDrive
#################################################
library(terra)
library(googledrive)

########################
### PARAMS
########################
# workdir
setwd("E:/MelinaStuff/BAM/NationalModelv4.1")

# Set extent 
rast1k <- rast(nrows=4527, ncols=7300, xmin=-4100000, xmax=3200000, ymin=1673000, ymax=6200000, crs = "EPSG:5072")

# Set output and download folder
out_f <- "./PredictionRasters/ABoVE"
if (!file.exists(out_f)) {
  dir.create(out_f, showWarnings = FALSE)
}
dwd.folder <- "./CovariateRasters/ABoVE"
if (!file.exists(dwd.folder)) {
  dir.create(dwd.folder, showWarnings = TRUE, recursive = TRUE)
}

## Access urls and extract yrlist
webData <- read.csv("./Scripts/webDataAccess.csv", fileEncoding="UTF-8-BOM")
above_lnd <- subset(webData, Dataset == "Landcover_ABoVE")
yrlist <- seq(as.numeric(strsplit(above_lnd$year, "[:]")[[1]][1]), as.numeric(strsplit(above_lnd$year, "[:]")[[1]][2]),by=1)
ABoVE_parentid <- above_lnd$url

########################################################
#Download ABoVE from GSDrive
dr <- drive_ls(as_id(ABoVE_parentid))
ABoVE_ls <- dr$name

for (i in 1:length(ABoVE_ls)){
  drive_download(as_id(dr$id[i]), path = file.path(dwd.folder, dr$name[i]))
}

setwd(dwd.folder)
r_files <- list.files(path = dwd.folder, pattern = ".tif")
# Read SpatRast
rasters <- lapply(file.path(dwd.folder,r_files), terra::rast)
#project
r_b1pj <- lapply(rasters, function(x) {project(x, rast1k, res = 1000, method = "near", align = TRUE)})
# Extend to NatMod extent
rast <- lapply(r_b1pj, function(x) {crop(x, rast1k, extend = TRUE)})
for (i in 1:length(r_files)) {names(rast[[i]]) <- paste0("lccABoVE1k_", substr(r_files[i],7,14))}
#save output
lapply(rast, function(x) {
  writeRaster(x, filename=file.path(out_f,names(x)), filetype="GTiff", overwrite=TRUE)
})




############################################################################
## Alternative source
############################################################################
## If TIFF is downloaded from source(https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=1691)
# Rasters have 31 band representing the years (1984-2014)
n_band<-c(1:31)
yearlist <- c(1984:2014)
download.file(url, destfile =file.path(dwd.folder,"Annual_Landcover_ABoVE_1691.zip"),  method="libcurl")
dwdname <- list.files(dwd.folder)
unzip(file.path(dwd.folder, dwdname), exdir = dwd.folder)
r_files <- list.files("./Annual_Landcover_ABoVE_1691/data", pattern = "ABoVE_LandCover_Bh", full.names = TRUE)
for (n in n_band) {
  #Create SpatRast
  rasters <- lapply(r_files, terra::rast)
  # select the band
  r_b1 <- lapply(rasters, function(x){x[[n]]})
  #Create SpatRastCollection
  raster_sprc <- terra::sprc(r_b1)
  # Mosaic rast
  raster_mosaic<- terra::mosaic(raster_sprc, fun = "max") 
  # Write just in case it crashed
  #writeRaster(raster_mosaic, filename=paste0("mosaic_b1_", as.character(yearlist[n]),".tif"), overwrite=TRUE)
  # Reproject
  r_b1pj <- project(raster_mosaic, rast1k, res = 1000, method = "near", align = TRUE)
  # Extend to NatMod extent
  rast <- terra::extend(r_b1pj, rast1k, fill = NA)
  #save output
  writeRaster(rast, filename=file.path(out.folder,paste0("lccABoVE1k_", as.character(yearlist[n]),".tif")), overwrite=TRUE)
}
  
  
  
  
  
  
  


