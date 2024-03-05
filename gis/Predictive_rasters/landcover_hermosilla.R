################################################3
#  BAM NAM 5.0- Producing Prediction Rasters for landcover Hermosilla data 
#  The script download from GEE, reproject and rescale to 1km
#  Temporal raster are saved on disk and transfered to Google SharedDrive
#################################################

#install.packages('terra', repos='https://rspatial.r-universe.dev')
library(terra)
library(sp)
library(rgee)
library(reticulate)

# Initialize Earth Engine
ee_check()
ee_Initialize()

# workdir
setwd("E:/MelinaStuff/BAM/NationalModelv5.0")

# Set extent 
rast1k <- rast(nrows=4527, ncols=7300, xmin=-4100000, xmax=3200000, ymin=1673000, ymax=6200000, crs = "EPSG:5072")

# Set output and download folder
out_f <- "./PredictionRasters/VLCE"
if (!file.exists(out_f)) {
  dir.create(out_f, showWarnings = FALSE)
}

## Access urls and extract yrlist
webData <- read.csv("./Scripts/webDataAccess.csv", fileEncoding="UTF-8-BOM")
hermo_lnd <- subset(webData, Dataset == "Landcover_hermosilla")
yrlist <- seq(as.numeric(strsplit(hermo_lnd$year, "[:]")[[1]][1]), as.numeric(strsplit(hermo_lnd$year, "[:]")[[1]][2]),by=1)
img_name <- hermo_lnd$url

##########################################################
##Set function
VLCE25 <- function(yr) {
  img_name <- paste0(img_name, as.character(yr))
  img <- ee$Image(img_name)$
    select('b1') 
  img_proj <- img$reproject(
    crs= "EPSG:5072", 
    scale= 1000
  )  
  img_1k <- ee_as_raster(
     img_proj,
      via = "drive",
      dsn = paste0("VLCE_", yr, ".tif"),
      scale = 1000,
      container = "rgee_backup",
      maxPixels = 4e+09,
      lazy = FALSE,
      public = TRUE,
      add_metadata = TRUE,
      timePrefix = TRUE,
      quiet = FALSE,
   )
  return(img_1k)
}

# Run function to extract
for (i in 1:length(yrlist)) {
  VLCE25extract <- VLCE25(yrlist[i])
  VLCE1k <- rast(VLCE25extract)
  VLCE_pj <- project(VLCE1k, rast1k, res = 1000, method = "near", align = TRUE)
  VLCE_cp <-crop(VLCE_pj, ext(rast1k), extend = TRUE)
  writeRaster(VLCE_cp, filename=file.path(out_f, paste0("VLCE1km_",as.character(yrlist[i]), ".tif")), overwrite=TRUE)
}



