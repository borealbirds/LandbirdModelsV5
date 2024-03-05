################################################3
#  BAM NAM 5.0 - Producing Prediction Rasters for SCANFI biomass data 
#  The script download source data from the sharedDrive
#  Temporal raster are saved on disk and transferred to Google SharedDrive
#################################################

library(googledrive)
library(terra)
library(stringr)

########################
### PARAMS
########################
# workdir
setwd("E:/MelinaStuff/BAM/NationalModelv5.0")

# Set extent 
EPSG.5072 <- "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"
rast1k <- rast(nrows=4527, ncols=7300, xmin=-4100000, xmax=3200000, ymin=1673000, ymax=6200000, crs = EPSG.5072)

# output folder
dwd_f <- "./CovariateRasters/SCANFI_LCC"
if (!file.exists(dwd_f)) {
  dir.create(dwd_f, showWarnings = FALSE)
}
out_f <- "./PredictionRasters/SCANFI_LCC"
if (!file.exists(out_f)) {
  dir.create(out_f, recursive = TRUE, showWarnings = FALSE)
}
## Access urls and extract yrlist
webData <- read.csv("./Scripts/webDataAccess.csv", fileEncoding="UTF-8-BOM")
SCANFI_lnd <- subset(webData, Dataset == "SCANFI_LCC")
SCANFI_parentid <- SCANFI_lnd$url


##########################
##SCANFI Landcover download - Download from GDrive
dr <- drive_ls(as_id(SCANFI_parentid))

for (i in 1:nrow(dr)){
  drive_download(as_id(dr$id[i]), path = file.path(dwd_f, dr$name[i]))
}

#Generate predictive raster
lcc_list <- list.files(path = dwd_f, full.names = TRUE)
for (j in 1:length(lcc_list)){
  outfile <- basename(lcc_list[j])
  lcc <-rast(lcc_list[j])
  rast1k_pj <- project(rast1k, crs(lcc))
  lcc_res <- terra::resample(lcc, rast1k_pj, method = "near")
  lcc_pj <- project(lcc_res, rast1k, res = 1000, method = "near", align = TRUE)
  SCANFI.lcc <-crop(lcc_pj, ext(rast1k), extend = TRUE)
  writeRaster(SCANFI.lcc, file.path(out_f, paste0("SCANFIlcc_", str_sub(outfile, start = -11, end = -8),".tif")), overwrite=TRUE)
}

