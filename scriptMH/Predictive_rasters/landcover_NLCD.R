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
dwd_f <- "./CovariateRasters/Landco_NLCD"
if (!file.exists(dwd_f)) {
  dir.create(dwd_f, showWarnings = FALSE)
}
out_f <- "./PredictionRasters/Landco_NLCD"
if (!file.exists(out_f)) {
  dir.create(out_f, recursive = TRUE, showWarnings = FALSE)
}
## Access urls and extract yrlist
webData <- read.csv("./Scripts/webDataAccess.csv", fileEncoding="UTF-8-BOM")
NLCD_lnd <- subset(webData, Dataset == "Landco_NLCD")
NLCD_parentid <- NLCD_lnd$url


##########################
##SCANFI Landcover download - Download from GDrive
dr <- drive_ls(as_id(NLCD_parentid))
# filter out index and non raster image
dr <- dr %>%
  dplyr::filter(!grepl("index", dr$name))

for (i in 1:nrow(dr)){
  drive_download(as_id(dr$id[i]), path = file.path(dwd_f, dr$name[i]))
}


#Generate predictive raster
nlcd_list <- list.files(dwd_f, full.names = TRUE, pattern="*.img")
for (j in 1:length(nlcd_list)){
  outfile <- basename(nlcd_list[j])
  year <- str_sub(outfile, start = 6, end = 9)
  nlcd <-rast(nlcd_list[j])
  rast1k_pj <- project(rast1k, crs(nlcd))
  NLCD.crop <-crop(nlcd, rast1k_pj)
  NLCD.res <- resample(NLCD.crop, rast1k_pj, method = "near")
  NLCD.pj <- project(NLCD.res, rast1k, res = 1000, method = "near", align = TRUE)
  NLCD.lcc <-crop(NLCD.pj, ext(rast1k), extend = TRUE)
  writeRaster(NLCD.lcc, file.path(out_f, paste0("NLCDlcc_", year,".tif")), overwrite=TRUE)
}
