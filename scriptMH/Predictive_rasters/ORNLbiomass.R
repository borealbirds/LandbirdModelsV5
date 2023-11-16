################################################3
#  BAM NM 5.0 - Producing Prediction Rasters for biomass ONRL carbon density data for the year 2010
#  The script download from GEE, reproject and rescale to 1km
#  Ouptput raster is saved on disk and transfered to Google SharedDrive
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
setwd("E:/MelinaStuff/BAM/NationalModelv4.1")

# Set extent 
rast1km <- rast(nrows=4527, ncols=7300, xmin=-4100000, xmax=3200000, ymin=1673000, ymax=6200000, crs = "EPSG:5072")

# Set output and download folder
out_f <- "./PredictionRasters/ORNL"
if (!file.exists(out_f)) {
  dir.create(out_f, showWarnings = FALSE)
}

## Access urls and extract yrlist
webData <- read.csv("./Scripts/webDataAccess.csv", fileEncoding="UTF-8-BOM")
carbon_lnd <- subset(webData, Dataset == "Biomass carbon density")
img_name <- carbon_lnd$url

##########################################################
##Set function
##########################################################
##Set function
carbonDen <- function(x) {
  ccc <- ee$ImageCollection(x)$ # read the data
    select('agb')$ # select landcover
    toBands();# filter the year
  return(ccc)
}

##########################################################
# Extract image on GEE per year
imageextract <- carbonDen(img_name)
temp <- ee_as_raster(
    image = imageextract,
    container = "ONRL",
    dsn = "ONRL_1km_",
    scale = 300,
    #crs = 'EPSG:5072',
    via = "drive",
    maxPixels= 4e10
  )
#List downloaded images fo a specific year
rastlist <- list.files(path = "./", pattern="ONRL_1km_", 
                       all.files=TRUE, full.names=FALSE)

flf <-lapply(rastlist, rast) # Load images as SpatRaster
raster_sprc <- terra::sprc(flf) # Create SpatRastCollection
raster_mosaic<- terra::mosaic(raster_sprc, fun = "mean")  # Mosaic all SpatRast together
r_abgpj <- project(raster_mosaic, rast1km, res = 1000, method = "bilinear", align = TRUE) # Reproject
rast <- crop(r_abgpj, ext(rast1km), extend = TRUE) # Crop to NatMod extent

ONRLfocal5k<-terra::focal(r_abgpj,w=matrix(1,5,5),fun=mean,na.rm=TRUE) #get 5k average at 1k resolution
ONRL5k <- terra::crop(ONRLfocal5k, rast1k, extend = TRUE) # Extend to NatMod extent

#write output
writeRaster(rast, filename=file.path(out_f,"ONRLbiomass.tif"), filetype="GTiff", overwrite=TRUE)
writeRaster(ONRL5k, filename=file.path(out_f,"ONRLbiomass5k.tif"), filetype="GTiff", overwrite=TRUE)






# Convert to SPatRast (terra format)
spat_rast <- rast(temp)
# Reproject
spat_rast.pj <- project(spat_rast, rast1km, res = 1000, method = "bilinear", align = TRUE)
# Crop to NatMod extent
rast.cp <- crop(spat_rast.pj, ext(rast1km), extend = TRUE)
  
#write output
writeRaster(rast.cp, filename=file.path(out_f,paste0("RAP_TCC1km_",as.character(yrlist[i]), ".tif")), filetype="GTiff", overwrite=TRUE)
file.remove(rastlist)


#get 5k average at 1k resolution
rastlist <- list.files(path = out_f, pattern="RAP_TCC1km_", 
                       all.files=TRUE, full.names=FALSE)
rasters_1km <-lapply(file.path(out_f, rastlist), rast)
rasters5k<-lapply(rasters_1km, function(x) {focal(x,w=matrix(1,5,5),fun=mean,na.rm=TRUE)})
for (i in 1:length(rastlist)) {names(rasters5k[[i]]) <- paste0("RAP_TCC5x5_", substr(rastlist[i],12,15), ".tif")}
lapply(rasters5k, function(i) {
  writeRaster(i, filename=file.path(out_f,names(i)), filetype="GTiff", overwrite = TRUE)
})




  img_1k <- ee_as_raster(
     img,
      via = "drive",
      dsn = "ORNLbiomass.tif",
      #scale = 1000,
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
carbonDen.extract <- carbonDen(yr)
carbonDen1k <- rast(carbonDen.extract)
carbonDen_pj <- project(carbonDen1k, rast1k, res = 1000, method = "near", align = TRUE)
carbonDen_cp <-crop(carbonDen_pj, ext(rast1k), extend = TRUE)
writeRaster(carbonDen_cp, filename=file.path(out_f, "ORNLbiomass1km.tif"), overwrite=TRUE)

writeRaster(carbonDen_pj, filename=file.path(out_f, "ORNLbiomass1km_notcrop.tif"), overwrite=TRUE)
writeRaster(carbonDen1k, filename=file.path(out_f, "ORNLbiomass1km_rast.tif"), overwrite=TRUE)



