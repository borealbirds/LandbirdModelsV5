################################################3
#  BAM NM 5.0 - Producing Prediction Rasters JRC Global Surface Water
#  Correcting 1km terrain roughness layer (coded "TRI" in models)
#  Wetland data downloaded from GEE to disk
#  Saved on disk and transfered to Google SharedDrive
#  Author Anna Drake
#################################################

require(raster)
require(terra)
require(dplyr)
require(sf)

root <- "G:/Shared drives/BAM_NationalModels5"

# Template layer
r<-terra::rast(file.path(root,"PredictionRasters","Topography","mTPI_1km.tif"))

### Wetland layers are taken down from GEE at 350m resolution (e.g.):
#var occur = ee.Image('JRC/GSW1_4/GlobalSurfaceWater')
#.select(['occurrence']);  
#var region = ee.Geometry.BBox(-169.57, 71.22,-50.39, 38.09);
#Export.image.toDrive({image: occur, description: 'Wetland_Occur',region: region,
#scale: 350, maxPixels: 1e13,crs: 'EPSG:5072'});
# = ~ area of 200m radius circle, replicating covariate extraction scale
# Correction - GEE appears to adjust resolution on these variables using
# it's pyramid approach to categorical variables (likely due to these being temporal
# not spatial. 350 therefore appears to be a nearest-neighbour sample, not a mean.
# Alignment and proportional representation with the covariate values seems acceptable however
# To match 5x5 covariate values in occurrence, image was exported at 5000m scale and 
# nearest neighbour sampling was used to bring it back to 1km scale

# Occurrence

occur<-terra::rast(file.path("G:/My Drive/Wetland","Wetland_Occur.tif"))
occur_1k<-project(occur,r, method="near")

# Recurrence and Seasonality

recur<-terra::rast(file.path("G:/My Drive/Wetland","Wetland_recur.tif"))
recur_1k<-project(recur,r, method="near")

season<-terra::rast(file.path("G:/My Drive/Wetland","Wetland_seasonality.tif"))
season_1k<-project(season,r, method="near")

writeRaster(occur_1k,"occur_1kmFINAL.tif")
writeRaster(occur_5k,"occur_5kmFINAL.tif")
writeRaster(recur_1k,"recur_1kmFINAL.tif")
writeRaster(season_1k,"seasonality_1kmFINAL.tif")

# Roughness Topography - original BAM layer produced from AdaptWest DEM

roughness<-terra::rast(file.path(root,"CovariateRasters","Topography","roughness.tif"))

# Match the covariate sampling radius using focal then use nearest neighbor

gf<-focalWeight(roughness, 200, "circle", fillNA=TRUE) #~200m radius circle
roughness <- focal(roughness, w=gf, fun=mean, na.rm=T)
roughness_test1km<-project(roughness,r, method="near") #1km
writeRaster(roughness_test1km,"roughness_focalweightresampled.tif")

# LED - original layer by BEACONS, requires downsampling only

led <- terra::rast(file.path(root, "CovariateRasters", "Wetland", "led_250.tif"))
led_1k <- project(led, r, method="near")
writeRaster(led_1k, file.path(root, "PredictionRasters", "Wetland", "LakeEdge_1km.tif"))
