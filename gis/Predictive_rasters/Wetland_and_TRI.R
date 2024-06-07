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

# Template layer
r<-terra::rast(file.path(root,"PredictionRasters","Topography","mTPI_1km.tif"))

### Wetland layers are taken down from GEE at 350m resolution (e.g.):
#var occur = ee.Image('JRC/GSW1_4/GlobalSurfaceWater')
#.select(['occurrence']);  
#var region = ee.Geometry.BBox(-169.57, 71.22,-50.39, 38.09);
#Export.image.toDrive({image: occur, description: 'Wetland_Occur',region: region,
#scale: 350, maxPixels: 1e13,crs: 'EPSG:5072'});
# = ~ area of 200m radius circle, replicating covariate extraction scale

# Occurrence

occur<-terra::rast(file.path("G:/My Drive/Wetland","Wetland_Occur.tif"))
occur_1k<-project(occur,r, method="near")

## Get 5k scaling up from an averaged 1km resolution (as was done with Peatland)

occur_5k<-occur%>%project(.,r)%>%focal(., w=matrix(1,5,5), fun=mean, na.rm=T)

# Recurrence and Seasonality

recur<-terra::rast(file.path("G:/My Drive/Wetland","Wetland_recur.tif"))
recur_1k<-project(recur,r, method="near")

season<-terra::rast(file.path("G:/My Drive/Wetland","Wetland_seasonality.tif"))
season_1k<-project(season,r, method="near")

writeRaster(occur_1k,"occur_1kmFINAL.tif")
writeRaster(occur_5k,"occur_5kmFINAL.tif")
writeRaster(recur_1k,"recur_1kmFINAL.tif")
writeRaster(recur_1k,"recur_1kmFINAL.tif")


# Roughness Topography - original BAM layer produced from AdaptWest DEM

roughness<-terra::rast(file.path(root,"CovariateRasters","Topography","roughness.tif"))

# Match the covariate sampling radius using focal then use nearest neighbor

gf<-focalWeight(roughness, 200, "circle", fillNA=TRUE) #~200m radius circle
roughness <- focal(roughness, w=gf, fun=mean, na.rm=T)
roughness_test1km<-project(roughness,r, method="near") #1km
writeRaster(roughness_test1km,"roughness_focalweightresampled.tif")
