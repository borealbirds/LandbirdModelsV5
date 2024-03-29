################################
### Feathering BCR Rasters  ####
### Anna Drake, March 2024  ####
################################
# Uses approach by Dominic Roye for distance calculation (dominicroye.github.io)

# Process calculates the distance of each grid to the edge of the prediction area - divides by max to standardize
# Creates weight raster with decreasing weight toward the edge and writes out currently
# Next steps:
# Clip out the extrapolated areas at the bcr level <- set to 0
# Create Canada-wide weight raster by mosaicing and *summing* the weight rasters
# Final mosaic:
# Multiply the prediction raster by the distance raster (eliminates areas where extrapolation occurs, downweighted edge)
# Merge raster stacks using *sums*
# Divide the final raster by the Canada-wide weight raster (this undoes the initial edge weighting where there is no overlap of rasters)

#packages
library(sf)
library(tidyverse)
library(RColorBrewer)
library(terra)

#1. Preamble
# NAD83(NSRS2007)/Conus Albers projection (epsg:5072)
crs<-"+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"
root<-"G:/Shared drives/BAM_NationalModels/NationalModels5.0" #PC link

#2. Running US BCRs - Get BCR boundaries ----
usa <- read_sf(file.path(root, "Regions", "USA_adm", "USA_adm0.shp")) %>% 
  st_transform(crs=crs) # US boundary

bcr.us <- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel.shp")) %>% 
  dplyr::filter(subUnit!=1) %>% # remove BCR 1
  st_transform(crs=crs) %>%
  st_intersection(usa) %>% 
  mutate(country="us") # Restrict to USA BCRs

#3. Produce weighted distance-to-edge raster

# Extract and buffer BCR
for (i in (1:nrow(bcr.us))) {
  bcr.i <-  bcr.us %>% 
    dplyr::filter(row_number()==i) %>% 
    st_buffer(100000)%>%  # Filter & buffer shapefile
    st_intersection(., usa) #crop at USA border

# Make fishnet (10km2 to speed processing - could make finer)
  grid <- st_make_grid(bcr.i, cellsize = 10000, what = "centers")
  grid <- st_intersection(grid, bcr.i)  
  plot(grid)
  # Polygon -> line
  bcr.i <- st_cast(bcr.i, "MULTILINESTRING")

# Calculate distance from the edge for each point
  dist <- st_distance(bcr.i, grid)

# Turn into dataframe as proportion of max distance
  df <- data.frame(dist = as.vector(dist)/as.vector(max(dist)),
                   st_coordinates(grid))

# Rasterize, match resolution
  ext <- ext(bcr.i)
  r <- terra::rast(resolution = 10000, ext = ext, crs = crs)
  dist_sf <- st_as_sf(df, coords = c("X", "Y")) %>% st_set_crs(crs)
  dist_raster <- terra::rasterize(dist_sf, r, "dist", fun = mean) %>% terra::disagg(.,fact=10, method="bilinear")
  writeRaster(dist_raster,paste("EdgeWeighting_USA_BCR_", bcr.i$subUnit, ".tif", sep=""))
} 


# Import extrapolation Rasters 
  
root<-"G:/Shared drives/BAM_NationalModels/NationalModels5.0/Gap Analysis/MaskingRasters"

f <- lapply(list.files(root, paste("BCR_",bcr.i$subUnit,sep=""), full.names = TRUE),
              terra::rast)  

# Need to finalize extrapolation rasters before masking
