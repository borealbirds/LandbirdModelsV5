###############################################
# title: "Range Limitation"
# author: Elly Knight
# date: December 12, 2025
##################################################

#This code uses presence/absence from BAM model training data, augmented by any complete ebird record up to 12hrs long and travelling up to 20km, to inform the probable maximal range extents of BAM-modeled species.

#PREAMBLE############################

#1. Load packages----
library(tidyverse) #basic data wrangling
library(sf) #read in & handle study area
library(terra) #raster management
library(data.table) #other data wrangling
library(SpatialKDE) #KDE
library(dggridR) #grid for spatial thinning

#2. Set root path for data on google drive----
root <- "G:/Shared drives/BAM_NationalModels5"

#3. Load data package with WildTrax data ----
load(file.path(root, "data", "AllDetections.Rdata"))

#4. List of species ----
spp <- read.csv(file.path(root,"data","Lookups","lu_species.csv")) |> 
  dplyr::filter(species_code %in% colnames(birdlist[,-1])) |> 
  mutate(species_common_name = case_when(species_common_name=="House Wren" ~ "Northern House Wren",
                                         species_common_name=="Le Conte's Sparrow" ~ "LeConte's Sparrow",
                                         # species_common_name=="Yellow Warbler" ~ "Northern Yellow Warbler",
                                         # species_common_name=="Warbling Vireo" ~ "Eastern Warbling Vireo",
                                         !is.na(species_common_name) ~ species_common_name))

#MAKE THE LIMITATION RASTERS########

#1. Create grid for thinning----
grid <- dgconstruct(spacing = 100, metric=TRUE)

#2. Set up loop ----
for(i in 1:nrow(spp)){
  
  spp.i <- spp$species_code[i]
  
  #3. Wrangle the data ----
  #reduce to just unique locations to smooth
  dat.i <- dat |> 
    dplyr::filter(species_code==spp[i]) |> 
    mutate(lat = round(lat, 3),
           lon = round(lon, 3)) |> 
    dplyr::select(lat, lon) |>
    unique()
  
  #4. Thin ----
  dat.thin <- dat.i |> 
    mutate(cell = dgGEO_to_SEQNUM(grid, lon, lat)$seqnum) |> 
    group_by(cell) |> 
    sample_n(1) |> 
    ungroup() |> 
    st_as_sf(coords=c("lon", "lat"), crs=4326) |> 
    st_transform(crs=5072)
  
  #5. KDE ----
  #we use a large band_width to get something smooth
  kde.i <- kde(dat.thin, cell_size = 100000, band_width = 500000)
  
  #6. Rasterize ----
  #let's use 100 km cells for now
  r <- rast(extent = ext(kde.i) + 5000, resolution = 10000, crs = crs(kde.i))
  kde.r <- rasterize(kde.i, r, field="kde_value")
  
  #7. Make an edge raster ----
  kde.edge <- kde.r
  values(kde.edge) <- ifelse(is.na(values(kde.r)), 1, 0)
  nb.sum <- focal(kde.edge, w= matrix(1,3,3), fun=sum, na.policy="omit")
  edge <- (kde.edge==1) & (nb.sum < 9)
  edge[edge==0] <- NA
  plot(edge)
  
  #9. Convert to points ----
  pts <- as.points(edge, values=FALSE)
  plot(pts)
  
  #10. Binarize the raster ----
  #the isopleths below the percentile of any observation
  dat.min <- min(terra::extract(kde.r, vect(dat.thin))$kde_value)
  
  kde.0 <- kde.r
  values(kde.0) <- ifelse(values(kde.r > dat.min), NA, 1)
  plot(kde.0)
  
  kde.1 <- kde.r
  values(kde.1) <- ifelse(values(kde.r > dat.min), 1, NA)
  plot(kde.1)
  
  #11. Now make a distance raster ----
  kde.dist <- distance(kde.1)
  plot(kde.dist)
  
  #12. Get the distance for each point ----
  pts.dist <- terra::extract(x=kde.dist, y=pts)
  
  #13. Get the nearest point for each cell of the raster
  kde.pts <- as.points(kde.0)
  kde.nearest <- nearest(kde.pts, pts)
  
  
  #7. Binarize ----
0

  
  #8. Edge of raster ----
  kde.edge <- rast(kde.10)
  kde.edge[] <- NA
  nr <- nrow(kde.edge)
  nc <- ncol(kde.edge)
  kde.edge[1, 1:nc] <- 1
  kde.edge[nr, 1:nc] <- 1
  kde.edge[1:nr, 1] <- 1
  kde.edge[1:nr, nc] <- 1
  

  

  
  #10. Make a raster of just the boundary points ----
  kde.0 <- kde.r 
  values(kde.0) <- ifelse(values(kde.0 > dat.min), NA, 1)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

  
  kde.0 <- kde.r
  values(kde.0) <- ifelse(values(kde.0 > dat.min), NA, values(kde.0))
  
  #8. Rescale ----
  kde.scale <- (kde.0 - minmax(kde.0)[1])/(minmax(kde.0)[2] - minmax(kde.0)[1])
  
  #8. Make a shp of the outside ----
  kde.NA <- kde.scale
  values(kde.NA) <- ifelse(values(is.na(kde.NA)), 1, NA)
  shp.NA <- as.polygons(kde.NA, dissolve=TRUE)
  
  #8. Create distance-based decay of 0s ----
  kde.dist <- distance(kde.1, shp.NA)
  kde.dist1 <- distance(kde.1) |> 
    mask(kde.r)
  kde.distNA <- distance(kde.NA) |> 
    mask(kde.r)
  kde.dist <- kde.dist1/(kde.dist1 + kde.dist0)
  plot(kde.dist)
  
  #9. Put together ----
  kde.mask <- kde.dist
  values(kde.mask) <- ifelse(values(is.na(kde.1)), values(kde.r), values(kde.r))
  
  #10. Rescale ----
  kde.out <- kde.mask
  kde.out <- 1 - (kde.mask - minmax(kde.mask)[1])/(minmax(kde.mask)[2] - minmax(kde.mask)[1])
  
  
  #7. Scale to 100 ----
  kde.100 <- (kde.r - minmax(kde.r)[1])/(minmax(kde.r)[2] - minmax(kde.r)[1])
  
  
  #check
  plot(kde.100)
  plot(dat.thin, add=TRUE, colour="black")
  dat.r <- terra::extract(kde.100, vect(dat.thin))
  
  #7. Convert to a mask ----
  kde.mask <- kde.100
  values(kde.mask) <- ifelse(values(kde.mask > 0.05), NA, 
                             values(kde.mask))
  
  
  
}
