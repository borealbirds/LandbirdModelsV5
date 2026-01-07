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

#3. Load data ----
load(file.path(root, "data", "AllDetections.Rdata"))
load(file.path(root, "data", "04_NM5.0_data_stratify.Rdata"))
rm(bcrlist, bird, bootlist, cov, covlist, offsets, visit)

#4. List of species ----
spp <- read.csv(file.path(root,"data","Lookups","lu_species.csv")) |> 
  dplyr::filter(species_code %in% colnames(birdlist[,-1])) |> 
  mutate(species_common_name = case_when(species_common_name=="House Wren" ~ "Northern House Wren",
                                         species_common_name=="Le Conte's Sparrow" ~ "LeConte's Sparrow",
                                         !is.na(species_common_name) ~ species_common_name)) |> 
  arrange(species_code)

options(scipen=9999)

#5. Filter the days further ----
use <- dat |> 
  mutate(doy = yday(date_time)) |> 
  dplyr::filter(doy > 155,
                doy < 181)

#MAKE THE LIMITATION RASTERS########

#1. Create grid for thinning----
grid <- dgconstruct(spacing = 100, metric=TRUE)

#2. Set up loop ----
for(i in 1:nrow(spp)){
  
  spp.i <- spp$species_code[i]
  
  #3. Wrangle the data ----
  #reduce to just unique locations to smooth
  dat.i <- use |> 
    dplyr::filter(species_code==spp$species_code[i]) |> 
    mutate(lat = round(lat, 3),
           lon = round(lon, 3)) |> 
    dplyr::select(lat, lon) |>
    unique()
  
  #4. Thin ----
  suppressMessages(dat.thin <- dat.i |> 
    mutate(cell = dgGEO_to_SEQNUM(grid, lon, lat)$seqnum) |> 
    group_by(cell) |> 
    sample_n(1) |> 
    ungroup() |> 
    st_as_sf(coords=c("lon", "lat"), crs=4326) |> 
    st_transform(crs=5072))
  
  #5. KDE ----
  #we use a large band_width to get something smooth
  suppressMessages(kde.i <- kde(dat.thin, cell_size = 100000, band_width = 1000000))
  
  #6. Rasterize ----
  r <- rast(extent = ext(kde.i) + 20000, resolution = 10000, crs = crs(kde.i))
  kde.r <- rasterize(kde.i, r, field="kde_value")
  #plot(kde.r)
  
  #7. Binarize the raster ----
  #the isopleths below any observation
  dat.min <- quantile(terra::extract(kde.r, vect(dat.thin))$kde_value, 0)
  
  kde.0 <- kde.r
  values(kde.0) <- ifelse(values(kde.r) > dat.min, NA, 1)
  #plot(kde.0)
  
  kde.1 <- kde.r
  values(kde.1) <- ifelse(values(kde.r) > dat.min, 1, NA)
  #plot(kde.1)
  
  kde.01 <- kde.r
  values(kde.01) <- ifelse(values(kde.r) > dat.min, 1, 0)
  #plot(kde.01)
  
  #8. Make an exterior raster ----
  kde.edge <- kde.r
  values(kde.edge) <- ifelse(is.na(values(kde.edge)), 0, 1)
  #plot
  
  #9. Now make a distance raster ----
  kde.dist <- distance(kde.1)
  #plot(kde.dist)
  
  #10. Make an boundary raster ----
  nb.sum <- focal(kde.edge, w=matrix(1,3,3), fun=sum, na.policy="omit")
  edge <- (kde.edge==1) & (nb.sum < 9) & (nb.sum > 0)
  edge[edge==0] <- NA
  #plot(edge)
  
  #11. Convert to points ----
  pts <- as.points(edge, values=FALSE)
  #plot(kde.0, col=c("grey50", "grey90"))
  #plot(pts, add=TRUE, alpha=0.2)
  
  #13. Get the nearest edge point for each cell of the border raster ----
  kde.pts <- as.points(kde.0)
  kde.nearest <- nearest(x=kde.pts, y=pts, method="haversine")
  
  #14. Get the distances ----
  edge.distance <- terra::extract(x=kde.dist, y=pts)
  border.distance <- terra::extract(kde.dist, kde.pts)
  
  #14. Add it back to the full border raster ----
  kde.distance <- st_as_sf(kde.pts) |> 
    dplyr::select(geometry) |> 
    cbind(st_as_sf(kde.nearest) |> 
            st_drop_geometry()) |> 
    dplyr::select(from_id, to_id, geometry) |> 
    left_join(edge.distance |> 
                rename(to_id = ID,
                       maxdistance = kde_value),
              by="to_id") |> 
    left_join(border.distance |> 
                rename(from_id = ID,
                       distance = kde_value),
              by="from_id")
  #plot(kde.distance)
  
  #15. Scale ----
  kde.scale <- kde.distance |> 
    mutate(scale1 = distance/maxdistance,
           scale = 1-(scale1 - min(scale1))/(max(scale1) - min(scale1)))
  
  #16. Re-rasterize ----
  kde.scaler <- rasterize(kde.scale, r, field="scale", fun=mean)
  #plot(kde.scaler)
  
  #17. Put together ----
  kde.sum <- (kde.01 + ifel(is.na(kde.scaler), 0, kde.scaler))
  
  #18. Smooth ----
  kde.out <- focal(kde.sum, w=matrix(1,25,25), fun=median)
  names(kde.out) <- "range"
  #plot(kde.out)
  
  #18. Save ----
  writeRaster(kde.out, filename = file.path(root, "gis", "ranges", paste0(spp.i, ".tif")), overwrite=TRUE)
  
  #19. Plot for visualization ----
  kde.df <- kde.out |> 
    aggregate(10) |> 
    as.data.frame(xy=TRUE)
  
  ggplot(kde.df) +
    geom_raster(aes(x=x, y=y, fill=range), show.legend = FALSE) +
    scale_fill_viridis_c() +
    theme_minimal() +
    xlim(c(-6000000, 4500000)) +
    ylim(c(500000, 8000000))
  
  ggsave(filename = file.path(root, "gis", "ranges", "plots", paste0(spp.i, ".tif")), height = 4, width = 6, units="in")
  
  cat(i, " ")
  
}
