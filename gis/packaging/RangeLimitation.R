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
library(classInt) #jenks break for outliers

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

#6. Get the BCRs for plotting ----
bcr <- read_sf(file.path(root, "gis", "Subregions_unbuffered.shp")) |> 
  st_simplify(dTolerance = units::set_units(100, "km")) |> 
  st_transform(crs=3978) |> 
  st_make_valid()

#7. Union for cropping ----
sa <- st_union(bcr)

#8. Species that we dont truncate for outliers ----
#determined by visual inspection, typically generalist species with wide distributions that we don't think have vagrants, or high north species with sparse data points
# spp0 <- c("ALFL", "AMCR", "ATTW", "BCCH", "BHCO", "BRCR", "CAJA", "CEDW", "CORA", "DEJU", "DOWO", "DUFL", "DUNL", "GCKI", "HAFL", "HAWO", "HETH", "HOLA", "LISP", "MAWR", "MOBL", "MODO", "NESP", "NOFL", "NOWA", "OCWA", "OSFL", "OVEN", "PHVI", "PISI", "PIWO", "PUFI", "RBNU", "RCKI", "REVI", "RUBL", "RUGR", "RWBL", "SAVS", "SOGR", "SOSA", "SOSP", "SPSA", "SWTH", "TEWA", "TRES", "VATH", "VESP", "WAVI", "WEWP", "WIPT", "WISN", "WWCR", "YBFL", "YEWA", "YRWA")

spp2 <- c("AMPI", "ATSP", "ATTW", "BANS", "BARS", "BAWW", "BBCU", "BBWA", "BBWO", "BGGN", "BHVI", "BLBW", "BLPW", "BOWA", "BRBL", "BRCR", "BRTH", "BTBW", "BTNW", "BWWA", "CCSP", "CLSW", "CMWA", "COGR", "CONW", "CSWA", "DUFL", "EABL", "EATO", "EAWP", "FISP", "FOSP", "GCSP", "GCTH", "GRSP", "GRYE", "GWWA", "HOLA", "INBU", "LALO", "LCSP", "LEYE", "MAWA", "MOBL", "MOWA",  "NAWA", "NOCA", "NOPA", "PAWA", "PHVI", "PIGR", "PIWA", "RBGR", "RBWO", "RECR", "RHWO", "SCTA", "SEWR", "TOSO", "TOWA", "UPSA", "WCSP", "WETA", "WITU", "WIWA", "YBCU", "YHBL")

#MAKE THE LIMITATION RASTERS########

#1. Create grid for thinning----
grid <- dgconstruct(spacing = 10, metric=TRUE)

#2. Set up loop ----
for(i in 1:nrow(spp)){
  
  spp.i <- spp$species_code[i]
  
  #3. Wrangle the data ----
  #reduce to just unique locations to smooth
  dat.i <- use |> 
    dplyr::filter(species_code==spp.i) |> 
    mutate(lat = round(lat, 3),
           lon = round(lon, 3)) |> 
    dplyr::select(lat, lon) |>
    unique()
  
  #4. Thin ----
  suppressMessages(dat.thin1 <- dat.i |> 
    mutate(cell = dgGEO_to_SEQNUM(grid, lon, lat)$seqnum) |> 
    group_by(cell) |> 
    sample_n(1) |> 
    ungroup() |> 
    st_as_sf(coords=c("lon", "lat"), crs=4326, remove=FALSE) |> 
    st_transform(crs=3978) |> 
    st_intersection(bcr))
  
  #Manually remove west BWWA detections that are probably supposed to be a different sppp
  if(spp.i %in% c("BWWA", "BTBW", "NOPA")){
    dat.thin1 <- dplyr::filter(dat.thin1, lon > -100)
  }
  
  #5. KDE ----
  #we use a large band_width to get something smooth
  suppressMessages(kde1.i <- kde(dat.thin1, cell_size = 70000, band_width = 700000))

  #6. Rasterize ----
  r1 <- rast(extent = ext(kde1.i) + 20000, resolution = 10000, crs = crs(kde1.i))
  kde.r1 <- rasterize(kde1.i, r1, field="kde_value")
  
  #7. Get the isopleths ----
  iso <- terra::extract(kde.r1, vect(dat.thin1))$kde_value

  #7. Binarize the raster to remove outliers ----
  #the isopleths below the 1% quantile of points with data unless is a 0 species
  if(!spp.i %in% spp2){
    dat.q <- quantile(terra::extract(kde.r1, vect(dat.thin1))$kde_value, 0)
  } else {
    dat.q <- quantile(terra::extract(kde.r1, vect(dat.thin1))$kde_value, 0.01)
  }

  kde.q <- kde.r1
  values(kde.q) <- ifelse(values(kde.r1) > dat.q, 1, NA)
  #plot(kde.q)

  #8. Remove observations below that
  dat.thin <- terra::extract(kde.q, vect(dat.thin1), bind=TRUE) |>
    st_as_sf() |>
    dplyr::filter(!is.na(kde_value))
  
  # ggplot() +
  #   geom_point(data=dat.thin1, aes(x=lon, y=lat), colour="black") +
  #   geom_point(data=dat.thin, aes(x=lon, y=lat), colour="grey")
    
  #9. Redo the KDE ----
  suppressMessages(kde.i <- kde(dat.thin, cell_size = 70000, band_width = 700000))
  
  #10. Rasterize again ----
  r <- rast(extent = ext(kde.i) + 20000, resolution = 10000, crs = crs(kde.i))
  kde.r <- rasterize(kde.i, r, field="kde_value")
  # plot(kde.r1)
  # plot(kde.r)
  
  #11. Get new min quantile ----
  dat.min <- quantile(terra::extract(kde.r, vect(dat.thin))$kde_value, 0)
  
  #11. Binarize ----
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
  #plot(kde.edge)
  
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
           scale = 1-(scale1 - min(scale1))/(max(scale1) - min(scale1)),
           scale2 = ifelse(scale < 0.5, NA, (scale - 0.5)/(0.5)))
  
  #16. Take 50% and rescale ----
  kde.scale2 <- 
  
  #16. Re-rasterize ----
  kde.scaler <- rasterize(kde.scale, r, field="scale2", fun=mean)
  #plot(kde.scaler)
  
  #17. Put together ----
  kde.sum <- (kde.01 + ifel(is.na(kde.scaler), 0, kde.scaler))
  kde.sum2 <- ifel(kde.sum==0, NA, kde.sum) 
  
  #18. Smooth ----
  kde.out <- focal(kde.sum2, w=matrix(1,25,25), fun=median) |> 
    project("EPSG:3978") |> 
    crop(bcr, mask=TRUE)
  names(kde.out) <- "range"
  #plot(kde.out)
  
  #19. Plot for visualization ----
  kde.df <- kde.out |> 
    aggregate(2) |> 
    as.data.frame(xy=TRUE)
  
  ggplot(kde.df) +
    geom_tile(aes(x=x, y=y, fill=range), show.legend = FALSE) +
    geom_sf(data=dat.thin1, colour="black") +
    geom_sf(data=dat.thin, colour="grey") +
    geom_sf(data=bcr, fill=NA, colour="grey30", linewidth=0.2) +
    scale_fill_viridis_c() +
    theme_minimal()
  
  ggsave(filename = file.path(root, "gis", "ranges", "plots", paste0(spp.i, ".tif")), height = 4, width = 6, units="in")
  
  #18. Save ----
  writeRaster(kde.out, filename = file.path(root, "gis", "ranges", paste0(spp.i, ".tif")), overwrite=TRUE)
  
  cat(i, " ")
  
}
