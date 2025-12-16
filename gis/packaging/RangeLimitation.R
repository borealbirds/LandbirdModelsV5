###############################################
# title: "Range Limitation"
# author: Elly Knight
# date: December 12, 2025
##################################################

#This code uses presence/absence from BAM model training data, augmented by any complete ebird record up to 12hrs long and travelling up to 20km, to inform the probable maximal range extents of BAM-modeled species.

#PREAMBLE############################

#1. Load packages----
library(auk) #eBird data handling
library(tidyverse) #basic data wrangling
library(sf) #read in & handle study area
library(terra) #raster management
library(data.table) #other data wrangling
library(SpatialKDE) #KDE
library(dggridR) #grid for spatial thinning

#2. Set root path for data on google drive----
root <- "G:/Shared drives/BAM_NationalModels5"

#3. Load data package with WildTrax data ----
load(file.path(root, "data", "04_NM5.0_data_stratify.Rdata"))

#4. List of species ----
spp <- read.csv(file.path(root,"data","Lookups","lu_species.csv")) |> 
  dplyr::filter(species_code %in% colnames(birdlist[,-1])) |> 
  mutate(species_common_name = case_when(species_common_name=="House Wren" ~ "Northern House Wren",
                                         species_common_name=="Le Conte's Sparrow" ~ "LeConte's Sparrow",
                                         # species_common_name=="Yellow Warbler" ~ "Northern Yellow Warbler",
                                         # species_common_name=="Warbling Vireo" ~ "Eastern Warbling Vireo",
                                         !is.na(species_common_name) ~ species_common_name))

#GET EBIRD DATA################

#Only need to run this once

# #1. Get list of ebd objects to process----
# ebd.files <- list.files(file.path(root, "data", "eBird", "ebd_raw"), pattern="*relOct-2022")
# 
# #2. Set the sampling file path ----
# auk_sampling(file.path(root, "data", "eBird", "ebd_raw", "ebd_sampling_relJul-2024", "ebd_sampling_relJul-2024.txt"))
# 
# #2. Set up loop----
# for(i in 3:length(ebd.files)){
#   
#   #3. Set ebd path----
#   auk_set_ebd_path(file.path(root, "data", "eBird", "ebd_raw", ebd.files[i]), overwrite=TRUE)
#   
#   #4. Define filters----
#   filters <- auk_ebd(paste0(ebd.files[i], ".txt")) |> 
#     auk_species(c(spp$species_common_name)) |>
#     auk_protocol(c("Traveling", "Stationary")) |>
#     auk_year(c(1980, 2022)) |> # max/min of model training data
#     auk_date(c("*-05-15", "*-07-15")) |> #May 15-Jul 15
#     auk_duration(duration=c(0, 720)) |>   #<=12 hr long checklists
#     auk_complete()  # only complete checklists
#   
#   #5. Filter data----
#   #select columns to keep
#   filtered <- auk_filter(filters, file=file.path(root, "data", "eBird", "ebd_filtered_range", paste0(ebd.files[i], ".txt")), overwrite=TRUE)
#   
#   cat(i, " ")
#   
# }

#5. Check list of processed files----
ebd.files.done <- list.files(file.path(root, "data", "eBird", "ebd_filtered_range"), pattern="*.txt", full.names = TRUE)

#PUT DATA TOGETHER####

#1. Read the ebird data in ----
#Note this next line takes a long time to run (couple hours)
#test <- read_ebd(ebd.files.done[2])
raw.ebd.list <- purrr::map(.x=ebd.files.done, .f=~read_ebd(.))
raw.ebd <- do.call(rbind, raw.ebd.list)

#2. Wrangle the ebird data ----
dat.ebd <- raw.ebd |> 
  mutate(date_time = ymd_hms(paste0(observation_date, time_observations_started))) |> 
  dplyr::select(checklist_id, latitude, longitude, observation_count, date_time, common_name, observation_count) |> 
  rename(count = observation_count,
         id = checklist_id,
         lat = latitude,
         lon = longitude) |> 
  left_join(spp |> 
              rename(common_name = species_common_name) |> 
              dplyr::select(common_name, species_code)) |> 
  dplyr::select(-common_name)

#3. Wrangle the rest of the dataset ----
dat.bam <- as.matrix(bird) |> 
  as.data.frame() |> 
  rownames_to_column("id") |> 
  pivot_longer(-id, names_to = "species_code", values_to = "count") |> 
  mutate(id = as.numeric(id)) |> 
  left_join(visit |> 
              dplyr::select(id, lat, lon, date) |> 
              rename(date_time = date))

#4. Put together ----
dat <- rbind(dat.ebd |> 
               mutate(source = "eBird"),
             dat.bam |> 
               mutate(source = "BAM"))

rm(dat.bam, dat.ebd, raw.ebd)

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
  r <- rast(extent = ext(kde.i), resolution = 1000, crs = crs(kde.i))
  kde.r <- rasterize(kde.i, r, field="kde_value")
  
  #7. Binarize ----
  dat.min <- min(terra::extract(kde.r, vect(dat.thin))$kde_value)
  kde.1 <- kde.r
  values(kde.1) <- ifelse(values(kde.1 > dat.min), 1, NA)
  
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
