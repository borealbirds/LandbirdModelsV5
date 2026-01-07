###############################################
# title: "Range Limitation - Data prep"
# author: Elly Knight
# date: December 12, 2025
##################################################

#This code collates presence/absence from BAM model training data, augmented by any complete ebird record up to 12hrs long and travelling up to 20km.
#PREAMBLE############################

#1. Load packages----
library(auk) #eBird data handling
library(tidyverse) #basic data wrangling
library(sf) #read in & handle study area
library(terra) #raster management
library(data.table) #other data wrangling

#2. Set root path for data on google drive----
root <- "G:/Shared drives/BAM_NationalModels5"

#3. Load data package with WildTrax data ----
load(file.path(root, "data", "04_NM5.0_data_stratify.Rdata"))

#4. List of species ----
spp <- read.csv(file.path(root,"data","Lookups","lu_species.csv")) |> 
  dplyr::filter(species_code %in% colnames(birdlist[,-1])) |> 
  mutate(species_common_name = case_when(species_common_name=="House Wren" ~ "Northern House Wren",
                                         species_common_name=="Le Conte's Sparrow" ~ "LeConte's Sparrow",
                                         !is.na(species_common_name) ~ species_common_name))

#GET EBIRD DATA################

#Only need to run this section once

#1. Get list of ebd objects to process----
ebd.files <- list.files(file.path(root, "data", "eBird", "ebd_raw"), pattern="*relOct-2022")

#2. Set the sampling file path ----
auk_sampling(file.path(root, "data", "eBird", "ebd_raw", "ebd_sampling_relJul-2024", "ebd_sampling_relJul-2024.txt"))

#2. Set up loop----
for(i in 3:length(ebd.files)){

  #3. Set ebd path----
  auk_set_ebd_path(file.path(root, "data", "eBird", "ebd_raw", ebd.files[i]), overwrite=TRUE)

  #4. Define filters----
  filters <- auk_ebd(paste0(ebd.files[i], ".txt")) |>
    auk_species(c(spp$species_common_name)) |>
    auk_protocol(c("Traveling", "Stationary")) |>
    auk_year(c(1980, 2022)) |> # max/min of model training data
    auk_date(c("*-05-15", "*-07-15")) |> #May 15-Jul 15
    auk_duration(duration=c(0, 720)) |>   #<=12 hr long checklists
    auk_complete()  # only complete checklists

  #5. Filter data----
  #select columns to keep
  filtered <- auk_filter(filters, file=file.path(root, "data", "eBird", "ebd_filtered_range", paste0(ebd.files[i], ".txt")), overwrite=TRUE)

  cat(i, " ")

}

#PUT DATA TOGETHER####

#1. Check list of processed files----
ebd.files.done <- list.files(file.path(root, "data", "eBird", "ebd_filtered_range"), pattern="*.txt", full.names = TRUE)

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
  dplyr::filter(count > 0) |> 
  left_join(visit |> 
              dplyr::select(id, lat, lon, date) |> 
              rename(date_time = date))

#4. Put together ----
dat <- rbind(dat.ebd |> 
               mutate(source = "eBird"),
             dat.bam |> 
               mutate(source = "BAM"))

save(dat, file = file.path(root, "data", "AllDetections.Rdata"))
