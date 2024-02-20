# ---
# title: National Models 4.1 - dataset extraction
# author: Elly Knight
# created: November 2, 2023
# ---

#NOTES################################

#This script provides portions of the national model dataset to BAM team members

#This script was initially written to provide Dan Yip with an approximate list of BAM database spatial coverage for BC so that he can seek out existing historic point count datasets for uploading to WildTrax and use in the National Models.

#PREAMBLE############################

#1. Load packages----

library(tidyverse) #basic data wrangling
library(sf)
library(gridExtra)

#2. Set root path for data on google drive----
root <- "G:/Shared drives/BAM_NationalModels/NationalModels5.0"

#3. Load data----
load(file=file.path(root, "Data", "04_NM5.0_data_stratify.Rdata"))

#BC COVERAGE FOR DAN#####

#1. Read in Canada shapefile----
can <- read_sf(file.path(root, "Regions", "CAN_adm", "CAN_adm1.shp"))
bc <- can %>% dplyr::filter(NAME_1=="British Columbia")

#2. Filter visits to BC----
visit.bc <- visit.bcr %>% 
  st_as_sf(coords=c("lon", "lat"), crs=4326, remove=FALSE) %>% 
  st_intersection(bc)

#3. Select only relevant columns and round lat lons----
visit.dan <- visit.bc %>% 
  st_drop_geometry() %>% 
  mutate(lat_round = round(lat, 1),
         lon_round = round(lon, 1)) %>% 
  group_by(source, organization, project, sensor, year, lat_round, lon_round) %>% 
  summarize(visits = n()) %>% 
  ungroup()

write.csv(visit.dan, file.path(root, "Data", "BAM_visits_BC_rounded.csv"), row.names = FALSE)

#BOREAL DATASET WITH OFFSETS FOR TATI####

#1. Filter visits----
visit.use <- visit %>% 
  dplyr::filter(source=="WildTrax") %>% 
  dplyr::select(source:jday)

#2. Filter bird and offset data by visits----
bird.use <- bird %>% 
  dplyr::filter(id %in% visit.use$id)

offsets.use <- offsets %>% 
  dplyr::filter(id %in% visit.use$id)

#3. Rename----
visit <- visit.use
bird <- bird.use
offset <- offsets.use

#4. Save----
save(visit, bird, offset, file="G:/Shared drives/BAM_Core/DataStuff/Data.Reqs/UBC-TatiMichelleti-Nov2023/BAM_NationalModelDataset_5_WildTraxOnly.R")
