# ---
# title: National Models 4.1 - presentation statistics & figures
# author: Elly Knight
# created: November 2, 2023
# ---

#NOTES################################

#This script provides portions of the national model dataset to BAM team members

#This script was initially written to provide Dan Yip with an approximate list of BAM database spatial coverage for BC.

#PREAMBLE############################

#1. Load packages----

library(tidyverse) #basic data wrangling
library(sf)
library(gridExtra)

#2. Set root path for data on google drive----
root <- "G:/Shared drives/BAM_NationalModels/NationalModels4.1"

#3. Load data----
load(file=file.path(root, "Data", "03_NM4.1_data_stratify.R"))

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
