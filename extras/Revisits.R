# ---
# title: National Models 4.1 - revisits
# author: Elly Knight
# created: November 2, 2023
# ---

#NOTES################################

#This script identifies revisits between years in the BAM national model dataset

#PREAMBLE############################

#1. Load packages----

library(tidyverse) #basic data wrangling
library(sf)
library(gridExtra)

#2. Set root path for data on google drive----
root <- "G:/Shared drives/BAM_NationalModels/NationalModels5.0"

#3. Load data----
load(file=file.path(root, "Data", "03_NM5.0_data_stratify.R"))

#FILTERING##############

#1. Read in Canada shapefile----
can <- read_sf(file.path(root, "Regions", "CAN_adm", "CAN_adm1.shp"))
ab <- can %>% dplyr::filter(NAME_1=="Alberta")

#2. Filter visits to AB----
visit.ab <- visit.bcr %>% 
  st_as_sf(coords=c("lon", "lat"), crs=4326, remove=FALSE) %>% 
  st_intersection(ab)

#3. Filter to locations & project----
loc.ab <- visit.ab %>% 
  dplyr::select(organization, project, location, lat, lon) %>% 
  unique() %>% 
  st_as_sf(coords=c("lon", "lat"), crs=4326, remove=FALSE) %>% 
  st_transform(3403)

#4. Buffer by 30 m----
loc.buff <- loc.ab %>% 
  st_buffer(30)

#5. Identify points within buffers----
loc.int <- st_intersects(loc.buff, loc.ab)

#6. Assign group ids----
loc.ab$group <- NA
for(i in 1:nrow(loc.ab)){
  
  loc.i <- loc.int[[i]]
  
  for(j in 1:length(loc.i)){
    
    loc.j <- loc.i[[j]]
    
    if(is.na(loc.ab$group[loc.j])){
      
      loc.ab$group[loc.j] <- i
      
    }

  }
}

#7. Rejoin to visits----
visit.out <- visit.ab %>% 
  st_drop_geometry() %>% 
  left_join(loc.ab) %>% 
  dplyr::select(source, organization, project, sensor, tagMethod, location, lat, lon, group, year, date, observer, duration, distance) %>% 
  rename(revisit_group = group)

write.csv(visit.out, file.path(root, "Data", "BAM_Revisits.csv"), row.names = FALSE)

#8. Summarize by year----
year.out <- visit.out %>% 
  dplyr::select(source, organization, project, sensor, location, lat, lon, revisit_group, year) %>% 
  unique() %>% 
  group_by(source, organization, project, sensor, location, lat, lon, revisit_group) %>%
  summarize(years = n()) %>% 
  ungroup() %>% 
  dplyr::filter(years > 1) %>% 
  left_join(visit.out %>% 
              dplyr::select(source, organization, project, sensor, location, lat, lon, revisit_group, year) %>% 
              unique(),
            multiple="all")
