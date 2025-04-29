# ---
# title: National Models 5.0 - stack and crop rasters for subunit prediction
# author: Elly Knight
# created: May 6, 2024
# ---

#NOTES################################

#This script uses the bcr shapefile and canada/USA shapefiles to make a polygon for each subunit with a 100 km buffer. We also do this for regions that we have merged.

#The output from this script is used in multiple stages of the analysis pipeline.

#Note the original version of this script also created a couple alternative subunits that were evaluated after models were run and removed from the script:
#merging 82 (Newfoundland) with mainland 81 and comparing to only 82 with small sample size
#merging 2 (coastal AK) with 41423 and comparing to only 2 with small sample size

#PREAMBLE############################

#1. Load packages----
library(tidyverse) #basic data wrangling
library(terra) #basic raster handling
library(sf) #basic shapefile handling
library(nngeo) #fill holes in polygons

#2. Set root path for data on google drive----
root <- "G:/Shared drives/BAM_NationalModels5"

#SUBUNIT POLYGONS#####################

#1. Read in country shapefiles----
can <- read_sf(file.path(root, "Regions", "CAN_adm", "CAN_adm0.shp")) |> 
  st_transform(crs=5072) 
usa <- read_sf(file.path(root, "Regions", "USA_adm", "USA_adm0.shp")) |> 
  st_transform(crs=5072) 

#2. Read in BCR shapefile----
#Remove subunit 1
bcr <- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel.shp")) |> 
  dplyr::filter(subUnit!=1) |> 
  st_transform(crs=5072)

ggplot(bcr) +
  geom_sf(aes(fill=factor(subUnit)))

#3. Identify BCRs for each country----
bcr.ca <- bcr |> 
  st_intersection(can) |> 
  mutate(country="can") |> 
  dplyr::select(subUnit, country)

bcr.usa <- bcr |> 
  st_intersection(usa) |> 
  mutate(country="usa") |> 
  dplyr::select(subUnit, country)

#4. Make merged subunits----

#For subunits that have < 1000 observations

#USA: merge 42 and 3 with 41

bcr.usa41423 <- bcr.usa |>
  dplyr::filter(subUnit %in% c(41, 42, 3)) |>
  st_union() |>
  nngeo::st_remove_holes() |>
  st_sf() |>
  mutate(country="usa",
         subUnit = 41423)
st_geometry(bcr.usa41423) <- "geometry"

#polygons to remove
bcr.remove <- data.frame(country=c("usa", "usa", "usa"),
                         subUnit=c(3, 42, 41))

#5. Put together----
bcr.country <- rbind(bcr.ca, bcr.usa, bcr.usa41423) |>
  anti_join(bcr.remove)

write_sf(bcr.country, file.path(root, "Regions", "BAM_BCR_NationalModel_Unbuffered.shp"))

#6. Set up loop for BCR buffering----
bcr.out <- data.frame()
for(i in 1:nrow(bcr.country)){

  #7. Filter & buffer shapefile----
  bcr.buff <- bcr.country |>
    dplyr::filter(row_number() == i) |>
    st_buffer(100000)

  #8. Crop to international boundary----
  if(bcr.buff$country=="can"){ bcr.i <- st_intersection(bcr.buff, can)}
  if(bcr.buff$country=="usa"){ bcr.i <- st_intersection(bcr.buff, usa)}

  #9. Put together----
  bcr.out <- rbind(bcr.out, bcr.i |>
                     dplyr::select(country, subUnit))

  print(paste0("Finished bcr ", i, " of ", nrow(bcr.country)))

}

#10. Save----
write_sf(bcr.out, file.path(root, "Regions", "BAM_BCR_NationalModel_Buffered.shp"))

#11. Plot----
ggplot(bcr.out)  +
  geom_sf(aes(fill=factor(subUnit)), alpha = 0.5)
