# ---
# title: National Models 4.1 - calculate offsets
# author: Elly Knight
# created: December 22, 2022
# ---

#NOTES################################

# This script uses the qpad-offsets package, which requires downloading that R project from github (https://github.com/borealbirds/qpad-offsets) to your local and setting the working directory to it (line 29).

#PREAMBLE############################

#1. Load packages----

library(tidyverse) #basic data wrangling
library(lubridate) #temporal data wrangling
library(QPAD) #to get offset model estimates
library(raster) #to handle covariate extraction from qpad rasters
library(maptools) #required for QPAD
library(intrval) #required for QPAD
library(data.table) #collapse list to dataframe

#2. Set root path for data on google drive----
root <- "G:/.shortcut-targets-by-id/0B1zm_qsix-gPbkpkNGxvaXV0RmM/BAM.SharedDrive/RshProjs/PopnStatus/NationalModelsV4.1/PointCount/"

#A. LOAD QPAD#####

#1. Set WD to qpad-offsets package----
setwd("C:/Users/Elly Knight/Documents/BAM/Projects/QPAD/qpad-offsets")

#2. Load QPAD requirements----
#2a. Estimates----
load_BAM_QPAD(version = 3)

#2b. Raster data----
rlcc <- raster("./data/lcc.tif")
rtree <- raster("./data/tree.tif")
rtz <- raster("./data/utcoffset.tif")
rd1 <- raster("./data/seedgrow.tif")
crs <- proj4string(rtree)

#2c. Source functions----
source("functions.R")

#B. PREP DATA####

#1. Load data package from script 01----
load(file.path(root, "01_NM4.1_data_clean.R"))

#2. Select columns----
#add column for whether is local or utc time
bird.x <- bird %>% 
  mutate(dt = as.character(date(date)),
         tm = as.character(paste0(hour(date), ":", minute(date)))) %>% 
  rename(dur = duration, dis = distance) %>% 
  dplyr::filter(lon >= -164,
                lon <= -52,
                lat >= 39,
                lat <= 69) 

#3. Split into local and utc time zone objects----
bird.local <- bird.x %>% 
  dplyr::filter(source!="eBird")
bird.utc <- bird.x %>% 
  dplyr::filter(source=="eBird")

#Visualize to check
hist(hour(bird.local$date))
hist(hour(bird.utc$date))

#C. CALCULATE OFFSETS####

#1. Get list of species for loop----
spp.tz <- expand.grid(tz=c("local", "utc"), species = sort(unique(bird.x$species)))

#2. Set up loop----
bird.o <- list()
for(i in 1:nrow(spp.tz)){
  
  flush.console()
  
  #3. Set tz method----
  tz <- spp.tz$tz[i]
  
  #4. Filter to species----
  if(tz=="local"){
    bird.i <- bird.local %>% 
      dplyr::filter(species==spp.tz$species[i])
  }
  if(tz=="utc"){
    bird.i <- bird.utc %>% 
      dplyr::filter(species==spp.tz$species[i])
  }

  #5. Make dataframe for prediction----
  x <- make_x(bird.i, tz)
  
  #6. Calculate offsets----
  o <- make_off(spp.tz$species[i], x)
  
  #7. Put together and save to list----
  bird.o[[i]] <- cbind(bird.i, o)
  
  print(paste0("Finished loop ", i, " of ", nrow(spp.tz), " - ", spp.tz$species[i]))
  
}

#D. SAVE####

#1. Collapse to dataframe----
bird <- rbindlist(bird.o)

#2. Save----
save(visit, bird, file=file.path(root, "01_NM4.1_data_offsets.R"))
