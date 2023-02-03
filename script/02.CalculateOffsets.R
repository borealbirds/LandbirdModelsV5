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
root <- "G:/Shared drives/BAM_NationalModels/NationalModels4.1/BirdData"

#A. LOAD QPAD#####

#1. Set WD to qpad-offsets package----
setwd("C:/Users/Elly Knight/Documents/BAM/Projects/QPAD/qpad-offsets")

#2. Load QPAD requirements----
#2a. Estimates----
load_BAM_QPAD(version = 4)

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

#2. Format visit data for offset calculation---
visit.x <- visit %>% 
  mutate(time = str_sub(date, 12, 16),
         date = as.character(date(date))) %>% 
  rename(dur = duration, dis = distance, tagmeth = tagMethod)

#3. Split into local and utc time zone objects----
visit.x.local <- visit.x %>% 
  dplyr::filter(source!="eBird")
visit.x.utc <- visit.x %>% 
  dplyr::filter(source=="eBird")

#4. Format for offset calculation----
x.local <- make_x(visit.x.local, tz="local")
x.local$id <- visit.x.local$id
x.utc <- make_x(visit.x.utc, tz="utc")
x.utc$id <- visit.x.utc$id
x <- rbind(x.local, x.utc)

#C. CALCULATE OFFSETS####

#1. Get list of species----
spp <- getBAMspecieslist()

#2. Set up output----
offsets <- data.frame(id=x$id)

#3. Make OFF----
for (i in 1:length(spp)) {
  cat(spp[i], "\n")
  flush.console()
  o <- make_off(spp[i], x, useMethod="y")
  offsets[,i+1] <- round(o$offset,4)
}
colnames(offsets) <- c("id", spp)

#D. SAVE####

#1. Save----
save(visit, bird, offsets, file=file.path(root, "02_NM4.1_data_offsets.R"))
