# ---
# title: National Models 5.0 - calculate offsets
# author: Elly Knight
# created: December 22, 2022
# ---

#NOTES################################

# This script uses the qpad-offsets package, which requires downloading that R project from github (https://github.com/borealbirds/qpad-offsets) to your local and setting the working directory to it (line 29).

# This script calculates offsets for ALL species available in QPAD, some of which may not fully satisfy the assumptions of QPAD. Further filtering of the species list should be conducted.

# This script uses the traditional QPAD method for offset calculation. Future versions of the models should consider the joint estimation approach for more accurate density estimation.

#This script uses the development version (dev-v4) of the QPAD package and qpad-offsets repo.

#wt_qpad_offsets() is now available as a function in wildRtrax and should replace much of the code below in future versions.

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
root <- "G:/Shared drives/BAM_NationalModels5"

#A. LOAD QPAD#####

#1. Set WD to qpad-offsets package----
setwd("C:/Users/elly/Documents/BAM/QPAD/qpad-offsets")

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
load(file.path(root, "01_NM5.0_data_clean.R"))

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
x.local <- make_x(visit.x.local, tz="local", check_xy=FALSE, tm=FALSE)
x.local$id <- visit.x.local$id
x.utc <- make_x(visit.x.utc, tz="utc",  check_xy=FALSE, tm=FALSE)
x.utc$id <- visit.x.utc$id
x <- rbind(x.local, x.utc)

#5. Add temporal covs to visits object
visit <- visit %>% 
  left_join(x %>% 
              mutate(tssr = 24*TSSR,
                     jday = JDAY*365) %>% 
              dplyr::select(id, tssr, jday))

#D. CALCULATE OFFSETS####

#1. Get list of species----
spp <- getBAMspecieslist()

#2. Set up output----
offsets <- data.frame(id=x$id)

#3. Make OFF----
for (i in 1:length(spp)) {
  cat(spp[i], "\n")
  flush.console()
  o <- make_off(spp[i], x, useMethod="n")
  offsets[,i+1] <- round(o$offset,4)
}
colnames(offsets) <- c("id", spp)

#4. Fix GRAJ####
offsets <- rename(offsets, CAJA = GRAJ)

#E. SAVE####

#1. Save----
save(visit, bird, offsets, file=file.path(root, "02_NM5.0_data_offsets.R"))
