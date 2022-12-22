# ---
# title: National Models 4.1 - calculate offsets
# author: Elly Knight, Peter Solymos
# created: December 22, 2022
# ---

#NOTES################################

# This script uses the qpad-offsets package, which requires downloading that R project from github (https://github.com/borealbirds/qpad-offsets) to your local and setting the working directory to it (line 22).

library(tidyverse)
library(qpad)
library(QPAD)
library(maptools)
library(intrval)
library(raster)

#1. Load data package from script 01----
load("01_NM4.1_data_clean.RDS")

#2. Set WD to qpad-offsets package----
setwd("C:/Users/Elly Knight/Documents/BAM/Projects/QPAD/qpad-offsets")

#3. Load QPAD requirements----
#3a. Estimates----
load_BAM_QPAD(version = 4)

#3b. Raster data----
rlcc <- raster("./data/lcc.tif")
rtree <- raster("./data/tree.tif")
rtz <- raster("./data/utcoffset.tif")
rd1 <- raster("./data/seedgrow.tif")
crs <- proj4string(rtree)

#3c. Source functions----
source("functions.R")



## species of interest
spp <- "OVEN"

## date and time
## https://en.wikipedia.org/wiki/ISO_8601
dt <- "2019-06-07" # ISO 8601 in YYYY-MM-DD (0-padded)
tm <- "05:20" # ISO 8601 in hh:mm (24 hr clock, 0-padded)

## spatial coordinates
lon <- -113.4938 # longitude WGS84 (EPSG: 4326)
lat <- 53.5461 # latitude WGS84 (EPSG: 4326)

## point count duration 
## and truncation distance (Inf for unlimited)
dur <- 10 # minutes
dis <- 100 # meters




x <- make_x(dt, tm, lon, lat, dur, dis)
str(x)
##'data.frame':	1 obs. of  8 variables:
## $ TSSR  : num 0.0089
## $ JDAY  : num 0.43
## $ DSLS  : num 0.14
## $ LCC2  : Factor w/ 2 levels "Forest","OpenWet": 2
## $ LCC4  : Factor w/ 4 levels "DecidMixed","Conif",..: 3
## $ TREE  : num 2.55
## $ MAXDUR: num 10
## $ MAXDIS: num 1



o <- make_off(spp, x)
str(o)
##'data.frame':	1 obs. of  5 variables:
## $ p         : num 0.991
## $ q         : num 0.562
## $ A         : num 3.14
## $ correction: num 1.75
## $ offset    : num 0.559


SPP <- getBAMspecieslist()
OFF <- matrix(0, nrow(x), length(SPP))
rownames(OFF) <- rownames(x) # your survey IDs here
colnames(OFF) <- SPP

for (spp in SPP) {
  cat(spp, "\n")
  flush.console()
  o <- make_off(spp, x)
  OFF[,spp] <- o$offset
}
str(OFF)
##num [1, 1:151] 0.1365 0.9699 0.0643 0.5917 -0.3132 ...
## - attr(*, "dimnames")=List of 2
##  ..$ : chr "17"
##  ..$ : chr [1:151] "ALFL" "AMCR" "AMGO" "AMPI" ...
