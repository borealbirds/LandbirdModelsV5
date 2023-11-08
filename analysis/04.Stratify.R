# ---
# title: National Models 5.0 - stratify by BCR*country
# author: Elly Knight
# created: September 29, 2022
# ---

#NOTES################################

#In this script, we stratify & prepare the data for modelling. Steps include:

#1. BCR Attribution: Diving the data into BCRs for separate models. Each BCR is buffered by 100km. We do this so that we can feather predictions from adjacent regions together. The exception is the international boundaries between US and Canada, which we don't buffer because spatial layers for the covariates are different on either side of the border. In this case, we use a shapefiles of country boundaries to intersect with the buffered regions.

#2. Filtering the data. We remove data outside the study area, data before a minimum year, and data outside of dawn surveys. Additional dataset filtering steps can be added here as necessary.

#3. Thinning. First we assign each survey to a grid cell (1.9 km spacing) Then we randomly selecting one survey per cell per year for each bootstrap. Note on line 186-192, we use an alternative to sample_n(), which does provide a truly random sample and biases towards lines of data near the top of the dataframe. This step takes a while to run.

#4. Indicating species that should be modelled for each subunit based on a minimum detection threshold.

#Removal of non-dawn surveys is a secondary line of defense in case temporal patterns of QPAD offsets are unrealistic. In other words, we are restricting the dataset to times of day when p(availability) should be relatively high. This filter should be removed once QPAD offsets have been revisited. It should also be revisited if/when the national models include nocturnal species.

#PREAMBLE############################

#1. Load packages----

library(tidyverse) #basic data wrangling
library(terra) #basic raster handling
library(sf) #basic shapefile handling
library(fasterize) #fast rasterization of shapefiles
library(exactextractr) #fast & efficient raster extraction
library(dggridR) #grid for spatial thinning

#2. Set root path for data on google drive----
root <- "G:/Shared drives/BAM_NationalModels/NationalModels5.0"

#3. Load data packages with offsets and covariates----
load(file.path(root, "Data", "03_NM5.0_data_covariates.R"))

#4. Turn off scientific notation---
options(scipen=99999)

#BCR ATTRIBUTION#####################

#1. Make visit a vector object----
visit.v <- visit %>% 
  dplyr::select(id, lat, lon) %>% 
  st_as_sf(coords=c("lon", "lat"), crs=4326) %>% 
  st_transform(5072) %>% 
  vect()

#2. Read in country shapefiles----
can <- read_sf(file.path(root, "Regions", "CAN_adm", "CAN_adm0.shp")) %>% 
  st_transform(crs=5072) 
usa <- read_sf(file.path(root, "Regions", "USA_adm", "USA_adm0.shp")) %>% 
  st_transform(crs=5072) 

#3. Read in BCR shapefile----
bcr <- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel.shp")) %>% 
  mutate(area = as.numeric(st_area(read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel.shp")))),
         bcr=paste0("bcr", subUnit)) %>% 
  st_transform(crs=5072) %>% 
  dplyr::filter(subUnit != 0)

ggplot() +
  geom_sf(data=bcr, aes(fill=bcr))
  
#3. Set up loop for Canada BCR----
bcr.ca.list <- list(id=visit$id)
for(i in 1:nrow(bcr)){
  
  #4. Filter, buffer, & crop bcr shapefile to Canada----
  bcr.ca <- bcr %>% 
    dplyr::filter(row_number()==i) %>% 
    st_buffer(100000) %>% 
    mutate(bcr=1) %>% 
    st_intersection(can)
  
  if(nrow(bcr.ca) > 0){
    #5. Convert to raster for fast extraction----
    r <- rast(ext(bcr.ca), resolution=1000, crs=crs(bcr.ca))
    bcr.r <- rasterize(x=bcr.ca, y=r, field="bcr")
    
    #6. Extract raster value----
    bcr.ca.list[[i+1]] <- extract(x=bcr.r, y=visit.v)[,2]
    names(bcr.ca.list)[i+1] <- paste0("ca-", bcr$bcr[i])
  }
  
  print(paste0("Finished bcr ", i, " of ", nrow(bcr)))
  
}

#7. Set up loop for USA BCR----
bcr.usa.list <- list(id=visit$id)
for(i in 1:nrow(bcr)){
  
  #8. Filter, buffer, & crop bcr shapefile to Canada----
  bcr.usa <- bcr %>% 
    dplyr::filter(row_number()==i) %>% 
    st_buffer(100000) %>% 
    mutate(bcr=1) %>% 
    st_intersection(usa)
  
  if(nrow(bcr.usa) > 0){
    #9. Convert to raster for fast extraction----
    r <- rast(ext(bcr.usa), resolution=1000, crs=crs(bcr.usa))
    bcr.r <- rasterize(x=bcr.usa, y=r, field="bcr")
    
    #10. Extract raster value----
    bcr.usa.list[[i+1]] <- extract(x=bcr.r, y=visit.v)[,2]
    names(bcr.usa.list)[i+1] <- paste0("usa-", bcr$bcr[i])
  }

  print(paste0("Finished bcr ", i, " of ", nrow(bcr)))
  
}

#11. Collapse to dataframe and name----
bcr.df <- data.frame(do.call(cbind, bcr.ca.list)) %>% 
  left_join(data.frame(do.call(cbind, bcr.usa.list)))

#12. Remove columns (combos of BCR & country) that don't exist---
bcr.df <- bcr.df[colSums(!is.na(bcr.df))>0]

#FILTERING######################

#1. Set filters----

#Minimum year of data
minyear <- 1980

#Survey time
mintssr <- -1
maxtssr <- 6

#day of year
minday <- 135 #May 15
maxday <- 196 #July 15

#2. Remove points outside study area
bcr.df <- bcr.df[rowSums(!is.na(bcr.df[,c(2:ncol(bcr.df))]))>0,]

#3. Filter visits----
visit.bcr <- visit %>% 
  dplyr::filter(id %in% bcr.df$id,
                year >= minyear,
                tssr >= mintssr,
                tssr <= maxtssr,
                jday >= minday,
                jday <= maxday)

#14. Filter offset and bird objects by remaining visits----
offsets.bcr <- offsets %>% 
  dplyr::filter(id %in% visit.bcr$id)

bird.bcr <- bird %>% 
  dplyr::filter(id %in% visit.bcr$id)

#THIN FOR EACH BOOTSTRAP##############

#1. Get list of BCRs----
bcrs <- unique(colnames(bcr.df)[-1])

#2. Create grid for thinning----
grid <- dgconstruct(spacing = 2.5, metric=TRUE)

visit.grid <- visit.bcr %>% 
  mutate(cell = dgGEO_to_SEQNUM(grid, lon, lat)$seqnum)

length(unique(visit.grid$cell))

#3. Set number of bootstraps----
boot <- 100

#4. Set up loop----
bootstraps <- list()

for(i in 1:length(bcrs)){
  
  #5. Select visits within BCR----
  bcr.i <- bcr.df[,c("id", bcrs[i])] %>%
    data.table::setnames(c("id", "use")) %>%
    dplyr::filter(!is.na(use))
  
  #6. Filter visits----
  visit.i <- visit.grid %>% 
     dplyr::filter(id %in% bcr.i$id)
  
  #6. Set up bootstrap loop----
  out <- data.frame(id=bcr.df$id)
  for(j in 1:boot){
    
    #7. Set seed----
    set.seed(j)
    
    #8. Random sample----
    visit.j <- visit.i %>% 
      group_by(year, cell) %>% 
      mutate(rowid = row_number(),
             use = sample(1:max(rowid), 1)) %>% 
      ungroup() %>% 
      dplyr::filter(rowid==use)
    
    #9. Set up output----
    out[,(j+1)] <- bcr.df$id %in% visit.j$id
    
  }
  
  #10. Fix column names----
  x <- seq(1, boot, 1)
  colnames(out) <- c("id", paste0("X", x))
  bootstraps[[i]] <- out
  
  print(paste0("Finished thinning ", bcrs[i], ": ", i, " of ", length(bcrs)))
  
}

#11. Label by BCR----
names(bootstraps) <- bcrs

#UPDATE BIRD LIST BY BCR##############

#1. Set detection threshold for inclusion----
nmin <- 30

#2. Set up loop----
birdlist <- data.frame(id = bcr.df$id)

for(i in 1:length(bcrs)){
  
  #3. Select visits within BCR----
  bcr.i <- bcr.df[,c("id", bcrs[i])] %>%
    data.table::setnames(c("id", "use")) %>%
    dplyr::filter(!is.na(use))
  
  #4. Filter bird data----
  bird.i <- bird.bcr %>% 
    dplyr::filter(id %in% bcr.i$id)
  
  #5. Determine whether exceeds threshold----
  birdlist[i,c(2:nrow(bird.bcr))] <- colSums(bird.i) > nmin
  
}


#MAKE COVARIATE LOOKUP####

#SAVE#####

#1. Rename objects----
visit <- visit.bcr
bird <- bird.bcr
offsets <- offsets.bcr

save(visit, bird, offsets, bootstraps, birdlist, file=file.path(root, "04_NM5.0_data_stratify.Rdata"))
