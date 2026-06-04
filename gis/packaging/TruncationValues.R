###############################################
# title: "Truncation values"
# author: Elly Knight
# date: March 13, 2026
##################################################

#This code estimates two sets of species-specific truncation values for the final raster packaging:
#1. Upper truncation: max density we would expect based on the counts in the raw data. These are calculated from the raw data. All values above this are truncated to this value.
#2. Lower truncation: density values that account for 0.1 % of the population size (i.e., decimal dust). These are calculated from the These are converted to zero.

#PREAMBLE############################

#1. Load packages----
library(tidyverse) #basic data wrangling
library(QPAD) #offsets
library(terra)
library(sf)

#2. Set root path for data on google drive----
root <- "G:/Shared drives/BAM_NationalModels5"

#3. Load dataset ----
load(file.path(root, "data", "04_NM5.0_data_stratify.Rdata"))

#4. Load corrections ----
load(file.path(root, "data", "02_NM5.0_corrections.Rdata"))

#5. Get species list ----
spp <- colnames(birdlist |> dplyr::select(-bcr))

#GET UPPER TRUNCATION VALUES #########

#1. Some wrangling ----
bird.all <- as.data.frame(as.matrix(bird)) |> 
  rownames_to_column("id") |> 
  mutate(id = as.integer(id)) |> 
  pivot_longer(ALFL:YTVI, names_to="spp", values_to="count") 

corr.all <- corrections |> 
  pivot_longer(ALFL:YTVI, names_to="spp", values_to="correction") 

all <- inner_join(bird.all, corr.all) |> 
  left_join(visit |> 
              dplyr::select(id, location, year))

#2. Set up loop ----
q.out <- data.frame()
for(i in 1:length(spp)){
  
  #3. Get the species
  spp.i <- spp[i]
  
  all.i <- all |> 
    dplyr::filter(spp==spp.i)
  
  #4. Get species max count filter ----
  countmax <- stats::quantile(all.i$count, probs = 0.99, na.rm = TRUE)
  if(countmax==0){
    countmax <- stats::quantile(all.i$count, probs = 0.995, na.rm = TRUE)
  }
  if(countmax==0){
    countmax <- stats::quantile(all.i$count, probs = 0.999, na.rm = TRUE)
  }
  if(countmax==0){
    countmax <- stats::quantile(all.i$count, probs = 0.9999, na.rm = TRUE)
  }
  
  #5. Truncate counts & calculate density ----
  dat.i <- all.i |> 
    mutate(count = ifelse(count > countmax, countmax, count),
           density = count/correction)

  #6. Get species-specific truncation count ----
  #use varying qs depending on spp
  densmax.i <- stats::quantile(dat.i$density, probs = 0.99, na.rm = TRUE)
  thresh.i = 0.99
  if(densmax.i==0){
    densmax.i <- stats::quantile(dat.i$density, probs = 0.995, na.rm = TRUE)
    thresh.i = 0.995
  }
  if(densmax.i==0){
    densmax.i <- stats::quantile(dat.i$density, probs = 0.999, na.rm = TRUE)
    thresh.i = 0.999
  }
  if(densmax.i==0){
    densmax.i <- stats::quantile(dat.i$density, probs = 0.9995, na.rm = TRUE)
    thresh.i = 0.9995
  }
  if(densmax.i==0){
    densmax.i <- stats::quantile(dat.i$density, probs = 0.9999, na.rm = TRUE)
    thresh.i = 0.9999
  }
  
  #7. Output ----
  q.out <- rbind(data.frame(spp = spp.i,
                            thresh = thresh.i,
                            countmax = countmax,
                            densmax = densmax.i),
                 q.out)
  
  cat(i, " ")
  
}


#GET LOWER TRUNCATION VALUES #########

#1. Get list of mosaics ----
mosaiced <- data.frame(file = list.files(file.path(root, "output", "08_mosaics"), pattern="*.tif", recursive=TRUE)) |> 
  separate(file, into=c("region", "sppfolder", "spp", "year", "filetype"), remove=FALSE) |>  
  mutate(year = as.numeric(year)) |> 
  dplyr::select(-sppfolder, -filetype) |> 
  dplyr::filter(year==2020)

#2. Data limit mask ----
limit <- read_sf(file.path(root, "gis", "DataLimitationsMask.shp")) |> 
  st_transform("EPSG:3978") |> 
  vect()

#6. Get the water layer----
print("* Getting water layer *")
water_ca <- vect(read_sf(file.path(root, "gis", "WaterMask_Canada.shp")))
water_us <- vect(read_sf(file.path(root, "gis", "WaterMask_US.shp")))

#2. Set up dataframe ----
previous <- try(load(file.path(root, "data", "SpeciesPredictionTruncationValues.Rdata"))
)
if(inherits(previous, "try-error")){
  l.out <- data.frame()
  loop <- spp} else {
    loop <- spp[!spp %in% l.out$spp]
  }

#3. Take inventory ----


for(i in 1:length(loop)){
  
  spp.i <- loop[i]
  
  #3. Get truncation from above ----
  #Need this to get the pop total
  qsp <- q.out[q.out$spp==spp.i,]$q
  
  #4. Get the available mosaics ----
  mosaiced.i <- dplyr::filter(mosaiced, spp==spp.i)
  
  #3. Read in & sum the 2020 mosaic predictions----
  #we use the same methods as the packaging script before we take the final lower threshold
  if("Canada" %in% mosaiced.i$region){
    can.i <- rast(file.path(root, "output", "08_mosaics/Canada",
                                spp.i, paste0(spp.i, "_2020.tif"))) |> 
                   clamp(upper = qsp, values=TRUE) |> 
                   app(mean, na.rm=TRUE) |> 
      project("EPSG:3978")
    
    canq <- global(can.i, quantile, probs=0.999, na.rm=TRUE)[1,1]
    
    range.i <- rast(file.path(root, "gis", "ranges", paste0(spp.i, ".tif"))) |> 
      resample(can.i)
    
    canmask.i <- clamp(can.i, upper=q99, values=TRUE)*range.i |> 
      mask(limit) |> 
      mask(water, inverse=TRUE)
    
    canvals <- values(canmask.i, na.rm=TRUE)
  } else {
    canvals <- c(0)
    }

  if("Alaska" %in% mosaiced.i$region){
    ak.i <- rast(file.path(root, "output", "08_mosaics/Alaska",
                            spp.i, paste0(spp.i, "_2020.tif"))) |> 
      clamp(upper = qsp, values=TRUE) |> 
      app(mean, na.rm=TRUE) |> 
      project("EPSG:3978")
    
    akq <- global(ak.i, quantile, probs=0.999, na.rm=TRUE)[1,1]
    
    range.i <- rast(file.path(root, "gis", "ranges", paste0(spp.i, ".tif"))) |> 
      resample(ak.i)
    
    akmask.i <- clamp(ak.i, upper=q99, values=TRUE)*range.i |> 
      mask(limit) |> 
      mask(water, inverse=TRUE)
    
    akvals <- values(akmask.i, na.rm=TRUE)
    if(length(akvals)==0){akvals <- c(0)}
  } else {
    akvals <- c(0)
  }

  if("Lower48" %in% mosaiced.i$region){
    l48.i <- rast(file.path(root, "output", "08_mosaics/Lower48",
                            spp.i, paste0(spp.i, "_2020.tif"))) |> 
      clamp(upper = qsp, values=TRUE) |> 
      app(mean, na.rm=TRUE) |> 
      project("EPSG:3978")
    
    l48q <- global(l48.i, quantile, probs=0.999, na.rm=TRUE)[1,1]
    
    range.i <- rast(file.path(root, "gis", "ranges", paste0(spp.i, ".tif"))) |> 
      resample(l48.i)
    
    l48mask.i <- clamp(l48.i, upper=q99, values=TRUE)*range.i |> 
      mask(limit) |> 
      mask(water, inverse=TRUE)
    
    l48vals <- values(l48mask.i, na.rm=TRUE)
    if(length(l48vals)==0){l48vals <- c(0)}
  } else {
    l48vals <- c(0)
  }
  
  #6. Put the density values together, sort and calculate cumulative sum ----
  cumdensity <- data.frame(mean = canvals, region="Canada") |> 
    rbind(data.frame(mean = akvals, region="Alaska")) |> 
    rbind(data.frame(mean = l48vals, region="Lower48")) |> 
    rename(density = mean) |> 
    arrange(density) |> 
    mutate(cumulative = cumsum(density))
  
  #7. Caculate 
  #We use 0.5% because this is the CWS threshold for "rare but regular occurrence"
  popsum <- sum(cumdensity$density)
  popthresh <- popsum*0.0005
  densthresh <- dplyr::filter(cumdensity, cumulative < popthresh)
  
  #4. Get the output ----
  l.out <- rbind(data.frame(spp = spp.i) |> 
                   mutate(popsum = sum(cumdensity$density),
                          popthresh = popsum*0.0005,
                          denshthresh = max(densthresh$density)),
                 l.out) |> 
    arrange(spp) |> 
    unique()

  cat(i, " ")
  
  #7. Save ----
  save(q.out, l.out, file = file.path(root, "data", "SpeciesPredictionTruncationValues.Rdata"))
  
}


