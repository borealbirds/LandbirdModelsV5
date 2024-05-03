# title: "'Evidence of absence' for removing over-prediction"
# (aka. projected range maps from raw survey data)
# author: Anna Drake
# date: Apr. 30, 2024

#1.Preamble -------

# Load libraries ---

library(doBy)
library(dplyr)
library(rasterVis)
library(tidyverse)
library(terra)
library(sf)

# set colour scale
col <- colorRampPalette(c("beige","lightgrey","darkgreen"))
                        
# NAD83(NSRS2007)/Conus Albers projection (epsg:5072) ----
crs<-"+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"

# Set root path for data on google drive ----
root<-"G:/Shared drives/BAM_NationalModels/NationalModels5.0" #PC link

# Get region shapefile ----
bcr<- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel.shp")) %>% 
  dplyr::filter(subUnit!=1) %>% # remove BCR 1
  st_transform(crs=crs) 

# Load data package with survey data ----
load(file.path(root, "Data", "04_NM5.0_data_stratify.R"))
rm("offsets", "cov", "bcrlist")
gc()

# Get any layer to produce a blank raster ----
Canada <- readxl::read_excel(file.path(root,"NationalModels_V5_VariableList.xlsx"), sheet = "ExtractionLookup")

rast<-area<-terra::rast(file.path(root, Canada$PredictPath[1],"ERAMAP_1km_1999.tif")) %>% mask(.,bcr)
rast[]<-0 # blank raster

#2. Calculate regional sampling intensity ----

Samp<-visit%>%.[11:12]%>% 
  st_as_sf(., coords = c("lon", "lat"),crs = 4326)%>% 
  st_transform(crs=crs) %>% 
  st_buffer(.,250000)

Samp<- rasterize(Samp,rast, fun="sum") # summarize surveys per 250km

#Opt: write out sampling intensity raster
#writeRaster(Samp,"Surveys_250km.tif")

#3. Get species occurrence data and determine probable occurrence ------------

Spp<-colnames(bird)

# set directory for output
setwd(file.path(root,"MosaicWeighting/ProbableRange/"))

# Run through species list  -------

for (i in c(19:length(Spp))){  #open species loop
  
out<-bird[,Spp[i]]
out<-as.data.frame(out)
out$id <- row.names(out)  
coord<-merge(visit[c(1,11:12)],out, by="id")

# extract sites where spp was seen and put 250 km buffer around sightings ------
Ypnts<-subset(coord,coord$out==1) 

Ypnts<-Ypnts%>%.[2:3]%>% 
  st_as_sf(., coords = c("lon", "lat"),crs = 4326)%>% 
  st_transform(crs=crs) %>% 
  st_buffer(.,250000)

Y <- rasterize(Ypnts,rast, fun="sum") # total sightings within 250km radius
Y[is.na(Y)]<-0
rm("Ypnts")
gc()

# Set a threshold sampling intensity at which we will accept coverage is sufficient to conclude a species is "not present"
# Values are set to exploit the "max" function in terra 

Thres<-Samp
Thres[Thres<1000]<-0.5 #less than 1000 surveys contributing
Thres[Thres>=1000]<-(-10)

#4. Calculate the proportion of sightings across study area and set a lenient "possibly present" threshold at 0.01% occurence -----

Weight<-Y/Samp
Weight[Weight<0.001]<-0  # 99.9% surveys do not have this spp: evidence suggests not present
Weight[Weight>=0.001]<-1 # 0.01% occurrence or greater: entertain the idea that the species is there

#5. Correct conclusions for low survey intensity regions -------
# In low survey intensity areas a detection likely represents a >0.01% occurrence rate
# However, in the absence of a detection, we don't have the power to say the sp. is *not* there as there hasn't been sufficient coverage
# Therefore, remove areas where conclusion is "not present" and coverage is <1000 surveys, but retain "present" 

Final<-c(Weight,Thres)
Final<- app(Final,fun="max") #0.5 (<1000 surveys sites) exceeds 0 but not 1; -10 (>1000 surveys sites) should not occur as 0/1 should overlay these areas
Final[Final==0.5]<-NA #set 0.5 to NA
Final<-mask(Final, bcr)

# Write out raster ------
writeRaster(Final, paste(Spp[i], "ProbRange.tiff"))

# Write out plot
Final[is.na(Final)]<-0.5
Final<-mask(Final, bcr)
png(filename=paste(Spp[i], "ProbRange.png"), width=5000, height=3000,res=300)
plot(Final, col=col(n=3), main=paste(Spp[i], "probable range"))
dev.off()
} #close species loop
