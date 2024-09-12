###############################################
# title: "Defining maximal range extents from WildTrax, 
# augmented by 1980-2022 eBird observations"
# author: Anna Drake
# date: August 15, 2024
##################################################

# This code uses presence/absence from BAM model training data, augmented by any complete 
# ebird record up to 12hrs long and travelling up to 20km, to inform the probable maximal 
# range extents of BAM-modeled species. It is intended to be lenient, as the BAM model 
# will determine density based on suitable habitat within this range limit.
# Pixels are informed by point counts within a 250km radius. 
# We use two occurrence thresholds 0.2% and 0.1% of observations to say a species 
# could be present. E.g. at the 0.2% threshold, if:
# >99.8% of surveys *do not* have a sp: not present. Coded as 0.
# >=0.2% of surveys have a sp: we entertain the idea the species could be present. Coded as 1.
# Areas where there are <450 records + 'not present' are re-coded insufficient to determine absence (NA).
# Areas where there are <450 records + 'present' are left as 'present' 
# This is because:
# false positives at <450 records are rare (~2% if true occurrence rate=0.1% or 1/1000 (our "absent" threshold)); 
# false negatives at <450 records are substantial (~19% if true occurrence rate =1% or 1/100 (a rare "present" species))
# Finally, smoothed range extents are produced as shapefiles using the product rasters at both thresholds
# These shapefiles are simplified by: (1) joining regions within 100km of each other (2) filling holes
# equivalent to 100km diameter or less and, (3) smoothing contours using the 'ksmooth' function.

# Note:
# 1. originally we intended to use model-independent records but the inclusion WildTrax records greatly increases 
# areas in the north with sufficient coverage
# 2. These ranges may perform badly for species occupying patchy habitats that are poorly surveyed. For example
# species occupying high elevations (e.g. HOLA) will be absent from low elevation surveys that fall within 250km.
# If there is insufficient high elevation sampling they will be considered absent when the peaks should be coded as NA. 
# Here reducing the sampling radius and/or P/A threshold may help ID undersampled patches and correct range.

### PREAMBLE ######
# Load packages ------
library(devtools)
library(auk)
library(terra)
library(tidyverse)
library(dplyr)
library(sf)
library(sp)
library(raster)
require(smoothr)

`%notin%` <- Negate(`%in%`)

#1. Set colour scale ----------
col <- colorRampPalette(c("beige","lightgrey","darkgreen"))

#2. NAD83(NSRS2007)/Conus Albers projection (epsg:5072) ----
crs<-"+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"

#3. Set file paths ----
root <- "G:/Shared drives/BAM_NationalModels5"
output <-"C:/Users/andrake/Documents/eBird Output" #saved ebird output locally

#4. Get graphing layers -----
# Region shapefile ----
bcr<- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel.shp")) %>%
  dplyr::filter(subUnit!=1) %>% # remove BCR 1
  st_transform(crs=crs)

bcr_s<- st_simplify(bcr, preserveTopology = FALSE, dTolerance = 10000) #simplified BCR for plots

# Ocean layer ---
ocean<-read_sf("C:/Users/andrake/Downloads/ne_50m_ocean/ne_50m_ocean.shp")
ocean<- ocean %>% st_transform(crs(bcr)) %>% st_crop(bcr)

# Lakes layer ---
water<-read_sf("G:/Shared drives/BAM_NationalModels5/Regions/Lakes_and_Rivers/hydrography_p_lakes_v2.shp")
water<- water %>% st_transform(crs(bcr)) %>% st_crop(bcr) %>% dplyr::filter(TYPE==16|TYPE==18)

#5. Set up blank raster ---
rast<-rast(file.path(root,"PredictionRasters/AnnualClim/ERAMAP_1km_2019.tif"))
rast[]<-0

#6. Load data package with WildTrax data ----
load(file.path(root,"Data", "04_NM5.0_data_stratify.R"))
rm(offsets, bcrlist, birdlist,cov) # clean workspace
list<-colnames(bird)  # get all species 4-letter codes
WTrax<-visit %>% filter(source!="eBird") # extract non-ebird observations

### Set up a lat-lon-date index to check duplications between ebird and WildTrax records--------
ch<-WTrax[11:14]
ch$date<-substr(ch$date, 1,10)
ch$ID<-paste(ch$lat,ch$lon,ch$date, sep=".")

#7. Load 4-letter code to species name lookup ----
Spp<-read.csv(file.path(root,"data","Lookups","lu_species.csv"))%>%
  dplyr::select(species_code,species_name,species_common_name)
spp<-Spp %>% dplyr::filter(species_code %in% list) %>% .$species_common_name # get common name
spp<-spp[spp %notin% c("NONE", "Mew Gull", "Le Conte's Sparrow","Northwestern Crow","Pacific-slope Flycatcher", "Northern Goshawk")] # remove mismatch spp

#######################################################################################
# Part 1: Produce and eBird + WildTrax sampling intensity raster
########################################################################################
# Sampling file ----
sampling <- file.path(root,"/data/eBird/ebd_raw/ebd_sampling_relJul-2024/ebd_sampling_relJul-2024.txt")
# 'dummy' ebd file-----
ebd <- file.path(root,"data/eBird/ebd_raw/ebd_US-AK_relOct-2022/ebd_US-AK_relOct-2022.txt")

setwd(output)
AK_out_ebd <- "AKcleanEBD.txt"
Clean_sample <- "FullSampleClean.txt"

# Isolate modeling region using eBird "states" -----
CA_states<- dplyr::filter(ebird_states, country =="Canada") %>% .$state_code
US_states<-c("US-AK","US-ID","US-IA","US-CT","US-IL",
             "US-IN","US-MA","US-ME","US-MI","US-MN",
             "US-MT","US-NE","US-ND","US-NH","US-NY","US-OH",
             "US-OR","US-PA","US-RI","US-SD","US-VT","US-WA","US-WI")
statelist<-c(CA_states, US_states)

auk_ebd(ebd,sampling)

# Set up filter ------
ebd_filter<- auk_ebd(ebd, sampling) %>%
  auk_species(c(spp)) %>%
  auk_state(statelist)%>%
  auk_protocol(c("Traveling", "Stationary")) %>%
  auk_year(c(1980,2022))%>% # max/min of model training data
  auk_date(c("*-05-15", "*-07-15")) %>% #May 15-Jul 15 (to reduce non-breeding/migrant sightings)
  auk_duration(duration=c(0,720)) %>%   #<=12 hr long checklists
  #auk_distance(distance = c(0, 25))%>% #20 km or less - not using as this eliminates NAs and results in many dropped records
  auk_complete()  # only complete checklists

# Filter -----
auk_filter(ebd_filter, file = AK_out_ebd,
                           file_sampling = Clean_sample, overwrite=T)

# eBird sampling cleaned using filter -----
Samp<-read_tsv(file.path(output,"FullSampleClean.txt"))
Samp<-subset(Samp,Samp$`EFFORT DISTANCE KM`<20|is.na(Samp$`EFFORT DISTANCE KM`)) # add traveling limit w/ NAs

# Merge data sources --------

# Remove eBird records that are in WildTrax ----
Samp$ID<-paste(Samp$LATITUDE,Samp$LONGITUDE,Samp$`OBSERVATION DATE`,sep=".")
Samp<-Samp%>%filter(ID %notin% ch$ID) #928641 obs

# Create a sampling ID list - ebd files seem to have more observations than in sampling record ----
SampleID<-c(Samp$ID,ch$ID) #1.7 million observations

# Reduce to sampling intensity raster with 250km buffer ----
WTrax<-WTrax%>%.[11:12]%>%
  st_as_sf(., coords = c("lon", "lat"),crs = 4326)%>%
  st_transform(crs=crs)

Samp<-Samp%>%.[15:16] %>%
  st_as_sf(., coords = c("LONGITUDE", "LATITUDE"),crs = 4326)%>%
  st_transform(crs=crs)

# Add a buffer of 250 km around each point (150 km dramatically reduces northern coverage) ----
All<-rbind(Samp,WTrax)%>% st_buffer(.,250000)

# Stack and sum buffers - 'rasterize' is memory intensive, so split up the file ----
All2<-list()
rows<-c(seq(1,nrow(All), by=round(nrow(All)/10,0)))
rows[11]<-nrow(All) # include the last few rows in group 10

for (b in c(1:10)){
All.b<-All[c(rows[b]:rows[b+1]),] # subset
All2[b]<- rasterize(All.b,rast, fun="sum") # summarize surveys per 250km
}

# Put together write out sampling intensity raster -----
Total<-rast(All2)
Total <- app(Total, fun = "sum", na.rm = TRUE)
writeRaster(Total, paste0(output,"BackgroundSamplingIntensity_eBirdSup.tif", overwrite=T))

#################################################################################
# Part 2: Get species sighting data & determine presence/absence/not enough data
#################################################################################

#1. First eBird records ----------

# Get list of ebd objects to process----
ebd.files <- list.files(file.path(root,"data/eBird/ebd_raw"), pattern="ebd_*")
ebd.files<-ebd.files[-c(2)] # remove sampling file

# Set up species loop ----
for(i in c(1:length(ebd.files))){
  
# Set ebd path----
auk_set_ebd_path(file.path(root,"data/eBird/ebd_raw", ebd.files[i]), overwrite=TRUE)

# Filter data----

  ebd_filter<- auk_ebd(paste0(ebd.files[i],".txt")) %>%
    auk_species(c(spp)) %>%
    #auk_state(statelist)%>% #only relevant states are downloaded
    auk_protocol(c("Traveling", "Stationary")) %>%
    auk_year(c(1980,2022))%>% # max/min of model training data
    auk_date(c("*-05-15", "*-07-15")) %>% #May 15-Jul 15 (to reduce non-breeding/migrant sightings)
    auk_duration(duration=c(0,720)) %>%   #<=12 hr long checklists
    #auk_distance(distance = c(0, 25))%>% #25km or less - not using here as it eliminates NAs which results in many dropped records
    auk_complete()  # only complete checklists
  
#select columns to keep ----
    auk_filter(ebd_filter, file=file.path(output,paste0(ebd.files[i], "filtered.txt")), overwrite=TRUE,
                         keep = c("group identifier","sampling event identifier","observer id","scientific name", "common_name", "observation_count", "latitude", "longitude", "observation_date","effort_distance_km"))
}

#2. Check list of processed files----
ebd.files.done <- list.files(file.path(output), pattern="*filtered.txt")

#3. Extract eBird species detections -----
# Memory intensive: active removal and clearing needed to get R to drop used memory 

for (j in c(1:length(spp))){ # get species data ----  

Detections<-list()  
for (i in c(1:length(ebd.files.done))) { # get each region file

# Additional processing step ----
Samp<-read_tsv(file.path(output,paste0(ebd.files.done[i])), show_col_types = FALSE)
Samp<-subset(Samp,Samp$`EFFORT DISTANCE KM`<20|is.na(Samp$`EFFORT DISTANCE KM`)) # traveling limit applied to sampling layer
  
# Remove eBird records that are in WildTrax ----
Samp$ID<-paste(Samp$LATITUDE,Samp$LONGITUDE,Samp$`OBSERVATION DATE`,sep=".")
Samp<-Samp%>%filter(ID %notin% ch$ID)

# Remove eBird records that are not in eBird sampling document - unclear why a mismatch exists here ---
Samp<-Samp%>%filter(ID %in% SampleID)

Samp<-Samp%>%filter(`COMMON NAME`==spp[j])%>%  
    .[4:5] %>% st_as_sf(., coords = c("LONGITUDE", "LATITUDE"),crs = 4326)%>%
    st_transform(crs=crs)%>% st_buffer(.,250000) # get the lat-lon for all ids of the focal species

Detections[i]<- rasterize(Samp,rast, fun="sum") # summarize surveys per 250km
rm(Samp)
gc()
}

Species<-rast(Detections)%>%
  app(., fun = "sum", na.rm = TRUE)

rm(list=setdiff(ls(), c("i","j","Species","Samp","ch","spp",
                        "root","Detections","crs","rast",
                        "output","ebd.files.done","%notin%","SampleID")))

writeRaster(Species,paste0("Ebird_obs_",spp[j],".tiff"), overwrite=T)
print(paste("done",j,spp[j],sep=" "))
gc()
}

#4. Add in Wildtrax records ----

load(file.path(root,"Data", "04_NM5.0_data_stratify.R"))
rm(offsets, bcrlist, birdlist,cov) # clean workspace
WTrax<-visit %>% filter(source!="eBird") # extract non-ebird observations

# WildTrax sightings, removing eBird records ----
bird.df <- as.data.frame(as.matrix(bird))
bird.df$id <- as.numeric(row.names(bird.df))
bird.df<-bird.df%>%filter(id %in% WTrax$id) # only non-eBird records

rm("bird")

# Bring in sampling intensity raster ----
intensity<-rast("BackgroundSamplingIntensity_eBirdSup.tif")
intensity[is.na(intensity)]<-0

# extract sp ---
Spp<-read.csv(file.path(root,"data","Lookups","lu_species.csv"))%>%
  dplyr::select(species_code,species_name,species_common_name)

#7. Finally, produce extent maps -----

for (j in c(1:length(spp))){
  code<-Spp%>% filter(species_common_name==spp[j])%>%.$species_code
  obs<-cbind(bird.df$id,bird.df[code])%>%filter(.[code]>=1)
  colnames(obs)<-c("id","spp")

# if any records exist ----
if(length(obs>0)){
  loc<-visit%>% filter(id %in% obs$id) %>% .[11:12] %>%
    st_as_sf(., coords = c("lon", "lat"),crs = 4326) %>%
    st_transform(crs=crs)%>% st_buffer(.,250000) # get the lat-lon for all ids of the focal species

WSpp<- rasterize(loc,rast, fun="sum") 

# bring in eBird layer and add together
eBird<-rast(paste0("Ebird_obs_",spp[j],".tiff"))
r<-c(WSpp,eBird)
r <- app(r, fun = "sum", na.rm = TRUE)
r[is.na(r)]<-0

#8. Calculate the proportion of sightings across study area and set a lenient "possibly present" threshold 
# at 0.01% or 0.02% occurrence -----

Weight<-r/intensity  

if (minmax(Weight)[2]>1) {
  print(paste("error in ", spp[j], ": more counts than observations"))} # error check

else{
#Run for 0.1% and 0.2% ------
Weight[Weight<0.001]<-0 # >99.9% or >999/1000 surveys *do not* have this spp: evidence suggests not present
Weight[Weight>=0.001]<-1  # >=0.1% = 1/1000 or more surveys have occurrence: entertain the idea that the species is there

#9. Set a threshold sampling intensity at which we will accept coverage is sufficient to conclude a species is "not present"
# Values are set to exploit the "max" function in terra

Thres<-intensity
Thres[Thres<450]<-0.5 # label for where <450 surveys contribute
Thres[Thres>=450]<-(-10) # label for suitably sampled areas

#10. Correct conclusions for low survey intensity regions -------
# In low survey intensity areas we don't have the power to say the sp. is *not* there as there hasn't been sufficient coverage
# Therefore, remove areas where conclusion is "not present" and coverage is <450 surveys, but retain "present"

Final<-c(Weight,Thres)
Final<-mask(Final, bcr)
FinalTiff<-Final<- app(Final,fun="max") #0.5 (<450 surveys sites) exceeds 0 but is not 1; -10 (>=450 surveys sites) should not occur as 0/1 should overlay all such areas

FinalTiff[FinalTiff==0.5]<-NA #set 0.5 to NA for tiff file

#11. Outputs ---------
# Write out raster ------
writeRaster(FinalTiff, paste0(root,"/MosaicWeighting/Range/",code, "_ProbRange_01per.tiff"), overwrite=T)

# Write out plot -------
png(filename=paste0(root,"/MosaicWeighting/Range/",code, "_ProbRange_01per.png"), width=5000, height=3000,res=300)

plot(Final, col=col(3), main=paste(spp[j], "maximal range extent"), 
     legend=FALSE, axes=F, cex.main=3)
plot(ocean$geometry, lwd=0.1,col="#ccdce3", border = "#436175", add=T)
plot(water$geometry,lwd=0.05,col="white",border = "#436175", add=T)
plot(bcr_s$geometry,lwd=0.2, border="grey20",add=T)
legend("bottomleft", inset=(0.055), title="0.1% thres",
       c("Absent","NA","Present"), fill=col(3), cex=1.8)
dev.off()

gc()
print(paste0("done ", code," ", j))
} # end of 'else' statement
  }  # end of region
} # end of species

################################################################
# PART 3: convert "presence" to smoothed shapefile boundary 
################################################################

bcr_m<-st_union(bcr)%>%fill_holes(.,2000000) #region boundary - AK has some invalidities

for (j in c(1:length(spp))){ #run each species

code<-Spp%>% filter(species_common_name==spp[j])%>%.$species_code

# Import the different thresholds ---------
input<-list()
input[[1]]<-rast(paste0(root,"/MosaicWeighting/Range/Rasters/",code, "_ProbRange_01per.tiff"))
input[[2]]<-rast(paste0(root,"/MosaicWeighting/Range/Rasters/",code, "_ProbRange_2per.tiff"))

poly_done<-list()
for (k in c(1:2)){  
print(paste0("smoothing ", k))
  
poly <-ifel(input[[k]]==1, input[[k]], NA) %>% as.polygons() %>% 
  st_as_sf()%>% fill_holes(., 7900000)%>% #fill in 100 km diameter holes
  st_buffer(50000) %>% st_union() %>% st_as_sf() %>% st_transform(crs=crs) #buffer and merge areas 100km apart

poly_smooth <- smooth(poly,"ksmooth",smoothness=200) %>% st_make_valid()
poly_done[[k]]<-st_intersection(poly_smooth,bcr_m)
}

Boundary <- do.call(rbind, poly_done)
Boundary$Limit<-c("0.1% limit","0.2% limit")

#save shapefile ----
st_write(Boundary,paste0(root,"/MosaicWeighting/Range/",code,"_rangelimit.shp"))
print(paste0("Done ", code))
gc()
}

##################### END OF CODE ##################################
