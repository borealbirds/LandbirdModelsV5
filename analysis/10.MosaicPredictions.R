#-----------------------------------------------------
# Title: Mosaicking prediction rasters by bootstrap with extrapolation mask & border blend
# Author: Anna Drake
# Date: March 2024           
#-----------------------------------------------------------------------------

# Extrapolation analysis is run for each species-BCR-bootstrap. Because of time-varying
# prediction rasters (e.g. biomass variables), year needs to be specified as well.

#------------------------------------------------
# This code does the following:
# 1) Identify areas of a BCR where the BCR-, species- & year-specific suite of predictor variables
# will cause model extrapolation (Type I and Type II) 
# 2) Identify which of these areas have alternate prediction layers that are 
# *not* extrapolated (originating from the buffer predictions of neighbouring BCRs)
# 3) Where a better option exists, clip out the extrapolated areas at the BCR level 
# 4) Where no better option exists retain these areas - for non-overlapping areas this will be
# the extrapolation of 1 model, for overlap areas where all models have no non-extrapolated
# data, it will be the weighted average output
# 5) Mosaic and output species-and-year-specific region-wide "Extrapolated area" raster that 
# shows areas where the only option was to retain extrapolated predictions
# 6) Mosaic and output the species-specific region-wide predictions using distance weighting
# to blend BCR predictions and using the clipping from (3) to retain only non-extrapolated
# predictions in overlap zones that have multiple outputs

# --------------------
# This code requires:
# 1) All functions in dsmextra_fn.R (or downloading the "dsmextra" package)
# 2) USA/CAN BCR overlap rasters and individual BCR weighting rasters produced 
# in "09.WeightingRasters.R"
# -------------------

# Load packages -------- 
library(sf)
library(tidyverse)
library(terra)
library(leaflet)
`%notin%` <- Negate(`%in%`)

#1. NAD83(NSRS2007)/Conus Albers projection (epsg:5072) ----
crs<-"+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"

#2. Set file paths ----
root<-"G:/Shared drives/BAM_NationalModels/NationalModels5.0" #PC link
setwd(file.path(root,"MosaicWeighting"))

#3. Canada and USA BCR boundaries ----

usa <- read_sf(file.path(root, "Regions", "USA_adm", "USA_adm0.shp")) %>% 
  st_transform(crs=crs) # US boundary

can <- read_sf(file.path(root, "Regions", "CAN_adm", "CAN_adm0.shp")) %>% 
  st_transform(crs=crs) # Canadian boundary

bcr.us <- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel.shp")) %>% 
  dplyr::filter(subUnit!=1) %>% # remove BCR 1
  st_transform(crs=crs) %>%
  st_intersection(usa) %>% 
  mutate(country="us") # Restrict to USA BCRs

bcr.ca <- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel.shp")) %>% 
  dplyr::filter(subUnit!=1) %>% # remove BCR 1
  st_transform(crs=crs) %>%
  st_intersection(can) %>% 
  mutate(country="ca") # Restrict to Canadian BCRs

#4.Functions for extrapolation analysis: pull appropriate rasters ---------------

#(1) Function to extract raster list 
# Returns matrix of raster pathways [1,] + variable names [2,]

Get_paths<- function(Variables, YR){
  
  Path_list<-Name_list<-c()
  Group.j<-dplyr::filter(Lookup,Lookup$Label %in% Variables)
  
  for (k in (1:nrow(Group.j))){ # Get variables within each group
    
    Period<-Static<-Name<-NULL
    
    # Return the closest year to desired period
    if(Group.j$TemporalResolution[k]=="match") {
      
      years<-dir(file.path(root,Group.j$PredictPath[k]))%>%
        stringr::str_sub(.,-8,-5)%>%as.numeric(.)%>%suppressWarnings()  #pull all years
      
      Match<-which(abs(years - YR) == min(abs(years - YR),na.rm=T))%>%
        years[.]%>%unique(.)  #ID closest year (minimum value when two options)
      
      if(Group.j$YearMatch[k]==0){
        
        Period<-paste(Group.j$PredictPath[k],"/",Group.j$Raster[k],"_",Match,".tif",sep="")
      }    
      
      else if (Group.j$YearMatch[k]==(-1)){
        
        Lag<-min(which(abs(years[years!=Match]- (YR-1)) == min(abs(years[years!=Match] - (YR-1)),na.rm=T)))%>%
          years[.]%>%unique(.)
        
        Period<-paste(Group.j$PredictPath[k],"/",Group.j$Raster[k],"_",Lag,".tif",sep="")
      }
      
    } else {
      Static<-paste(Group.j$PredictPath[k],"/",Group.j$Raster[k],".tif",sep="")
    }
    
    Path_list<-c(Path_list,Period,Static) #early and static list
    Name_list<-c(Name_list,Group.j$Label[k]) 
    
  } # close group
  
  return(list(Path_list,Name_list))
}

#(2) Function to stack rasters, crop by BCR, and convert to data frame 
# Uses "Get_paths" function to extract rasters for stacking

create_df <- function(Variables,YR,bcr.i){
  
  Stack<-list()  
  e<-ext(bcr.i)
  
  call<-BCR_lookup%>%dplyr::filter(subUnit==bcr.i$subUnit)%>%.$Names #match to bcr name
  covariates<-covlist%>%dplyr::filter(bcr==call[1])%>%.[,.==TRUE]%>% colnames()#get co-variates, assumes that BCR covariates don't differ from pooled BCRs (81,82 and 41,42)
  
  Var<-Get_paths(Variables,YR=YR) #get raster paths
  
  Names<-Var[2]%>%unlist(.)
  Retain<-Names %in% covariates  # ID ones we want by names [2,]
  
  Path<-data.frame(unlist(Var[1]),Retain)%>%dplyr::filter(Retain==TRUE)%>%.[,1]
  Label<-data.frame(Names,Retain)%>%dplyr::filter(Retain==TRUE)%>%.[,1]
  
  for (k in 1:length(Path)) { #stack rasters
    Stack[[k]]<-terra::rast(file.path(root, Path[k])) %>% 
      project(.,crs)  %>%
      crop(.,e) %>%
      mask(.,bcr.i) 
  }
  
  DF<-rast(Stack[1:length(Stack)]) %>% #turn into dataframe
    as.data.frame(., xy=TRUE) 
  
  names(DF)<-c("x","y",Label)
  
  return(DF)
}

#5. Load lookup tables ----------

#Variable look-up table:
Lookup <- readxl::read_excel(file.path(root,"NationalModels_V5_VariableList.xlsx"), 
                             sheet = "ExtractionLookup") 

Lookup$Raster<-Lookup$Label
Lookup$Raster[Lookup$Raster=="ERAPPTsmt_1km"]<-"ERAPPTsm_1km" # remove the "t" label for pulling raster (t vs t-1)
Lookup$Raster[Lookup$Raster=="ERATavesmt_1km"]<-"ERATavesm_1km" # remove the "t" label for pulling raster (t vs t-1)

#Load data package with covariates and BCR lists
load(file.path(root, "Data", "04_NM5.0_data_stratify.R"))
rm(bird,birdlist,gridlist,offsets)
gc()

#Limit covariates to those with continuous values
cov_clean<-cov%>%dplyr::select(where(is.numeric))
cov_clean<-cov_clean[names(cov_clean)!="hli3cl_1km"] #remove hli - it is categorical

#Produce look-up table: "BCR label"->"spatial sf subunit": **Probably want to just make this an object in stratify
BCR_Names <- colnames(bcrlist[-1])
BCR_lookup <- data.frame(matrix(nrow=length(BCR_Names),ncol=0))
BCR_lookup$Names<-BCR_Names

BCR_lookup$subUnit<-as.numeric(gsub("\\D", "", BCR_Names))
BCR_lookup$Region<-gsub("[^a-zA-Z]", "", BCR_Names) # can vs. usa

First<-Second<-Third<-subset(BCR_lookup,BCR_lookup$subUnit>100) #set multi. component BCRs to link to merged BCRs
First$subUnit<-as.numeric(substr(First$subUnit,1,2))
Second$subUnit<-as.numeric(substr(Second$subUnit,3,4))
Third$subUnit<-as.numeric(substr(Second$subUnit,5,6))
BCR_lookup<-rbind(subset(BCR_lookup,BCR_lookup$subUnit<100),
                  First,Second,subset(Third,!is.na(Third$subUnit))) 
rm(First,Second,Third,BCR_Names) #cleanup

#6. Extrapolation Analysis ----

### OK select species and year and Region ("can" vs "usa")
SPP<-"OSFL"
YR <-1985
Region<-"can"

if (Region=="can"){
  BCRsource="bcr.ca"
  BCR=bcr.ca 
  BORDER=can
}

if (Region=="usa"){
  BCRsource="bcr.us"
  BCR=bcr.us
  BORDER=usa
}

# Run through each BCR...---------------------------------------------

for(j in c(1:nrow(BCR_lookup))){  # select BCR  
  
  # Get files for region
  if(BCR_lookup$Region[j]!=Region){ 
    cat("Skipping", BCR_lookup$Region[j], "BCRs...")
    next
  }
  
  file<-list.files(file.path(root, "Output","bootstraps"),pattern=paste(SPP,"_",BCR_lookup$Names[j],".+",sep=""))
  
  # Confirm all bootstraps available
  if(length(file)!=10){ 
    cat("missing bootstraps for",BCR_lookup$Names[j], ":skipping...\n")
    next
  }
  
  # Check if species-BCR-year has already run in totality
  ck<-list.files(file.path(root,"MosaicWeighting","MaskingRasters"),
                pattern=paste("ExtrapolatedArea_",BCR_lookup$Names[j],".*",SPP,"_",YR,".tif", sep=""))
  if(length(ck)==10){
    cat("extrapolation complete, skipping...\n")
    next
  }
  
  #Note: covariates are the same across boots, retrieve from one...
  
  boots=load(file.path(root, "Output","bootstraps",file[1]))
  gbm = get(boots[1])
  
  Variables<-gbm$var.names %>% 
    subset(.,. %in% names(cov_clean)) # get continuous covariates
  
  Variables<-subset(Variables,Variables %notin% "TRI_1km") #temporary fix for TRI
  
  rm(boots,gbm, out.i) # clean up

#Buffer BCR  -----------------------
  bcr.i <-  BCR %>% 
    dplyr::filter(subUnit==BCR_lookup$subUnit[j]) %>% 
    st_buffer(100000)%>%  # Filter & buffer shapefile
    st_intersection(., BORDER) #crop at border
  
# Extract BCR specific data frames for each group. Year=2020  ---
  cat("Retrieving",SPP,BCR_lookup$Names[j], "prediction values for", YR, "...")  
  Dataframes<-create_df(Variables, YR=YR, bcr.i=bcr.i)
  
# Run through the bootstraps -------------------------
  
  Raster<-list()
  cat("Running bootstraps \n")
  
  for (k in c(1:10)){  #run through each bootstrap
    
    boots=load(file.path(root, "Output","bootstraps",file[k]))
    ids=visit.i$id
    
    # Get sample data for BCR -------------------------
    sample.i<-cov_clean%>%subset(.,.$id %in% ids) 
    
    # Calculate extrapolation --
    target=Dataframes
    
    # Ensure match training-predictor variables ------------------------
    sample<-sample.i%>% dplyr::select(one_of(names(target))) #drop any in ref not present in target - will get an xy warning
    target<-target%>%dplyr::select(c("x","y",names(sample))) #drop any in target not present in ref
    
    xp = names(sample)
    
    Extrapol <- compute_extrapolation(samples = sample,
                                      covariate.names = xp,
                                      prediction.grid = target,
                                      coordinate.system = crs,
                                      verbose=F)

#7.Output #######
    
    #Produce binary raster of extrapolation area
    
    Raster[[k]]<-as(Extrapol$rasters$mic$all, "SpatRaster")
    Raster[[k]][Raster[[k]]>0] <- 1  # extrapolated locations in each boot
    Raster[[k]][Raster[[k]]!=1|is.na(Raster[[k]])]<-0 #eliminate NA areas
    Raster[[k]]<-mask(Raster[[k]],bcr.i) # crop to BCR + buffer
    
    writeRaster(Raster[[k]],file.path(root,
                "MosaicWeighting","MaskingRasters",
                paste("ExtrapolatedArea_",BCR_lookup$Names[j],"_boot",k,"_",SPP,"_",YR,".tif", sep="")), 
                overwrite=T)
    
    # Print progress
    cat("Finished", k, "of", 10, "bootstraps \n")
    
  } #end of bootstraps - all extrapolated areas
  
  #if (verbose) {
  cat("Processed", j, "of",nrow(BCR_lookup), "BCRs \n")
  #}
} # end of BCRs


#8. Mosaic together predictions using edge weighting & removing extrapolation areas where alternate predictions exist

#Determine region

if(Region=="can") {
  MosaicOverlap<-terra::rast(file.path(root,"MosaicWeighting","ModelOverlap_Can.tif"))
}

if(Region=="usa") {
  MosaicOverlap<-terra::rast(file.path(root,"MosaicWeighting","ModelOverlap_US.tif"))
}

#if (verbose) {
cat("Mosaicking", Region, YR, "predictions for", SPP,"\n")
# }

# Import BCR extrapolation rasters by bootstrap -------

for (k in c(1:10)){ # open bootstraps
  
f <- lapply(list.files(path=file.path(root,"MosaicWeighting","MaskingRasters"), 
            pattern=paste("_boot",k,"_",SPP,"_",YR,"*", sep=""), full.names = TRUE),
            terra::rast)  
f <- sprc(f)
f.mosaic <- mosaic(f, fun="sum") # sum of overlapping extrapolation 
f.mosaic<-crop(f.mosaic, ext(MosaicOverlap))
MosaicOverlap<-crop(MosaicOverlap, ext(f.mosaic)) #deal with any extent mismatch

# Divide extrapolation layers by available prediction layers to determine if alternate exists -----
Overlay<-f.mosaic/MosaicOverlap # If ==1 there is no suitable raster, if ==0.25-0.75 then a suitable raster exists
plot(Overlay)

# Produce region-wide extrapolation and correction raster ----------

Extrapolation<-Correction<-Overlay
Extrapolation[Extrapolation<1]<-0  #Areas we retain extrapolation because we have no alternative == 1
writeRaster(Extrapolation,paste("ExtrapolatedArea_","boot",k,"_",Region,"_",SPP,"_",YR,".tif",sep=""), overwrite=T) #Save 

Correction[Correction==1]<-0 #ignore areas where nothing *or* everything is missing 
Correction[Correction>0]<- 1 #flag areas where non-extrapolated predictions exist in any layer

# Correct extrapolation rasters at BCR level ----------

MosaicStack<-list() # hold the weighting rasters
SppStack<-list() #hold the species prediction rasters

for(j in c(1:nrow(BCR_lookup))){  # select BCR  
  
  if(BCR_lookup$Region[j]!=Region) {next} # skip over us/can depending on region
  
  c<-list.files(file.path(root, "MosaicWeighting","MaskingRasters"),
                paste("ExtrapolatedArea_",BCR_lookup$Names[j],"_boot",k,"_",SPP,"_",YR,".tif", sep="")) 
  
  if(length(c)==0){ 
    cat("Missing",BCR_lookup$Names[j], "extrapolation output for", YR, SPP, "boot",k,"\n")
    next}
  
  c<-terra::rast(file.path(root, "MosaicWeighting","MaskingRasters",c)) %>%crop(.,Correction) #correct edges (Newfoundland is off)
  
  cor<-Correction%>%crop(.,c)%>%mask(.,c)
  cor<-cor*c #if missing and substitute exists (1*1=1), if not missing or no substitute (0*1/1*0=0)
  cori<-classify(cor,cbind(1,0), others=1) #invert values to give regions with substitute values no weight
  
  # Modify weighting raster to account for removed sections ------------
  
  w<-terra::rast(list.files(path=file.path(root,"MosaicWeighting","CAN_BCR_Weighting"),
                            pattern=paste("BCR_",BCR_lookup$subUnit[j], sep=""), full.names = TRUE))
  
  w<-resample(w, cori) #different origin, re-sample the weight raster
  w<-w*cori  # 0-out areas of extrapolated values where alternate predictions exist
  MosaicStack[[j]]<-w # stack to produce raster by which we will divide output

  # Apply BCR weighting raster to Prediction raster 
  plist<-list.files(path=file.path(root,"Output","predictions"),
                    pattern=paste(SPP,"_",BCR_lookup$Names[j],"_",k,"_",YR,".tiff", sep=""), full.names = TRUE) # 
  if (length(plist)==0) {cat("Missing bootstrap", k, "predictions for", SPP, BCR_lookup$Names[j], "in", YR,"\n")
    next} 
  
  p<-terra::rast(file.path(plist))
  p<-crop(p,w) #ensure match
  pw<-p*w  #apply weighting to the predictions
  SppStack[[j]]<-pw
} # end of BCRs

# !sapply... deals with missing layers, can remove ultimately
CANwideW<-MosaicStack[!sapply(MosaicStack,is.null)]%>%sprc(.)%>%mosaic(., fun="sum") #sum weighting (divisor)
CANwideSp<-SppStack[!sapply(SppStack,is.null)]%>%sprc(.)%>%mosaic(., fun="sum") # sum weighted predictions

#the following is just to deal with missing layers, once run is complete it shouldn't need
CANwideSp<-crop(CANwideSp,CANwideW)
CANwideW<-crop(CANwideW,CANwideSp)

# correct the weighting -------------
FinalOut<-CANwideSp/CANwideW #== weighted average
FinalOut<-mask(FinalOut,bcr.ca)
writeRaster(FinalOut,paste("PredictionMosaic_boot",k,"_",Region,"_",SPP,"_",YR,".tif",sep=""), overwrite=T)
plot(FinalOut)
} # end of bootstraps

########## END OF CODE ########################
