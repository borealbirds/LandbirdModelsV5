# ---
# title: National Models 5.0 - Gap analysis of all covariates
# author: Anna Drake
# created: June 30, 2024
# ---

#NOTES################################
# This code uses dsextra to quantify extrapolation
# This package calls the "ExDet" function which needs it's
# tolerances adjusted to avoid singular matrices (rounding low values down)

#PREAMBLE############################

#1. Load packages----

# Import dsmextra
#if (!require("remotes")) install.packages("remotes")
#remotes::install_github("densitymodelling/dsmextra")
#-----------------------------
#Tweak ExDet function in dsmextra by increasing 
#the tolerance of the Mahalanobis function: 
#getAnywhere(ExDet)

pacman::p_load(dsmextra, sp, tidyverse, sf, ggplot2, 
               rlang, leaflet, dplyr, terra, smoothr, stringi,
               readr, ecospat, mapview)

# extra function -----
`%notin%` <- Negate(`%in%`)

#2. NAD83(NSRS2007)/Conus Albers projection (epsg:5072) ----
crs <- "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"

#3. Set file paths ----
root<-"G:/Shared drives/BAM_NationalModels5"

#4. Load data package with covariates and BCR lists ----
load(file.path(root, "Data", "04_NM5.0_data_stratify.R"))
rm(bird,birdlist,gridlist,offsets)

cov_clean<-cov%>%dplyr::select(where(is.numeric)) #limit to numeric covariates

#5. Import BCR boundaries with political boundaries + the political boundaries alone ----
bcr <- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel_Unbuffered.shp")) %>% 
  st_transform(crs=crs)

can <- read_sf(file.path(root, "Regions", "CAN_adm", "CAN_adm0.shp")) %>% 
  st_transform(crs=crs) # Canadian boundary

us <- read_sf(file.path(root, "Regions", "USA_adm", "USA_adm0.shp")) %>% 
  st_transform(crs=crs) # US boundary

#6. Set up a bunch of formatting and bring in water layers for plots ------

# Ocean layer ---
ocean<-read_sf("C:/Users/andrake/Downloads/ne_50m_ocean/ne_50m_ocean.shp")
ocean<- ocean %>% st_transform(crs(bcr)) %>% st_crop(bcr)

# Lakes layer ---
water<-read_sf("G:/Shared drives/BAM_NationalModels5/Regions/Lakes_and_Rivers/hydrography_p_lakes_v2.shp")
water<- water %>% st_transform(crs(bcr)) %>% st_crop(bcr) %>% dplyr::filter(TYPE==16|TYPE==18)

# Create labeling system for categories (using codes introduced below) -------

label<-c("No extrapolation","Climate (Annual)", #0,1
         "Vegetation","Vegetation + Climate (Annual)", #10,11
         "Climate (Normal)","Climate (Annual + Normal)", #100,101
         "Vegetation + Climate (Normal)", "Vegetation + Climate (Annual + Normal)",#110,111
         "Topography","Topography + Climate (Annual)", #1000,1001
         "Topography + Vegetation","Topography + Climate (Annual) + Vegetation",#1010,1011
         "Topography + Climate (Normal)","Topography + Climate (Annual + Normal)",#1100,1101
         "Topography + Climate (Normal) + Vegetation", "Topography + Climate (Annual + Normal) + Vegetation", #1110,1111
         "Wetland", "Wetland + Climate (Annual)", #10000,10001
         "Wetland + Vegetation",  "Wetland + Climate (Annual) + Vegetation", #10010,10011
         "Wetland + Climate (Normal)","Wetland + Climate (Annual + Normal)", #10100, 10101
         "Wetland + Climate (Normal) + Vegetation","Wetland + Climate (Annual + Normal) + Vegetation", #10110,10111,
         "Wetland + Topography", "Wetland + Topography + Climate (Annual)",#11000,11001
         "Wetland + Topography + Vegetation","Wetland + Topography + Climate (Annual) + Vegetation",#11010,11011,
         "Wetland + Topography + Climate (Normal)", "Wetland + Topography + Climate (Annual + Normal)", # 11100,11101,
         "Wetland + Topography + Climate (Normal) + Vegetation",
         "Wetland + Topography + Climate (Annual + Normal) + Vegetation")#11110,11111

cls <- data.frame(id=c(0,1,10,11,100,101,110,111,
                       1000,1001, 1010,1011, 1100,1101, 1110,1111,
                       10000,10001, 10010,10011, 10100,10101, 10110,10111,
                       11000,11001, 11010,11011, 11100,11101,11110,11111), cover=label)

# Create contrasting colour scheme for easier visualization ----------

c32<-c("lightgrey","gold1","limegreen","chartreuse1","darkgoldenrod2", #climate normal
          "darkgoldenrod4","darkseagreen","darkgreen","darkred","darkorchid1",#topography+climate(annual)
          "salmon2","pink","purple2","purple4","hotpink",#topography + normal + vegetation
          "deeppink3","dodgerblue4","dodgerblue","yellow","tan",#wetland+ annual + vegetation
          "blue","navy","grey30","black","turquoise1",#wetland + topography
          "yellow3","turquoise","magenta","yellow4","#633", # Wetland + Topography + Climate (Annual + Normal)
          "magenta3","magenta4")

#Improve order of variables in cls
levels <- label[c(1,2,5,6,3,4,7,8,9,11,10,13,14,12,15,16,17,19,
                  18,21,22,20,23,24,25,27,26,29,30,28,31,32)] #rearrange classes
levels(cls$cover)<-levels #apply

# Reorder cls
cls<-cls[c(1,2,5,6,3,4,7,8,9,11,10,13,14,12,15,16,17,19,
      18,21,22,20,23,24,25,27,26,29,30,28,31,32),] 

# Do the same to colour coding
c32<-c32[c(1,2,5,6,3,4,7,8,9,11,10,13,14,12,15,16,17,19,
           18,21,22,20,23,24,25,27,26,29,30,28,31,32)]

### Simplify BCR for plots-----------------
bcr_s<-st_simplify(bcr, preserveTopology = FALSE, dTolerance = 100)

#7. Load prediction raster categories -----
Lookup <- readxl::read_excel(file.path(root,"NationalModels_V5_VariableList.xlsx"), sheet = "ExtractionLookup")%>% 
  dplyr::select(Category,Label,YearMatch, TemporalResolution, PredictPath)

#Modify t-1 notation for ERA
Lookup$Raster<-Lookup$Label
Lookup$Raster[Lookup$Raster=="ERAPPTsmt_1km"]<-"ERAPPTsm_1km" # remove the "t" label for pulling raster (t vs t-1)
Lookup$Raster[Lookup$Raster=="ERATavesmt_1km"]<-"ERATavesm_1km" # remove the "t" label for pulling raster (t vs t-1)

#Remove hli3cl and SCANFI_1km - categorical, not continuous 
Lookup<-subset(Lookup, Lookup$Label %notin% c("hli3cl_1km","SCANFI_1km"))

#Get grouping categories and remove LCC
Lookup$Category<-ifelse(Lookup$Category=="Road","Disturbance",ifelse(Lookup$Category=="Greenup","Annual Climate", Lookup$Category)) ### need to finish aligning categories
Grouping<-unique(Lookup$Category) #Groupings
Grouping<- Grouping[! Grouping %in% c('Landcover','LCC_MODIS')] #LCC is NA for this

#MAIN ###################

#8. Run functions to retrieve and prep data for analysis -----

#(1) Function to extract raster list - returns matrix of raster pathways [1,] + variable names [2,]

Get_paths<- function(Grouping, YR){
  
  Path_list<-Name_list<-c()
  Group.k<-dplyr::filter(Lookup,Lookup$Category==Grouping)
 
  for (k in (1:nrow(Group.k))){ # Get variables within each group
  
    Matched<-Static<-Name<-NULL
    
    fdir<-dir(file.path(root,Group.k$PredictPath[k])) # full directory
    
    # Return the closest year to desired period
    if(Group.k$TemporalResolution[k]=="match") {
    
      years<-fdir[grep(Group.k$Raster[k], fdir)] %>%  #isolate variable of interest
      stringi::stri_sub(-8,-5)%>% as.numeric()%>% suppressWarnings()  #pull all years
    
    if(Group.k$YearMatch[k]==0){
    
      Match<- which(abs(years - YR) == min(abs(years - YR),na.rm=T)) %>% min(.) %>% #where equidistant years occur use the lower
      years[.]    #ID closest year 
      
      Matched<-paste(Group.k$PredictPath[k],"/",Group.k$Raster[k],"_",Match,".tif",sep="")
      }    
      
      else if (Group.k$YearMatch[k]==(-1)){
        
      Lag<-min(which(abs(years- (YR-1)) == min(abs(years- (YR-1)),na.rm=T))) %>% min(.) %>% # as above
      years[.]  # ID lagged closest year
        
      Matched<-paste(Group.k$PredictPath[k],"/",Group.k$Raster[k],"_",Lag,".tif",sep="")
      }
      
    } else {
      Static<-paste(Group.k$PredictPath[k],"/",Group.k$Raster[k],".tif",sep="")
    }
    
    Path_list<-c(Path_list,Matched,Static) #Join matched and static list
    Name_list<-c(Name_list,Group.k$Label[k]) 
    
  } # close group
  
  return(list(Path_list,Name_list))
}

#(2) Function to stack rasters, crop by BCR, and convert to data frame 
# this uses the above "Get_paths" function to extract appropriate rasters for stacking

create_df <- function(Grouping,YR,bcr.i){
  
  Stack<-list()  
  e<-ext(bcr.i)

  #ID BCR appropriate covariates --- 
  call<-paste(bcr.i$country,bcr.i$subUnit, sep="") #match to bcr name
  covariates<-covlist %>% dplyr::filter(bcr==call)%>%.[,.==TRUE]%>% colnames()#get co-variates
  
  #get Grouping raster paths ----
  Var<-sapply(Grouping,Get_paths,YR=YR)
  
  # filter by "covariates" list above
  Names<-Var[2,]%>%unlist(.)
  Retain<-Names %in% covariates  # ID ones we want by names [2,]
  
  Path<-data.frame(unlist(Var[1,]),Retain)%>%dplyr::filter(Retain==TRUE)%>%.[,1]
  Label<-data.frame(Names,Retain)%>%dplyr::filter(Retain==TRUE)%>%.[,1]
  
  for (k in 1:length(Path)) { #stack rasters
    Stack[[k]]<-terra::rast(file.path(root, Path[k])) %>% 
      project(.,crs)  %>%
      crop(.,e) %>%
      mask(.,bcr.i) 
  }
  
if (Grouping=="Disturbance" & i==33){
 Stack[[2]]<-crop(Stack[[2]],Stack[[4]])
 Stack[[3]]<-crop(Stack[[3]],Stack[[4]])} #patch for AK rasters extent mismatch
  
  DF<-rast(Stack[1:length(Stack)]) %>% #turn into dataframe
    as.data.frame(., xy=TRUE) 
  
  names(DF)<-c("x","y",Label)
 
  return(DF)
}

#9. Analysis ----
for (b in c(1985,1990,1995,2000,2005,2010,2015,2020)){ #open year 
BCR_list<-list() # for storing output by BCR

for(i in c(1:nrow(bcr))){  # select BCR
  
  #Buffer BCR  ---  
  bcr.i <-  bcr %>% 
    dplyr::filter(row_number()==i) %>% 
    st_buffer(100000)  # Filter & buffer shapefile
    
  if (bcr.i$country=="can"){
    bcr.i<-st_intersection(bcr.i, can) }#crop at Canadian border

  if (bcr.i$country=="usa"){
    bcr.i<-st_intersection(bcr.i, us) }#crop at US border
  
#Get BCR reference sample - pool samples for the pooled BCRs  --- 
call<-paste(bcr.i$country,bcr.i$subUnit, sep="") #match to bcr.i to bcr name in bcrlist
ids <- bcrlist %>% .[, c("id",call)] %>% subset(.,.[,2]==TRUE)%>%.$id
  
  if(length(call)==2) {
    id2<- bcrlist %>%.[, c("id",call)]%>%subset(.,.[,3]==TRUE)%>%.$id  
    ids<-unique(c(ids,id2))  }

# Get sample data for BCR ---
sample.i<-cov_clean%>%subset(.,.$id %in% ids) 
sample.i<-sample.i %>%.[,colSums(.,na.rm=T)!=0] # remove covariates with no data
  
#Convert sample locations to SPDF (for Leaflet plots)  ---  
#pnts<-visit%>% subset(.,.$id %in% ids)%>%.[11:12]%>% #sample lon, lat
#    st_as_sf(., coords = c("lon", "lat"),crs = 4326)%>% 
#    st_transform(crs=crs) %>%
#    as(.,"Spatial") 

# Extract BCR specific data frames for each group. pick year of interest e.g Year=2020  ---  

Dataframe<-purrr::map(Grouping,create_df, YR=b, bcr.i=bcr.i,.progress = T)  

#Calculate extrapolation ---

Ext<-Raster<-list()
  
 for (j in c(1:length(Grouping))){ #run through each group
  
target=Dataframe[[j]] #the region of interest
  
sample<-sample.i%>% dplyr::select(one_of(names(target))) #drop any in ref not present in target - will get an xy warning
target<-target%>%dplyr::select(c("x","y",names(sample))) #drop any in target not present in ref

xp = names(sample)

Ext[[j]] <- compute_extrapolation(samples = sample,
                                  covariate.names = xp,
                                  prediction.grid = target,
                                  coordinate.system = crs,
                                  verbose=F)
  
#10. Produce binary raster of where extrapolation areas are flagged for each grouping ----------
#NOTE values >1= combinatorial extrapolation; <0= uni-variate; 0-1 = analogue
  
Raster[[j]]<-Ext[[j]]$rasters$ExDet$all
#MIC[[j]]<-EarlyExt[[j]]$rasters$mic$all

Raster[[j]][Raster[[j]]>=0 & Raster[[j]]<=1] <-0  #set analogue to 0
Raster[[j]][Raster[[j]]<0]<-1 # univariate to 1
Raster[[j]][Raster[[j]]>1]<-1 # combinatorial to 1 # could code to 2 for other applications

cat("Finished Group", j, "\n")
} # end of groupings

BCR_list[[i]]<-Raster
cat("Finished", i, "of", 34, "\n")

} #end of BCRS

save(BCR_list,file = paste("Extrapolation", b,".RData",sep=""))  

#If starting here to rework plot appearance -------
#for(b in c(1985,1990,1995,2000,2005,2010,2015,2020)){
#load(paste("Extrapolation", b,".RData",sep=""))

#11. Mosaic BCRs into region-wide raster ------------------
# Reorganize, merge, and look at available data in overlap zones

#Import the areas of BCR overlap -----
M<-terra::rast(file.path(root,"gis","ModelOverlap.tif"))

# Set up lists ------
merge<-mergeb<-Final<-list()

# Reorder nesting: j in i rather than i in j ------

for (j in c(1:length(Grouping))) { 
for (i in c(1:34)){
  merge[[i]]<-rast(BCR_list[[i]][[j]]) # this is the original values
  mergeb[[i]]<-rast(BCR_list[[i]][[j]])
  mergeb[[i]][mergeb[[i]]==2]<-1 # this layer allows us to ID areas where data exists in overlapping BCR buffers
}
  
# Mosaic region for each group ----- 
f <- sprc(merge) %>% mosaic(.,fun="sum") %>% crop(.,ext(M))
M <- crop(M, ext(f)) # deal with any extent mismatch
f<-f/M #get rid of the additive areas

# Find areas where overlap contributes non-extrapolated data ----
b.spat <- sprc(mergeb)%>% mosaic(.,fun="sum") %>% crop(.,ext(M))
Overlay<-b.spat/M # Divide total extrapolation layers by available prediction layers to determine if alternate exists -----
Overlay[Overlay<1]<-0  #If all layers extrapolated == 1, otherwise == 0

Final[[j]]<-f*Overlay #mask out the areas where alt data exists
Final[[j]][Final[[j]]>2]<-NA #get rid of tiny mismatch in layers
}

#12. Categorize for stacking Groupings ------------------

Stacked<-Final
cat<-c(1,10,100,NA,1000,10000) #Disturbance has no extrapolation so drop 

for (i in c(1:3,5:6)){
  Stacked[[i]][Stacked[[i]]>=1]<-cat[i]
}

Stacked<-sprc(Stacked)%>% mosaic(.,fun="sum") %>% as.factor()

# Assign categories to stacks ----
levels(Stacked) <- cls #set factors in stacked

# Write out raster -----
writeRaster(Stacked,paste("Potential_Extrapolation_by_Class_",b,"_BAM_V6.tif", sep=""), overwrite=T)

## Survey locations:
#loc <- visit %>% 
#  dplyr::select(id, lat, lon) %>% 
#  unique() %>% 
#  st_as_sf(coords=c("lon", "lat"), crs=4326, remove=FALSE) %>% 
#  st_transform(crs=5072)

#14. Save final plots ---------

png(paste("ExtrapolationClass",b,".png", sep=""),width =15,
    height=5,units = "in",res = 1200)
plot(Stacked, col=c32, all_levels=T, cex=1, main=paste("Extrapolation: all classes",b, sep=" "))
plot(ocean$geometry, lwd=0.1,col="#ccdce3", border = "#436175", add=T)
plot(water$geometry,lwd=0.05,col="white",border = "#436175", add=T)
#plot(loc$geometry,col=Pnt_col,cex=0.001, add=T)
plot(bcr_s$geometry,lwd=0.2, border="#434162",add=T)
dev.off()
gc()  
cat("Finished", b)} #end of year

################## End of Code ######################
