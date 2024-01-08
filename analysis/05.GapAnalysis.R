# ---
# title: National Models 5.0 - Gap Analysis
# author: Anna Drake
# created: November 2023
# ---

#NOTES################################

# This code uses dsextra to quantify extrapolation
# The package calls the ExDet function which needs it's
# tolerances adjusted to avoid singular matrices (rounding low values down)

# The first analysis uses all sampling data (-) to assess coverage of the predictor 
# variable environment
# The second analysis will look at temporal variation in coverage (assuming year * predictor
# interactions occur). We use sampling of 10-year windows straddling 5-year 
# predictor rasters (2000,2005,2010,2015,2020).

# extra function
`%notin%` <- Negate(`%in%`)

# NAD83(NSRS2007)/Conus Albers projection (epsg:5072)
crs<-"+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"

#PREAMBLE############################

#1. Load packages----

# Import dsmextra
  #if (!require("remotes")) install.packages("remotes")
  #remotes::install_github("densitymodelling/dsmextra")

pacman::p_load(dsmextra, sp, tidyverse, sf, ggplot2, lattice,
               rlang, rgeos, tools, dplyr, plyr, terra, smoothr, utils, 
               knitr, ecospat, magrittr, utils, readr, mapview, htmlwidgets)

# Tweak ExDet function in dsmextra: increase the tolerance of the Mahalanobis function 

#getAnywhere(ExDet) 

ExDet<-function (ref, tg, xp) 
{
  tg <- as.matrix(tg)
  ref <- as.matrix(ref)
  a <- apply(ref, 2, min, na.rm = TRUE)
  b <- apply(ref, 2, max, na.rm = TRUE)
  minref <- matrix(a, nrow = nrow(tg), ncol = ncol(tg), byrow = TRUE)
  maxref <- matrix(b, nrow = nrow(tg), ncol = ncol(tg), byrow = TRUE)
  nt1.df <- data.frame(apply(array(data = c(tg - minref, maxref - 
                                              tg, rep(0, nrow(tg) * ncol(tg))), dim = c(dim(tg), 3)), 
                             c(1, 2), min)/(maxref - minref))
  names(nt1.df) <- xp
  nt1 <- rowSums(nt1.df)
  mic_nt1 <- apply(nt1.df, 1, FUN = function(x) base::which.min(x))
  univ.rge <- which(nt1 == 0)
  mic_nt1[univ.rge] <- NA
  tg.univ <- matrix(tg[univ.rge, ], ncol = ncol(tg))
  aa <- apply(ref, 2, mean, na.rm = TRUE)
  bb <- stats::var(ref, na.rm = TRUE)
  mah.ref <- stats::mahalanobis(x = ref, center = aa, cov = bb, tol=1e-20) # increase tolerance
  mah.pro <- stats::mahalanobis(x = tg.univ, center = aa, cov = bb, tol=1e-20)# increase tolerance
  mah.max <- max(mah.ref[is.finite(mah.ref)])
  nt2 <- mah.pro/mah.max
  nt1[univ.rge] <- nt2
  if (length(xp) == 1) 
    cov.combs <- matrix(1)
  if (length(xp) > 1) 
    cov.combs <- utils::combn(x = 1:ncol(tg.univ), m = length(xp) - 
                                1)
  cov.combs <- as.list(data.frame(cov.combs))
  if (length(xp) == 1) {
    cov.aa <- cov.combs %>% purrr::map(., ~apply(as.matrix(ref[, 
                                                               .]), 2, mean))
    cov.bb <- cov.combs %>% purrr::map(., ~var(as.matrix(ref[, 
                                                             .])))
  }
  else {
    cov.aa <- cov.combs %>% purrr::map(., ~apply(as.matrix(ref[, 
                                                               .]), 2, mean, na.rm = TRUE))
    cov.bb <- cov.combs %>% purrr::map(., ~var(as.matrix(ref[, 
                                                             .]), na.rm = TRUE))
  }
  if (nrow(tg.univ) < 2) {
    warning("Only one prediction point within analogue conditions. Mahalanobis distances cannot be calculated.")
    mah_nt2 <- vector(mode = "list", length = length(cov.combs))
  }
  else {
    mah_nt2 <- purrr::pmap(.l = list(cov.combs, cov.aa, cov.bb), 
                           .f = function(a, b, c) stats::mahalanobis(x = as.matrix(tg.univ[, 
                                                                                           a]), center = b, cov = c, tol=1e-20))# increase tolerance
  }
  mah_nt2 <- mah_nt2 %>% purrr::set_names(., xp)
  mah_nt2 <- mah_nt2 %>% purrr::map_df(., cbind)
  mah_nt2 <- as.matrix(mah_nt2)
  mic_nt2 <- 100 * (mah.pro - mah_nt2)/mah_nt2
  mic_nt2 <- apply(mic_nt2, 1, FUN = function(x) base::which.max(x))
  results <- tibble::tibble(ExDet = nt1, mic_univariate = mic_nt1, 
                            mic_combinatorial = NA)
  if (nrow(tg.univ) > 1) 
    results$mic_combinatorial[univ.rge] <- mic_nt2
  results <- results %>% dplyr::mutate(mic_combinatorial = ifelse(ExDet >= 
                                                                    0 & ExDet <= 1, NA, mic_combinatorial))
  results <- results %>% dplyr::mutate(mic = rowSums(.[2:3], 
                                                     na.rm = TRUE))
  return(results)
}


#2. Set root path for data on google drive ----
root <- "/Volumes/GoogleDrive/Shared drives/BAM_NationalModels/NationalModels5.0"

#3. Load data package with covariates and BCR lists ----
load(file.path(root, "Data", "04_NM5.0_data_stratify.R"))
rm(bird,birdlist,gridlist,offsets)

#clean sample covariates
cov_clean<-cov%>%dplyr::select(where(is.numeric))

#temp value correction for negative biomass values
cov_clean[c(10:23,30:53)][cov_clean[c(10:23,30:53)]<0]<- NA

#4. Get BCR boundaries ----
can <- read_sf(file.path(root, "Regions", "CAN_adm", "CAN_adm0.shp")) %>% 
  st_transform(crs=crs) # Canadian boundary

bcr.ca <- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel.shp")) %>% 
  dplyr::filter(subUnit!=1) %>% # remove BCR 1
  st_transform(crs=crs) %>%
  st_intersection(can) %>% 
  mutate(country="ca") # Restrict to Canadian BCRs

#5. Produce look-up table: "BCR label"-"spatial sf subunit"   ----
BCR_Names <- colnames(bcrlist)%>%.[grep('ca',.)]
BCR_lookup <- data.frame(matrix(nrow=length(BCR_Names),ncol=0))
BCR_lookup$Names<-BCR_Names

BCR_lookup$subUnit<-as.numeric(gsub("\\D", "", BCR_Names))

Left<-Right<-subset(BCR_lookup,BCR_lookup$subUnit>100) #set both component BCRs to link to merged BCR
Left$subUnit<-as.numeric(substr(Left$subUnit,1,2))
Right$subUnit<-as.numeric(substr(Right$subUnit,3,4))

BCR_lookup<-rbind(subset(BCR_lookup,BCR_lookup$subUnit<100),Left,Right)
rm(Left,Right,BCR_Names)

#6. Load prediction raster look-up, get Canadian prediction variables ----

Canada <- readxl::read_excel(file.path(root,"NationalModels_V5_VariableList.xlsx"), 
                             sheet = "ExtractionLookup") %>%
                             dplyr::filter(Extent %in% c("Canada","global","north america","Canada-forested"))

Canada$Raster<-Canada$Label
Canada$Raster[Canada$Raster=="ERAPPTsmt_1km"]<-"ERAPPTsm_1km" # remove the "t" label for pulling raster (t vs t-1)
Canada$Raster[Canada$Raster=="ERATavesmt_1km"]<-"ERATavesm_1km" # remove the "t" label for pulling raster (t vs t-1)

Grouping<-unique(Canada$Category) #Groupings
Grouping<-Grouping[-c(1,3,6,7,9)] ##**dropping annual climate for now with labelling issue unresolved

#7. Functions to get and prep data for analysis -----

#(1) Function to extract raster list - returns matrix of raster pathways [1,] + variable names [2,]
Get_paths<- function(Grouping, YR){
  
  Path_list<-Name_list<-c()
  Group.j<-dplyr::filter(Canada,Canada$Category==Grouping)

for (k in (1:nrow(Group.j))){ # Get variables within each group
    
Period<-Static<-Name<-NULL
    
# Return the closest year to desired period
if(Group.j$TemporalResolution[k]=="match") {
  
  years<-dir(file.path(root,Group.j$PredictPath[k]))%>%
    stringr::str_sub(.,-8,-5)%>%as.numeric(.)%>%suppressWarnings()  #pull all years
    
  Match<-which(abs(years - YR) == min(abs(years - YR),na.rm=T))%>%
    years[.]%>%unique(.)  #ID closest year (minimum value when multiple options)

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
# uses "Get_paths" function to extract rasters for stacking

create_df <- function(Grouping,YR,bcr.i){
  
Stack<-list()  
e<-ext(bcr.i)

call<-BCR_lookup%>%dplyr::filter(subUnit==bcr.i$subUnit)%>%.$Names #match to bcr name
covariates<-covlist%>%dplyr::filter(bcr==call[1])%>%.[,.==TRUE]%>% colnames()#get co-variates, assumes that BCR covariates don't differ from pooled BCRs (81,82 and 41,42)

Var<-sapply(Grouping,Get_paths,YR=2000) #get raster paths

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
  
DF<-rast(Stack[1:length(Stack)]) %>% #turn into dataframe
    as.data.frame(., xy=TRUE) 

names(DF)<-c("x","y",Label)

return(DF)
}

#8. Compile raster lists for each category + early and late periods ----

BCR_list<-list() # for storing output by BCR
                   
# Analysis ----
for(i in 1:nrow(bcr.ca)){  # select BCR

#Buffer BCR  ---  
bcr.i <-  bcr.ca %>% 
    dplyr::filter(row_number()==i) %>% 
    st_buffer(100000)%>%  # Filter & buffer shapefile
    st_intersection(., can) #crop at Canadian border
  
# Extract BCR specific data frames for each group. Year=2000  ---  
Early_Dataframes<-purrr::map(Grouping,create_df, YR=2000, bcr.i=bcr.i,.progress = T)
#Late_Dataframes<-purrr::map(Grouping,create_df, YR=2020, bcr.i=bcr.i,.progress = T)

#Get BCR reference sample - pool samples for the pooled BCRs  ---  
call<-BCR_lookup%>%dplyr::filter(subUnit==bcr.i$subUnit)%>%.$Names #match to bcr.i to bcr name
ids <- bcrlist %>%.[, c("id",call)]%>%subset(.,.[,2]==TRUE)%>%.$id

if(length(call)==2) {
id2<- bcrlist %>%.[, c("id",call)]%>%subset(.,.[,3]==TRUE)%>%.$id  
ids<-unique(c(ids,id2))  }

# Get sample data for BCR ---
sample.i<-cov_clean%>%subset(.,.$id %in% ids) 
sample.i<-sample.i %>%.[,colSums(.,na.rm=T)!=0] # remove covariates with no data

#Convert sample locations to SPDF (for Leaflet plots)  ---  
pnts<-visit%>% subset(.,.$id %in% ids)%>%.[10:11]%>% #sample lon, lat
  st_as_sf(., coords = c("lon", "lat"),crs = 4326)%>% 
  st_transform(crs=crs) %>%
  as(.,"Spatial") 

# Calculate extrapolation ---
EarlyExt<-LateExt<-list()

# Run through variable grouping for Early period -------

for (j in 1:length(Early_Dataframes)){ #run through each group

target=Early_Dataframes[[j]]

# Ensure match training-predictor variables
sample<-sample.i%>% dplyr::select(one_of(names(target))) #drop any in ref not present in target
target<-target%>%dplyr::select(c("x","y",names(sample))) #drop any in target not present in ref


xp = names(sample)

sample[sapply(sample, is.infinite)] <- NA #covariance returns infinite...

EarlyExt[[j]] <- compute_extrapolation(samples = sample,
              covariate.names = xp,
              prediction.grid = target,
              coordinate.system = crs,
              verbose=F)

#9. Plot output ---

Map1M<-map_extrapolation(map.type="mic",sightings= pnts, 
                         extrapolation.object=EarlyExt[[j]],
                         verbose=F)

saveWidget(Map1M, file=paste("BCR_", bcr.i$subUnit,
                             "_2000_",
                             Grouping[j],
                             "_MIC",".html", sep=""))

if(j!=4){ # skip issue with plotting road data: unable to find an inherited method for function ‘Which’ for signature ‘"logical"’
Map1E<-map_extrapolation(map.type="extrapolation",sightings= pnts, 
                         extrapolation.object=EarlyExt[[j]],
                         verbose=F)

saveWidget(Map1E, file=paste("BCR_", bcr.i$subUnit,
                             "_2000_",
                             Grouping[j],
                             "_Ext",".html", sep=""))
}
  
} # end of Grouping
  
BCR_list[[i]]<- EarlyExt
saveRDS(BCR_list, file = "/Volumes/GoogleDrive/Shared drives/BAM_NationalModels/NationalModels5.0/GapAnalysis/Extrapolation_object.rds", overwrite=T)

} # end of BCR
                   
#################### End of Code ######################
