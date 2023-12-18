# ---
# title: National Models 5.0 - Gap Analysis
# author: Anna Drake
# created: November 2023
# ---

#NOTES################################

# This code uses dsextra to quantify extrapolation
# The package calls the ExDet function which needs it's
# tolerances adjusted to avoid singular matrices (it rounds low values down to 0)

# The first analysis uses all sampling data (-) to assess coverage of the predictor 
# variable environment
# The second analysis looks at temporal variation in coverage (assuming year * predictor)
# interactions occur. We use sampling of 10-year windows straddling 5-year 
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

pacman::p_load(dsmextra, sp, tidyverse, sf, ggplot2,
               rlang, tools, dplyr, plyr,terra, smoothr, utils, 
               knitr, ecospat, magrittr, utils, readr)

# Tweak ExDet function in dsmextra:
# increase the tolerance of the Mahalanobis function, 
# Remove conditional computation for single variables: causing combitorial calculation to fail. Haven't had time to troubleshoot, so I've just cut it out... 

#getAnywhere(ExDet) 
ExDet<-function (ref, tg, xp) {
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
  mah.ref <- stats::mahalanobis(x = ref, center = aa, cov = bb, tol=1e-20 ) # increased tolerance
  mah.pro <- stats::mahalanobis(x = tg.univ, center = aa, cov = bb, tol=1e-20 )# increased tolerance
  mah.max <- max(mah.ref[is.finite(mah.ref)])
  nt2 <- mah.pro/mah.max
  nt1[univ.rge] <- nt2
  
  #if (length(xp) == 1)  # treating xp as never == 1 as throwing error 
  #  cov.combs <- matrix(1)
  #if (length(xp) > 1) 
  cov.combs <- utils::combn(x = 1:ncol(tg.univ), m = length(xp) - 1)
  cov.combs <- as.list(data.frame(cov.combs))
  #if (length(xp) == 1) {
  #  cov.aa <- cov.combs %>% purrr::map(., ~apply(as.matrix(ref[,.]), 2, mean))
  #  cov.bb <- cov.combs %>% purrr::map(., ~var(as.matrix(ref[,.])))
  #}
  #else {
  cov.aa <- cov.combs %>% purrr::map(., ~apply(as.matrix(ref[,.]), 2, mean, na.rm = TRUE))
  cov.bb <- cov.combs %>% purrr::map(., ~var(as.matrix(ref[,.]), na.rm = TRUE))
  #}
  if (nrow(tg.univ) < 2) {
    warning("Only one prediction point within analogue conditions. Mahalanobis distances cannot be calculated.")
    mah_nt2 <- vector(mode = "list", length = length(cov.combs))
    
  } else {
    mah_nt2 <- purrr::pmap(.l = list(cov.combs, cov.aa, cov.bb), 
                           .f = function(a, b, c) stats::mahalanobis(x = as.matrix(tg.univ[,a]), center = b, cov = c, tol=1e-20 ))
  }
  mah_nt2 <- mah_nt2 %>% purrr::set_names(., xp)
  mah_nt2 <- mah_nt2 %>% purrr::map_df(., cbind)
  mah_nt2 <- as.matrix(mah_nt2)
  mic_nt2 <- 100 * (mah.pro - mah_nt2)/mah_nt2
  mic_nt2 <- apply(mic_nt2, 1, FUN = function(x) base::which.max(x))
  results <- tibble::tibble(XDet = nt1, mic_univariate = mic_nt1, mic_combinatorial = NA) #changed ExDet to XDet as creates confusion with function
  if (nrow(tg.univ) > 1) 
    results$mic_combinatorial[univ.rge] <- mic_nt2
  results <- results %>% dplyr::mutate(mic_combinatorial = ifelse(XDet >=0 & XDet <= 1, NA, mic_combinatorial))
  results <- results %>% dplyr::mutate(mic = rowSums(.[2:3], na.rm = TRUE))
  return(results)
}

#2. Set root path for data on google drive----
root <- "/Volumes/GoogleDrive/Shared drives/BAM_NationalModels/NationalModels5.0"

#3. Load data packages with offsets and covariates ----
load(file.path(root, "Data", "04_NM5.0_data_stratify.Rdata"))

#4. Load prediction raster lookup ----
Lookup <- readxl::read_excel(file.path(root, "NationalModels_V5_VariableList.xlsx"), sheet = "ExtractionLookup")

Lookup$Label[3]<-"ERAPPTsm_1km" # remove the "t" label for this (t vs t-1)
Lookup$Label[6]<-"ERATavesm_1km"# remove the "t" label for this (t vs t-1)

#list for stripping paths down to just variable names
cleanup = unique(paste(c(Lookup$PredictPath,"/",".tif"), collapse = "|")) 

#5. Get Canadian prediction variables for each variable grouping ----
Canada<-Lookup %>% dplyr::filter(Extent %in% c("Canada","global","north america","Canada-forested"))
Grouping<-unique(Canada$Category) #Groupings

#6. Get BCR boundaries----

can <- read_sf(file.path(root, "Regions", "CAN_adm", "CAN_adm0.shp")) %>% 
  st_transform(crs=5072) # Political boundaries

bcr <- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel.shp")) %>% 
  dplyr::filter(subUnit!=1) %>% 
  st_transform(crs=5072) #BCR shapefiles & remove subunit 1

bcr.ca <- bcr %>% 
  st_intersection(can) %>% 
  mutate(country="ca") # Restrict to Canadian BCRs

#7. Functions to get and prep data for analysis -----

#(1) Function to extract raster list - returns matrix of early [1,] and late [2,] variables
Get_paths<- function(Grouping){
  Early_list<-Late_list<-c()
  Group.j<-dplyr::filter(Canada,Canada$Category==Grouping)
  
  for (k in (1:nrow(Group.j))){ # Get variables within group
    
    Early<-Late<-Static<-NULL
    
    if(Group.j$TemporalResolution[k]=="match"){
      Early<-paste(Group.j$PredictPath[k],"/",Group.j$Label[k],"_",Group.j$Early[k]+Group.j$YearMatch[k],".tif",sep="")
      Late<-paste(Group.j$PredictPath[k],"/",Group.j$Label[k],"_",Group.j$Late[k]+Group.j$YearMatch[k],".tif",sep="")
      
    }else{
      Static<-paste(Group.j$PredictPath[k],"/",Group.j$Label[k],".tif",sep="")
    }
    
    Early_list<-c(Early_list,Early,Static) #early and static list
    
    if (unique(Group.j$Category) %notin% c("Climate Normals","Topography","Wetland"))
    Late_list<-c(Late_list,Late,Static) #late list only where temporally variable
    
  } # close group
  
  return(list(Early_list,Late_list))
}

#(2) Function to stack rasters, crop by BCR, and convert to data frame 
create_df <- function(Var){
  
  Names<-str_remove_all(Var,cleanup)
  Stack<-list()  
 
if(!is.null(Var)) {  
  for (k in 1:length(Var)) { #stack rasters
   Stack[[k]]<-terra::rast(file.path(root, Var[k])) %>% 
        project(.,crs) %>%
        crop(x = ., y = ext(bcr.i)) %>%
        mask(x = ., mask = bcr.i)
    }
    
   DF<-rast(Stack[1:length(Stack)]) %>% #turn into dataframe
       as.data.frame(., xy=TRUE) %>% 
       .[rowSums(.[,-c(1:2)], na.rm=T) > 0,] #get rid of fully masked rows
    
    names(DF) <- c("x","y",Names)
    
}else{
    DF<-NULL 
  }
  return(DF)
}

#8. Compile raster lists for each category + early and late periods----
Variables<-mapply(Grouping,Get_paths)

# Split early and late periods...
Early<-Variables[1,]
Late<-Variables[2,]

#9. Run analysis across each BCR ----

BCR<-list()

for(i in 1:nrow(bcr.ca)){  

 bcr.buff <- bcr.ca %>% 
    dplyr::filter(row_number()==i) %>% 
    st_buffer(100000) # Filter & buffer shapefile

 bcr.i <- st_intersection(bcr.buff, can) #crop at Canadian border

#Extract data frame for each grouping within BCR
Early_Dataframes<-purrr::map(Early,create_df)
Late_Dataframes<-purrr::map(Late,create_df)

# Calculate extrapolation
EarlyExt<-LateExt<-list()

# Early period
for (j in 1:length(Early_Dataframes)){

#create a toy sample
ref=Early_Dataframes[[j]][c(1:1000,4000:5000,70000:80000,100000:200000),] 

#eliminate variables with no training data
nodata<-colSums(ref, na.rm=T)%>%subset(.==0)%>% names(.) #which variables have no data

target=Early_Dataframes[[j]]%>%.[, which(names(.) %notin% nodata)] #drop those in target
ref=ref[, which(names(ref) %notin% nodata)] #drop those in training

xp = names(ref[,-c(1,2)]) 
samples = ref[,-c(1,2)]

EarlyExt[[j]] <- compute_extrapolation(samples = samples,
              covariate.names = xp,
              prediction.grid = target,
              coordinate.system = crs)
} # end of early variable group

# Late period
for (j in 1:length(Late_Dataframes)){
   
   #create a toy sample
   ref=Late_Dataframes[[j]][c(1:1000,4000:5000,70000:80000,100000:200000),] 
   
   #eliminate variables with no training data
   nodata<-colSums(ref, na.rm=T)%>%subset(.==0)%>% names(.) #which variables have no data
   
   target=Late_Dataframes[[j]]%>%.[, which(names(.) %notin% nodata)] #drop those in target
   ref=ref[, which(names(ref) %notin% nodata)] #drop those in training
   
   xp = names(ref[,-c(1,2)]) 
   samples = ref[,-c(1,2)]
   
   LateExt[[j]] <- compute_extrapolation(samples = samples,
                                          covariate.names = xp,
                                          prediction.grid = target,
                                          coordinate.system = crs)
 } # end of late variable group
 
}

BCR[[i]]<- list(EarlyExt,LateExt)
names(BCR)<-bcr.ca$subUnit
#} # end of BCRs


#11. Plot output
pal.a <- colorRampPalette(c("lightgreen","darkgreen"))
pal.c <- colorRampPalette(c("lightpurple","violet"))
pal.u <- colorRampPalette(c("yellow","red"))
dev.off()
#Early 
par(mfrow=c(3,3))
i=1
for (i in 1:length(Grouping)){
  u<-EarlyExt[[1]]$rasters$ExDet$univariate
  c<-EarlyExt[[1]]$rasters$ExDet$combinatorial
  a<-EarlyExt[[1]]$rasters$ExDet$analogue
  
  if(!is.null(u)) terra::plot(u, main=paste(Grouping[i]," Early"), col=pal.u(300)) 
  if(!is.null(c)) terra::plot(c, col=pal.c(300), add=T) 
  if(!is.null(a)) terra::plot(a, col=pal.a(300), add=T) 
}
 
# Late
for (i in 1:length(Grouping)){
  u<-VarGroup[[i]][[2]]$rasters$ExDet$univariate
  c<-VarGroup[[i]][[2]]$rasters$ExDet$combinatorial
  a<-VarGroup[[i]][[2]]$rasters$ExDet$analogue
  
  if(!is.null(u)) terra::plot(u, main=paste(Grouping[i]," Late"),col=pal.u(300)) 
  if(!is.null(c)) terra::plot(c, col=pal.c(300), add=T) 
  if(!is.null(a)) terra::plot(a, col=pal.a(300), add=T) 
}

### end of code
