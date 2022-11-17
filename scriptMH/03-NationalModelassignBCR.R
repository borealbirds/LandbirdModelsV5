########################################################
###### Overlay PKEY locations with buffered BCR
########################################################
library(sf)
library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(rgdal)
library(sp)

#Set params
w <- "E:/MelinaStuff/BAM/NationalModelv4.1/runNationalModel/data"
out <- file.path(w,"outsample")


# projection 
LCC <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
LAEA <- CRS("+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
DD <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
TM10 <- CRS("+proj=tmerc +lat_0=0 +lon_0=-115 +k=0.9992 +x_0=500000 +y_0=0 +datum=NAD83 +units=m +no_defs")


##Read all Diana csv
PKEYannualcombo <-read.csv(file.path(out,"PKEYannualcombo.csv"))
  #head(PKEYannualcombo)
  #length(unique(PKEYannualcombo$SS))  #221735

SSannualcombo <-as.data.frame(unique(PKEYannualcombo[c("SS", "longitude","latitude")]))
  #head(SSannualcombo)
  #length(SSannualcombo$SS)

##############
# Check Map coordinate
##############
   # Map Point count
coordinates(SSannualcombo) <- c("longitude", "latitude")
proj4string(SSannualcombo) <- LAEA

pc <- st_as_sf(SSannualcombo, coords = c("longitude", "latitude"), crs = LAEA)
plot(st_geometry(pc))


  ## Map canada extent
f_canada <-"E:/MelinaStuff/BAM/NationalModelv4.1/runNationalModel/data/basemaps/CanadaLAEA.shp"
canada<- st_read(f_canada)
plot(st_geometry(canada), add = TRUE)

##############
# Extract BCR
##############
# Map BCR
bcrUnit <- c("bcr4", "bcr5","bcr60","bcr61","bcr70","bcr71","bcr80","bcr81","bcr82","bcr83","bcr9","bcr10","bcr11","bcr12",
             "bcr13","bcr14","bcrus2","bcrus4","bcrus5","bcrus10","bcrus11","bcrus12","bcrus13","bcrus14")

fo_bcr <-"E:/MelinaStuff/BAM/NationalModelv4.1/runNationalModel/data/basemaps/"

  ##Check
  #bcr4<- st_read(paste0(fo_bcr,"bcr4_100km.shp"))
  #bcr4_LAEA <-st_transform(bcr4, LAEA)
  ##plot(st_geometry(bcr4_LAEA), col="navy", add=TRUE)

# Create df version of SSannualcombo to append bcr identity
SSannual_df <- as.data.frame(SSannualcombo)

# Intersect pc with bcr
for (bcr in bcrUnit){
   poly <-st_read(paste0(fo_bcr,bcr,"_100km.shp"))
   poly_LAEA <-st_transform(poly, LAEA)
   poly_LAEA$bcr <-1
   
   # use OVER fct to clip. over needs SpatialPolygonsDataFrame format. Convert simple feature to SpatialPolygonsDataFrame.
   poly_sp <- as_Spatial(poly_LAEA)
   a.data <- over(SSannualcombo, poly_sp[,"bcr"])
   # If you don't want the NA
   a.data$bcr[is.na(a.data$bcr)] <- 0
   SSannual_df[bcr] <- a.data$bcr
   
}

## Check bcr assigment
ll<- lapply(bcrUnit, function(x){
        pc_all <- st_as_sf(SSannual_df, coords = c("longitude", "latitude"), crs = LAEA)
        plot(st_geometry(pc_all))

        bcr<- st_read(paste0(fo_bcr,x,"_100km.shp"))
        bcr_LAEA <-st_transform(bcr, LAEA)
        plot(st_geometry(bcr_LAEA), border= "red", add=TRUE)
        
        ptc <- SSannual_df[c("SS", "longitude","latitude", x)]
        pc_bcr <- na.omit(ptc)
        pc_map <- st_as_sf(pc_bcr, coords = c("longitude", "latitude"), crs = LAEA)
        plot(st_geometry(pc_map), col="blue", add=TRUE)
}) 

dd <-merge(PKEYannualcombo, SSannual_df, by = c("SS", "longitude","latitude"))

save(dd, file = file.path(out,"PKEYannualcombo_bcr.RData"))

#write.csv(SSannual_df, file="E:/MelinaStuff/BAM/NationalModelv4.1/runNationalModel/data/SSannualcombo_bcr.csv", row.names=FALSE)
temp_env <- new.env()
load("E:/MelinaStuff/BAM/NationalModelv4.1/runNationalModel/data/NatMod41-2021-12-10.RData", envir=temp_env)


