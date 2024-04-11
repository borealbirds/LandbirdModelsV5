#------------------------------# 
#  mTPI calculation for region #
# author: Anna Drake           #
# created: December 2022       #
#------------------------------#

# mTPI is a multi-scale topographic position value. Based upon the Global SRTM mTPI on GEE which does not cover northern Canada
# Values represent the absolute maximum difference between a pixel and the surrounding area at 5 radii (r=1.5km-115.5km)
# uses 1km resolution Earth Env median elevation layer. Based on USGS GMTED 2010 layer
#Citation: Amatulli, G., Domisch, S., Tuanmu, M.-N., Parmentier, B., Ranipeta, A., Malczyk, J., and Jetz, W. (2018) A suite of global, 
#cross-scale topographic variables for environmental and biodiversity modeling. Scientific Data volume 5, 
#Article number: 180040. DOI: doi:10.1038/sdata.2018.40

# Packages ----------
require(terra)

# Import and crop raster ------------
Global<-rast("/Users/owner/Downloads/elevation_1KMmd_GMTEDmd.tif")
e<-extent(-175,-52,35,71)
GlobalCrop<-crop(Global,e) #cut to Canada + northern USA

# function for moving window TPI ------------- 
tpiw <- function(x, w) {
  m <- matrix(1/(w^2-1), nc=w, nr=w)
  m[ceiling(0.5 * length(m))] <- 0
  f <- focal(x, m)
  x - f
}

# Mimic the radii used for SRTM mTPI: 115.8, 89.9, 35.5, 13.1, 5.6, 2.8, and 1.2km radius -----------
tpi3 <- tpiw(GlobalCrop, 3) #1.5km radius
tpi5 <- tpiw(GlobalCrop, 5) #2.5km
tpi11<- tpiw(GlobalCrop, 11) # 5.5
tpi27<- tpiw(GlobalCrop, 27) #13.5
tpi71<- tpiw(GlobalCrop, 71) #35.5
tpi179<-tpiw(GlobalCrop, 179) #89.5
tpi231<- tpiw(GlobalCrop, 231) #115.5

Stack<-c(tpi3,tpi5,tpi11,tpi27,tpi71,tpi179,tpi231)
rs_ord<-app(Stack, fun=function(X,na.rm){X[order(abs(X),decreasing=T)]})
mTPI<-rs_ord[[1]]

# Write out -----------------
writeRaster(mTPI,"mTPI_Canada.tif")
#----------------------------------------------------
#confirm pixel selection function works ------
#test<-rast(ncol=2,nrow=2)
#values(test)<-c(1,1,10,10)

#test2<-rast(ncol=2,nrow=2)
#values(test2)<-c(2,-2,-12,-9)

#testing<-c(test,test2)
#plot(testing)
#rs_ord<-app(testing, fun=function(X,na.rm){X[order(abs(X),decreasing=T)]})
#out<-rs_ord[[1]]
#plot(out) #correct
####### END OF CODE ###########
