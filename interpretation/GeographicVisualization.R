##########################################
# title: Plotting Mosaicked Outputs - working draft
# author: Anna Drake
# date: Sept 9, 2024
##########################################
# Script plots the mean output and SD of mosaicked V5 predictions by species
# Post-processing uses an 'Absence/Unk' raster underlay + a range extent shapefile (from Probable_Range_Extent.R)
# Predictions are averaged and then masked to the range extent. 
# Areas where we could infer absence from P/A data are coded 0 outside the range extent. 
# Where insufficient data available (<450 records within 250km), areas are coded NA
# Testing Canada Lambert projection for visualization

# Preamble -----
require(terra)
require(dplyr)
require(sf)
library(RColorBrewer)
require(data.table)

# Set root ----
#root1<-"G:/Shared drives/BAM_NationalModels5/output/mosaics/predictions"
root1<-"G:/Shared drives/BAM_NationalModels5/output/archive/mosaic_inf/predictions"
root <- "G:/Shared drives/BAM_NationalModels5"

# Custom pixel-level trend function -----
reg <-  function(val) { 
  if(sum(val)==0){NA} # return NA where no birds
  else{
  if (length(subset(val,val>=0))>2) { # if at least 3 predictions in stack
      start<-mean(val[1:10]) #calculate 1985 mean as baseline
      val<-val/start #set value as proportion of baseline
      rg<-lm(x~y,data = data.frame(x = val, y=rep(1:8, each=10))) #run lm (raster stacking is 10 boots x year)
  if (base::summary(rg)$coefficients[2,4]<=0.05) { # and if slope is significant
      base::summary(rg)$coefficients[2,1] # return slope
}
  else{NA} # otherwise return NA where slope is NS
  }

# Colour formatting --- 
colfunc<-colorRampPalette(c("lightyellow","darkgreen","navy")) # range colours
Low<-colorRampPalette(c("darkred","pink")) #negative trend colours
High<-colorRampPalette(c("lightblue","navy")) #positive trend colours

# BCR  ---
bcr<- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel.shp")) %>%
  dplyr::filter(subUnit!=1) %>% # remove BCR 1
  #st_transform(crs=crs(trend))
  st_transform(crs=crs("EPSG:3348"))

# Ocean layer ---
ocean<-read_sf("C:/Users/andrake/Downloads/ne_50m_ocean/ne_50m_ocean.shp")
ocean<- ocean %>% st_transform(crs("EPSG:3348")) %>% st_crop(bcr)

# Lakes layer ---
water<-read_sf("G:/Shared drives/BAM_NationalModels5/Regions/Lakes_and_Rivers/hydrography_p_lakes_v2.shp")
water<- water %>% st_transform(crs("EPSG:3348")) %>% st_crop(bcr) %>% dplyr::filter(TYPE==16|TYPE==18)

# Sufficient data layer ---
layer<-rast(paste0("G:/Shared drives/BAM_NationalModels5/MosaicWeighting/Range/","ABS_NA_Layer450.tiff"))
layer[layer<=0]<-(-0.01)
layer[layer>0]<-0

# Define crs and limits by underlying layer ---
ext<-ext(layer)
crs<-crs(layer)

### Bring in predictions -----
list<-dir(root1)
spp<-unique(substr(list,1,4))
yr<-c(1985,1990,1995,2000,2005,2010,2015,2020)

for(j in c(1:length(spp))) {

#1. Import species range limit - using 0.2% occurrence ---
range<-st_read(paste0("G:/Shared drives/BAM_NationalModels5/MosaicWeighting/Range/ShapeFile/",spp[j],"_rangelimit.shp"))
range<-range%>%dplyr::filter(Limit=='0.2% limit')%>% st_transform(crs)

# Create alternate projection for plotting ---
plotlimit<-st_transform(range,crs(bcr))

#2. Import all years to standardize colour range ---
mean<-sd<-list()
upper<-c()

for (i in c(1:length(yr))){ #need to write in bootstrap flag for <10 layers
lyr<-list%>%.[. %like% spp[j]]%>%.[.%like% yr[i]]  

if (length(lyr)<10){
  print(paste0("Missing bootstraps for ", spp[j], yr[i]))
}

#else{
stk<-rast(file.path(root1,lyr))
stk[is.infinite(stk)]<-NA
mean[[i]]<-app(stk, fun="mean",na.rm=TRUE) %>% project(crs) %>%extend(ext) 
sd[[i]]<-app(stk, fun = "sd", na.rm=TRUE)
upper[i]<-global(mean[[i]], quantile, probs=.95, na.rm=TRUE) #95% ile
}
#}
  
#3. Set up un-modified raster stack for trend analysis, below----
trend<-rast(mean)
  
#4. Use the highest 95% value in the time series to standardize range ---
top<-as.numeric(max(unlist(upper)))
width<-(top+0.01)/14 #chop up range into 14 segments
breaks <- c(-0.01,seq(0,top,width),1) #16 categories incl. NA and >top category

# Set species-specific colour scale ----
colours<-colfunc(15) 
cpal <- c('grey80',colours)

#5. Post-process and run off each year -----
for (i in c(1:length(yr))){
mean[[i]][mean[[i]]<0.0001]<-0 #<0.01 bird km2 = 0; could modify this threshold
mean[[i]][mean[[i]]>=top]<-top

#temporary fix as some predictions don't reach the range limit - this will be corrected at mosaic stage ---
layer1<-mask(layer,mean[[i]], inverse=T)
out<-sum(layer1,mean[[i]], na.rm=T)

#cut to range extent ----
layer2<-mask(layer,range[2], inverse=T) 
out2<-mask(out,range[2])

#join to background and reproject to Canada Lambert ----
final<-sum(layer2,out2, na.rm=T)
final<-terra::project(final,crs("EPSG:3348"))

#6. save out image ----
png(paste("Mean_",spp[j],"_",yr[i],".png", sep=""),width=12.2,
    height=9,units = "in",res = 1200)

plot(ch, zlim = c(0, top),breaks=breaks, axes=F, col=cpal, main=paste0(spp[j], " males/ha ", yr[i]), cex=1)
#plot(ch,col = terrain.colors(30, rev=-1), cex=1) # alternate colour scheme
plot(ocean$geometry, lwd=0.1,col="#ccdce3", border = "#436175", add=T)
plot(water$geometry,lwd=0.05,col="#ccdce3",border = "#436175", add=T)
plot(plotlimit$geometry,lwd=0.5, border="darkred",add=T)
plot(bcr$geometry,lwd=0.6, border="grey50",add=T)
dev.off()
  }

#7. Calculate pixel-level trends ----

slope <- app(trend, fun=reg)  # note this has a 5-year interval, so result is proportion change/5 yrs

# reproject to Canada Lambert ----
slope_map<-terra::project(slope_sm,crs("EPSG:3348"))  

# Center range on 0 for visualization -----
range<-global(slope_map, quantile, probs=c(0.01,0.99), na.rm=TRUE)%>%unlist()%>%as.vector()
step<-round((range[2]-range[1])/10,4) # assign interval of 10 breaks

lower<-round(seq(as.integer(range[1]/step)*step,(-(step)),step),4)
upper<-round(seq(step,as.integer(range[2]/step)*step,step),4)

diverge<-c(Low(length(lower)),"white",High(length(upper)))
brk<-c(lower,0,upper)
lim<-c((min(brk))-(step*0.5),(max(brk))+(step*0.5))
plg <-list(loc = "bottom",at = brk, cex=2)

# R plot drops out-of-range values so need to set to max/min -------
slope_map[slope_map<lim[1]]<-lim[1]
slope_map[slope_map>lim[2]]<-lim[2]

# 8. Plot trend output ------
  png(paste0("1985-2020_PredictedTrend_1k_",spp[j],".png"),width=15,
      height=10, units = "in",res = 1200)
  plot(bcr$geometry, col="grey80", main=paste0("Trend/ 5yrs ", spp[j]))
  plot(ocean$geometry, lwd=0.1,col="black", border = "black", add=T)
  plot(slope_map,col=diverge, range=lim, plg=plg, add=T)
  plot(water$geometry,lwd=0.05,col="black",border = "black", add=T)
  plot(plotlimit$geometry,lwd=0.7, border="darkgreen",add=T) 
  plot(bcr$geometry,lwd=0.7, border="black", add=T)
  dev.off()
  #} #close spp loop
##################### End of Code ###################

  #testing lower resolution - need to think on the correct approach
  #should we run regression on 9 values per layer (possibly pseudo-replication?) or average prior?
  # Produce a reference layer ----
  r2 <- trend[[1]]
  res(r2) <- 3000 #set to 3km2
  # Average boots over 9 pixels ----
  trend_sm<-resample(trend,method="average",r2) 
  # run regression on lower resolution result ----
  slope_sm<-app(trend_sm, fun=reg)
