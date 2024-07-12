#TO DO: MASK OUT WATER########
#TO DO: ADD SAMPLING PER SPECIES#######
#TO DO: MEANS AND SD######
#LIST OF PRODUCTS PER SPECIES:
#2 extents: boreal, regional
#Products in both: mean, sd, extrapolation, sampling


#LIST OF SAMPLING PRODUCT:
#All the extrapolation
#Sampling density
#Distance to sampling









#MOVE THIS TO LATER#######
#MASKING####

#1. Get the water layer----
water <- read_sf(file.path(root, "Regions", "Lakes_and_Rivers", "hydrography_p_lakes_v2.shp")) |> 
  dplyr::filter(TYPE %in% c(16, 18)) |> 
  st_transform(crs=crs(bcr.out)) |> 
  vect()

#2. Set up loop----
for(i in 1:nrow(files.stack)){
  
  #3. Get the stack----
  stack.i <- rast(file.path(root, "PredictionRasters", "SubunitStacks", "Unmasked", files.stack$file[i]))
  
  #4. Mask it----
  mask.i <- mask(stack.i, water, inverse=TRUE)
  
  #5. Save it----
  terra::writeRaster(mask.i, file.path(root, "PredictionRasters", "SubunitStacks", "Masked", paste0(files.stack$bcr[i], "_", files.stack$year[i], ".tif")), overwrite=TRUE)
  
  #6. Report----
  print(paste0("Finished raster ", i, " of ", nrow(files.stack)))
  
  
}