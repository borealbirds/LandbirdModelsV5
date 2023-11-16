brt.tune <- function(i){
  
  #1. Get model settings---
  bcr.i <- loop$bcr[i]
  spp.i <- loop$spp[i]
  boot.i <- loop$boot[i]
  lr.i <- loop$lr[i]
  id.i <- loop$id[i]
  
  #2. Get visits to include----
  visit.i <- visit[bootstraps[[bcr.i]][,boot.i+1]]
  
  #3. Get response data (bird data)----
  bird.i <- bird[bird$id==visit.i$id, spp.i]
  
  #4. Get covariates----
  #TO DO: UPDATE THIS!!!!
  covlist.i <- covlist[covlist$id==bcr.i,]
  cov.i <- visit.i[,colnames(visit.i) %in% covlist.i==TRUE]
  
  #5. Put together data object----
  dat.i <- rbind(bird.i, cov.i)
  
  #6. Get offsets----
  off.i <- offsets[offsets$id==visit.i$id, spp.i]
  
  #7. Run model----
  m.i <- dismo::gbm.step(data=dat.i,
                         gbm.x=c(2:ncol(dat.i)),
                         gbm.y=1,
                         offset=off.i,
                         tree.complexity = id.i,
                         learning.rate = lr.i,
                         family="poisson")
    
}