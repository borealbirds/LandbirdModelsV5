###############################################
# title: "Truncation values"
# author: Elly Knight
# date: March 13, 2026
##################################################

#This code 

#PREAMBLE############################

#1. Load packages----
library(tidyverse) #basic data wrangling
library(QPAD) #offsets

#2. Set root path for data on google drive----
root <- "G:/Shared drives/BAM_NationalModels5"

#3. Load dataset ----
load(file.path(root, "data", "04_NM5.0_data_stratify.Rdata"))

#4. Load QPAD ----
load_BAM_QPAD("3")

#GET SPECIES VALUES #########

#1. Get species list ----
spp <- colnames(birdlist |> dplyr::select(-bcr))

#2. Get global max count ----
bird.all <- as.data.frame(as.matrix(bird)) |> 
  pivot_longer(ALFL:YTVI, names_to="spp", values_to="count")
countmax <- stats::quantile(bird.all$count, probs = 0.995, na.rm = TRUE) 

#3. Set up loop ----
q.out <- data.frame()
for(i in 1:length(spp)){
  
  spp.i <- spp[i]
  
  #4. Get species-specific truncation count ----
  #use varying qs depending on sppp
  bird.i <- as.integer(bird[,spp.i])
  countmax.i <- stats::quantile(bird.i, probs = 0.99, na.rm = TRUE)
  thresh.i = 0.99
  if(countmax.i < countmax){
    countmax.i <- stats::quantile(bird.i, probs = 0.995, na.rm = TRUE)
    thresh.i = 0.995
  }
  if(countmax.i < countmax){
    countmax.i <- stats::quantile(bird.i, probs = 0.999, na.rm = TRUE)
    thresh.i = 0.999
  }
  if(countmax.i < countmax){
    countmax.i <- stats::quantile(bird.i, probs = 0.9999, na.rm = TRUE)
    thresh.i = 0.9999
  }
  
  #5. Get null QPAD correction ----
  if(spp.i=="CAJA"){
    cf0 <- exp(unlist(coefBAMspecies("GRAJ", 0, 0)))
  } else {
    cf0 <- exp(unlist(coefBAMspecies(spp.i, 0, 0)))
  }
  
  off.i <- exp(cf0[1])*exp(cf0[2])
  
  #6. Get density truncation ----
  q99 <- countmax.i/off.i
  
  #7. Output ----
  q.out <- rbind(data.frame(spp = spp.i,
                            thresh = thresh.i,
                            countmax = countmax.i,
                            off = off.i,
                            q = q99),
                 q.out)
  
  cat(i, " ")
  
}

#7. Save ----
write.csv(q.out, file.path(root, "data", "Lookups", "SpeciesPredictionTruncationValues.csv"), row.names = FALSE)
