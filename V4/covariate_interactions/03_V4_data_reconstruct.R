### this script reconstructs a data object `dat.i` corresponding to a bootstrap object `b.i`
### this data object can be used for estimating 2-way interactions via `gbm::interact.gbm`


library(tidyverse)
library(gbm)

# load large data object (333MB)
root1 <- "G:/Shared drives/BAM_NationalModels4/NationalModels4.0/"
load(file.path(root1, "data", "BAMdb-GNMsubset-2020-01-08.RData"))

# load bootstrap models 
# NOTE: need to create `DAT` for every bcr x spp tuple
root2 <- "G:/Shared drives/BAM_NationalModels4/NationalModels4.0/Feb2020/out/boot"
dirs <- "CAWA/BCR_60"
gbm_objs <- list.files(file.path(root2, dirs), full.names = TRUE)
 
i <- 1  
load(file=gbm_objs[i])

# define species and region
spp <- attr(out, "__settings__")$species
bcr <- attr(out, "__settings__")$region

# create an index indicating where in `dd` we will find the current BCR's data
# ss = subset
ss <- dd[, bcr] == 1L

# reconstruct DAT object (see: https://github.com/borealbirds/GNM/blob/master/R/04-boot.R)
# `yy` contains the counts data, and is subset using `ss` and `spp` 
DAT <- 
  data.frame(
  count = as.numeric(yy[ss, spp]),  
  offset = off[ss, spp],
  cyid = dd$cyid[ss],
  YEAR = dd$YEAR[ss],
  ARU = dd$ARU[ss],   
  dd2[ss, out$var.names[-c(1:2)]]  #-c(1:2) because we already have `YEAR` and `ARU` from `dd`
)

# keep just one observation (row) for each unique "cell x year" combination
DAT <- DAT[sample.int(nrow(DAT)),]
DAT <- DAT[!duplicated(DAT$cyid),]



# TESTING
# using `slice_sample` to subsample and speed up `interact.gbm`
# testing on variables #50 and #51



set.seed(1)
pts <- list()

if (sum(DAT$count) < 1) {
  pts[[interaction_index]] <- NA
} else {
  pts[[interaction_index]] <- 
    interact.gbm(x = out, 
                 data = dplyr::slice_sample(DAT[,-(1:3)], prop=0.20), 
                i.var = c(50,51))
} 



