### this script reconstructs a data object `dat.i` corresponding to a bootstrap object `b.i`
### this data object can be used for estimating 2-way interactions via `gbm::interact.gbm`


library(tidyverse)
library(gbm)

# load large data object (333MB)
root1 <- "G:/Shared drives/BAM_NationalModels4/NationalModels4.0/"
load(file.path(root1, "data", "BAMdb-GNMsubset-2020-01-08.RData"))

# load bootstrap models 
root2 <- "G:/Shared drives/BAM_NationalModels4/NationalModels4.0/Feb2020/out/boot"
dirs <- c("CAWA/BCR_60")
gbm_objs <- unlist(lapply(file.path(root2, dirs), list.files, full.names = TRUE)) 
 
i <- 1  
load(file=gbm_objs[i])

# define species and region
spp <- attr(out, "__settings__")$species
bcr <- attr(out, "__settings__")$region

# create data subset for BCR unit
ss <- dd[, bcr] == 1L

# reconstruct DAT object (see: https://github.com/borealbirds/GNM/blob/master/R/04-boot.R)
DAT <- data.frame(
  count = as.numeric(yy[ss, spp]),  
  offset = off[ss, spp],
  cyid = dd$cyid[ss],
  YEAR = dd$YEAR[ss],
  ARU = dd$ARU[ss],   
  dd2[ss, out$var.names[-c(1:2)]]  #already have `YEAR` and `ARU`
)

# keep just one observation (row) for each unique "cell x year" combination
DAT <- DAT[sample.int(nrow(DAT)),]
DAT <- DAT[!duplicated(DAT$cyid),]


# subsample to speed up `interact.gbm`
set.seed(1)
if (sum(DAT$count) < 1) {
  pts[[interaction_index]] <- NA
} else {
  pts[[interaction_index]] <- 
    interact.gbm(x = out, 
                 data = dplyr::slice_sample(DAT[,-(1:3)], prop=0.10), 
                i.var = c(50,51))
} 



