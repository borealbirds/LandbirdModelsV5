### this script reconstructs a data object `dat.i` corresponding to a bootstrap object `b.i`
### this data object can be used for estimating 2-way interactions via `gbm::interact.gbm`


library(tidyverse)
library(gbm)

root <- "G:/Shared drives/BAM_NationalModels5"

# load large data object (380MB)
load(file.path(root, "data", "04_NM5.0_data_stratify.R"))

# 1. load a specific GBM object (b.i) by index (e.g., i = 1)
gbm_objs <- list.files(file.path(root, "output", "bootstraps"))[1]
i <- 1  # change this index to match the desired bootstrap
load(file=file.path(root, "output", "bootstraps", gbm_objs[i]))

# 2. extract the covariates directly from b.i
  # this stores the covariates used in the model

# 3. Extract the relevant data (dat.i) and offsets (off.i) for the loaded GBM model
brt_boot_extract <- function(b.i, visit.i) {
  
  # 1. Get the visits that were included in the model (already in visit.i)
  # visit.i is already provided when loading the file, so no need to regenerate
  
  # 2. Extract the response data (bird counts) for the visits in `visit.i`
  bird.i <- bird[as.character(visit.i$id), out.i$spp]  # Assuming 'response.name' holds the species column
  
  # 3. Extract the covariates using only those present in b.i$var.names
  cov.i <- cov[cov$id %in% visit.i$id, colnames(cov) %in% b.i$var.names]
  
  # 4. Get the year data for the selected visits
  year.i <- visit[visit$id %in% visit.i$id, "year"]
  
  # 5. Combine the data into a single dataset for GBM (dat.i)
  dat.i <- cbind(bird.i, year.i, cov.i) %>%
    rename(count = bird.i)
  
  # 6. Extract the offsets for the selected visits (if applicable)
  off.i <- offsets[offsets$id %in% visit.i$id, out.i$spp]
  
  # Return the dataset (dat.i) and offsets (off.i)
  return(list(dat.i = dat.i, off.i = off.i))
}

# 4. Run the extraction for the specific GBM model that you just loaded
extracted_data <- brt_boot_extract(b.i, visit.i)

# Now you have the extracted data (dat.i and off.i) for use in `gbm::interact.gbm`
dat.i <- extracted_data$dat.i
off.i <- extracted_data$off.i



