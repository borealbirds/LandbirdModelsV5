library(gbm)

root <- "/home/mannfred/scratch"
gbm_objs <- list.files(file.path(root, "CONW", "BCR_81"), full.names = TRUE) 


for (i in 1:length(gbm_objs)) {
  result <- suppressWarnings(try(load(gbm_objs[i]), silent = TRUE))  # Capture result and suppress warnings
  
  if (inherits(result, "try-error") || length(result) == 0) {  # Check if load failed or resulted in no objects
    message("Could not load gbm_obj at index ", i, ": ", gbm_objs[i])
  } else {
    message("Successfully loaded gbm_obj at index ", i)
  }
}

# if corrupted, load and resave with `version=2`
root <- "C:/Users/mannf/Downloads"
load(file.path(root, "CONW", "BCR_81", "gnmboot-CONW-BCR_81-18.RData"))
save(out, file=file.path(root, "CONW", "BCR_81", "gnmboot-CONW-BCR_81-18.RData"), version = 2)

