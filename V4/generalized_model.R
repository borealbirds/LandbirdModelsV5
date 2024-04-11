# ---
# title: National Models 4.0 - generalized modelling process
# author: Peter Solymos, Elly Knight
# updated: April 5, 2024
# ---

#NOTES#####

# This script contains the function used for running the boosted regression models used in Version 4 of BAM's national bird species abundance models. The script is generalized to demonstrate the process and can be rerun using the data 

#1. Load libraries----
library(mefa4)
library(gbm)

#2. Load data----
load("NatModV4-2024-04-05.Rdata")

# Data contains:
# dd: visit x covariate dataframe
# yy: sparse matrix of visit x species abundances
# off: visit x species matrix of statistical offsets

# Row numbers of each object are unique primary keys that are matched across objects 

#3. List variables to be used as predictors----
# you can vary this by region
cn <- c("TPI", "TRI", "slope", "roughness",
  "ROAD", "AHM", "bFFP", "CMD", "DD_0", "DD_18", "DD18",
  "DD5", "eFFP", "EMT", "EXT", "FFP", "MAP", "MAT", "MCMT",
  "MSP", "MWMT", "NFFD", "PPT_sm", "PPT_wt", "SHM", "YEAR")

#4. Select visits for the region of interest----
ss <-  dd$bcr4 > 0

#5. Write function to run bootstrapped models----

#' Run 1 iteration for a species
#'
#' @param b bootstrap id (positive integer)
#' @param spp Species code
#' @param cn colnames from dd to use in BRT
#' @param ss subset of rows to use (region, BCR etc.)
#' @param nt number of trees
#' @param verbose print stuff (silent when running with Rscript)

run_brt_boot1 <- function(b, spp, cn, ss, nt=10^4, verbose=interactive()) {
  
    t0 <- proc.time()["elapsed"]
    
    #create data subset for BCR unit
    DAT <- data.frame(
        count=as.numeric(yy[ss, spp]),
        offset=off[ss, spp],
        cyid=dd$cyid[ss],
        dd[ss, cn])
    
    # subsample based on 2.5x2.5km^2 cell x year units
    if (b > 1)
        DAT <- DAT[sample.int(nrow(DAT)),]
    DAT <- DAT[!duplicated(DAT$cyid),]
    if (b > 1)
        DAT <- DAT[sample.int(nrow(DAT), replace=TRUE),]
    
    # 0 detection output
    if (sum(DAT$count) < 1) {
        out <- structure(
            sprintf("0 detections for %s", spp),
            class="try-error")
    } else {
        if (verbose)
            cat("\nFitting gbm::gbm for", spp, "in", ", b =", b, ", n =", nrow(DAT), "\n")
        out <- try(gbm::gbm(DAT$count ~ . + offset(DAT$offset),
                            data=DAT[,cn],
                            n.trees = nt,
                            interaction.depth = 3,
                            shrinkage = 1/nt,
                            bag.fraction = 0.5,
                            distribution = "poisson",
                            var.monotone = NULL,
                            keep.data = FALSE,
                            verbose = verbose,
                            n.cores = 1))
    }
    # package output
    attr(out, "__settings__") <- list(
        b=b, spp=spp, cn=cn, ss=ss, elapsed=proc.time()["elapsed"]-t0)
    out
}

#6. Run model for selected species and bootstrap number----

# Run this function using lapply to run across multiple bootstraps/regions/species
# You can set verbose=FALSE to hide the output during model fitting

o <- run_brt_boot1(b=1, spp="AMRO", cn=cn, ss=ss)

#7. Save output----

save(o, file="ModelRun.Rdata")