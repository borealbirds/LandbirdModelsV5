library(mefa4)
library(gbm)

## RData file to load
fn <- "E:/MelinaStuff/BAM/NationalModelv4.1/runNationalModel/data/NatMod41-2021-12-10.RData"
load(fn)

## some checks
ls()
dim(dd)
dim(yy)
dim(off)
keep <- intersect(dd$PKEY, rownames(yy))
dd <- dd[dd$PKEY %in% keep,]
dd <- nonDuplicated(dd, PKEY, TRUE)
yy <- yy[keep,]
off <- off[keep,]
stopifnot(all(colnames(off) == colnames(yy)))
stopifnot(all(rownames(off) == rownames(yy)))
stopifnot(all(rownames(dd) == rownames(yy)))

## variables to be used as predictors (you can vary this by region)
## Note: you have to include ARU, YEAR, ROAD etc
## I only take care of spp count and offset (plus the cyid)
cn <- c("TPI", "TRI", "slope", "roughness",
  "RoadDist", "lf", "AHM", "bFFP", "CMD", "DD_0", "DD_18", "DD18",
  "DD5", "eFFP", "EMT", "Eref", "EXT", "FFP", "MAP", "MAT", "MCMT",
  "MSP", "MWMT", "NFFD", "PAS", "PPT_sm", "PPT_wt", "RH", "SHM")

## subset to define region (logical or numeric)
ss <-  dd$bcr4 > 0
#ss <- which(ss)[1:1000]

#' Function to run 1 iteration for a species
#'
#' @param b bootstrap id (positive integer)
#' @param spp Species code
#' @param cn colnames from dd to use in BRT
#' @param ss subset of rows to use (region, BCR etc.)
#' @param nt number of trees
#' @param verbose print stuff (silent when running with Rscript)
run_brt_boot1 <- function(b, spp, cn, ss, nt=10^4, verbose=interactive()) {
    t0 <- proc.time()["elapsed"]
    
    ## create data subset for BCR unit
    DAT <- data.frame(
        count=as.numeric(yy[ss, spp]),
        offset=off[ss, spp],
        cyid=dd$cyid[ss],
        dd[ss, cn])
    ## subsample based on 2.5x2.5km^2 cell x year units
    if (b > 1)
        DAT <- DAT[sample.int(nrow(DAT)),]
    DAT <- DAT[!duplicated(DAT$cyid),]
    if (b > 1)
        DAT <- DAT[sample.int(nrow(DAT), replace=TRUE),]
    ## 0 detection output
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
    attr(out, "__settings__") <- list(
        b=b, spp=spp, cn=cn, ss=ss, elapsed=proc.time()["elapsed"]-t0)
    out
}


o <- run_brt_boot1(b=1, spp="AMRO", cn=cn, ss=ss)
o <- run_brt_boot1(b=2, spp="AMRO", cn=cn, ss=ss, verbose=FALSE)


# You can decide how to loop over regions/species/iterations
# but save often (every b/spp/region combination)

fo <- file.path("E:/MelinaStuff/BAM/NationalModelv4.1/runNationalModel/outBRT", paste0("brt-", spp, "-", bcr, "-", b, ".RData"))
save(o, file=fo)

#Explore
load("E:/MelinaStuff/BAM/NationalModelv4.1/runNationalModel/outBRT/brt-AMRO-bcr4-10.RData", AMRObcr4 <- new.env())
ls(AMRObcr4)
AMRObcr4$o
summary(AMRObcr4$o)

library(gbm)
plot.gbm(AMRObcr4$o)

