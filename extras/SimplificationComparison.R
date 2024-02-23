#1. Load packages----
library(tidyverse) #basic data wrangling
library(usdm) #vif

#2. Set root path for data on google drive----
root <- "G:/Shared drives/BAM_NationalModels/NationalModels5.0"

#3. Load data packages with offsets and covariates----
load(file.path(root, "Data", "04_NM5.0_data_stratify.R"))

#4. Turn off scientific notation---
options(scipen=99999)

#5. Set thresholds----
rmax <- 0.9
vifmax <- 10

#6. Figure out which species & BCR combos to look at----
spbcr <- data.frame(file = list.files("output/simplifying"),
                    path = list.files("output/simplifying", full.names = TRUE)) %>% 
  dplyr::filter(str_sub(file, -9, -5)=="local") %>% 
  separate(file, into=c("f1", "spp", "bcr", "f2"), sep="_", remove=FALSE) %>% 
  dplyr::select(-f1, -f2)

#7. Set up loop----
droplist <- data.frame()
r.i <- data.frame()
for(i in 1:nrow(spbcr)){
  
  #8. Get visits to include----
  set.seed(i)
  visit.i <- gridlist[bcrlist[,spbcr$bcr[i]],] %>% 
    group_by(year, cell) %>% 
    mutate(rowid = row_number(),
           use = sample(1:max(rowid), 1)) %>% 
    ungroup() %>% 
    dplyr::filter(rowid==use)
  
  #9. Get covariates----
  covlist.i <- covlist %>% 
    dplyr::filter(bcr==spbcr$bcr[i]) %>% 
    pivot_longer(ERAMAP_1km:mTPI_1km, names_to="cov", values_to="use") %>% 
    dplyr::filter(use==TRUE)
  cov.i <- cov[cov$id %in% visit.i$id, colnames(cov) %in% covlist.i$cov] %>% 
    select_if(is.numeric) %>% 
    data.frame()
  
  #10. Run vif----
  vif.i <- vifstep(cov.i, vifmax)
  
  #11. Get pairwise r and summarize total r per var----
  r.i <- cor(cov.i, use="pairwise.complete.obs") %>% 
  # r.i <- vif.i@corMatrix %>% 
    data.frame() %>% 
    mutate(var1 = row.names(.)) %>% 
    pivot_longer(colnames(.)[1]:colnames(.)[nrow(.)], names_to="var2", values_to="r") %>% 
    mutate(rabs = abs(r)) %>% 
    dplyr::filter(var1!=var2) %>% 
    arrange(-rabs) %>% 
    group_by(var1) %>% 
    mutate(totalr1 = sum(rabs, na.rm = TRUE)) %>% 
    ungroup() %>% 
    mutate(spp = spbcr$spp[i],
           bcr = spbcr$bcr[i]) %>% 
    rbind(r.i)
  
  #12. Take out by r----
  r.out <- r.i %>% 
    dplyr::filter(rabs > rmax) %>% 
    group_by(rabs) %>% 
    dplyr::filter(totalr1 == max(totalr1)) %>% 
    ungroup() %>% 
    dplyr::select(var1)
  
  #13. Read in the gbm.simplify results----
  simp.i <- read.csv(spbcr$path[i]) %>% 
    dplyr::filter(!is.na(order))
  
  #14. Put it all together----
  droplist <- covlist.i %>% 
    dplyr::select(-use) %>% 
    mutate(gbm = ifelse(cov %in% simp.i$preds, 0, 1),
           r = ifelse(cov %in% r.out$var1, 0, 1),
           vif = ifelse(cov %in% vif.i@excluded, 0, 1),
           spp = spbcr$spp[i]) %>% 
    rbind(droplist)
    
  print(paste0("Finished loop ", i))
  
}

#15. Interrogate the ones that gbm.simplify kept that vif dropped----
hmm <- droplist %>% 
  dplyr::filter(vif==0 & gbm==1) %>% 
  left_join(r.i %>% 
              rename(cov = var1, 
                     rval = r),
            multiple = "all") %>% 
  group_by(cov) %>% 
  dplyr::filter(rabs==max(rabs)) %>% 
  ungroup() 