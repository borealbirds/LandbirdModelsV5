# ---
# title: National Models 4.1 - presentation statistics & figures
# author: Elly Knight
# created: November 2, 2023
# ---

#NOTES################################

#This script calculates summary statistics about the national models and creates figures for presentations

#This script was initially written for a BAM update meeting on November 2, 2023 and compares the 4.0 version to the 4.1 version of the models.

#PREAMBLE############################

#1. Load packages----

library(tidyverse) #basic data wrangling
library(sf)
library(gridExtra)

#2. Set root path for data on google drive----
root <- "G:/Shared drives/BAM_NationalModels/NationalModels4.1"

#3. Load data----
e0 <- new.env()
load(file="G:/Shared drives/BAM_NationalModels/NationalModels4.0/data/BAMdb-GNMsubset-2020-01-08.Rdata", env=e0)

e1 <- new.env()
load(file=file.path(root, "Data", "03_NM4.1_data_stratify.R"), env=e1)

#4. Load BCRs----
unit1 <- read_sf(file.path(root, "Regions", "BAM_BCR_NationalModel.shp")) %>% 
  mutate(bcr=paste0("bcr", subUnit))

unit0 <- read_sf("G:/Shared drives/BAM_NationalModelsV(EMB)/NationalModelsV4.0/BCRunits2.shp") %>% dplyr::filter(!is.na(STATES)) %>% 
  st_transform(crs=terra::crs(bcr1))

#COMPARE 4.0 vs 4.1####################

#1. Sample size-----
nrow(e0$dd)
nrow(e1$visit.bcr)

#2. Temporal distribution----
visit.n <- e0$dd %>% 
  group_by(YEAR) %>% 
  summarize(n=n()) %>% 
  ungroup() %>% 
  rename(year = YEAR) %>% 
  mutate(version = "V4.0") %>% 
  rbind(e1$visit.bcr %>% 
          group_by(year) %>% 
          summarize(n=n()) %>% 
          ungroup() %>% 
          mutate(version = "V4.1")) %>% 
  dplyr::filter(year > 1993) %>% 
  pivot_wider(names_from=version, values_from=n, values_fill=0) %>% 
  pivot_longer(V4.0:V4.1, names_to="version", values_to="n")

plot.n <- ggplot(visit.n) +
  geom_col(aes(x=year, y=n, fill=version), position = position_dodge()) +
  scale_fill_manual(values=c("grey30", "grey70")) +
  theme(legend.position = "bottom")

ggsave(plot.n, filename=file.path(root, "Figures", "SampleSizeComparison.jpeg"), width=6, height = 4)

#3. Spatial coverage----
loc.0 <- e0$dd %>% 
  dplyr::select(X, Y) %>% 
  unique() %>% 
  mutate(version="V4.0")

loc.1 <- e1$visit.bcr %>% 
  rename(X=lon, Y=lat) %>% 
  dplyr::select(X, Y) %>% 
  unique() %>% 
  mutate(version="V4.1")

map <- map_data("world", region=c("Canada", "USA"))

plot.loc0 <- ggplot(loc.0) +
  geom_polygon(data=map, aes(x=long, y=lat, group=group), fill="grey80", colour="gray90", linewidth=0.3) +
  geom_hex(aes(x=X, y=Y), bins=200) +
  scale_fill_viridis_c() +
  coord_sf(xlim=c(-170, -50), ylim=c(38, 75), expand=FALSE, crs=4326) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "right") +
  ggtitle("V4.0")
plot.loc0

plot.loc1 <- ggplot(loc.1) +
  geom_polygon(data=map, aes(x=long, y=lat, group=group), fill="grey80", colour="gray90", linewidth=0.3) +
  geom_hex(aes(x=X, y=Y), bins=200) +
  scale_fill_viridis_c() +
  coord_sf(xlim=c(-170, -50), ylim=c(38, 75), expand=FALSE, crs=4326) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "right") +
  ggtitle("V4.1")
plot.loc1

ggsave(grid.arrange(plot.loc0, plot.loc1, nrow=2),
       filename=file.path(root, "Figures", "DataCoverage.jpeg"), width=6, height=7)

#4. Study area----
sa1 <- unit1 %>% 
  summarize() %>% 
  mutate(version = "4.1")

sa0 <- unit0 %>% 
  summarize() %>% 
  mutate(version = "4.0")

sa <- rbind(sa1, sa0)

plot.sa <- ggplot() +
  geom_sf(data=sa, aes(fill=version)) +
  scale_fill_manual(values=c("grey30", "grey70")) +
  theme(legend.position = "bottom")
plot.sa

ggsave(plot.sa, filename=file.path(root, "Figures", "StudyAreaComparison.jpeg"), width=6, height = 4)

#5. BCRs----
bcr1 <- unit1 %>% 
  mutate(bcr = factor(subUnit)) %>% 
  group_by(bcr) %>% 
  summarize() %>% 
  ungroup()

bcr0 <- unit0 %>% 
  mutate(bcr = factor(BCR)) %>% 
  group_by(bcr) %>% 
  summarize() %>% 
  ungroup()

plot.bcr <- ggplot() +
  geom_sf(data=bcr0, aes(fill=bcr)) +
  geom_sf(data=bcr1, fill="grey30", alpha=0.3, colour="black", linewidth=0.7) +
  scale_fill_discrete(name="BCR version 4.0") +
  theme(legend.position = "bottom")

ggsave(plot.bcr, filename=file.path(root, "Figures", "BCRComparison.jpeg"), width=6, height = 5)
