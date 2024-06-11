rm(list = ls()) 
library(rnaturalearth)
library(rnaturalearthdata)
library(scatterpie)
library(sf)
library(tidyverse)
library(cowplot)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Xantusia/")
world <- ne_countries(scale = "medium", returnclass = "sf")

# get colors
source("./scripts/gps_colors.R")

# get metadata
sp = read.csv("metadata/xantusia_samples_v10.csv")

# get in the best data
cl = read.csv("data/admixture/best_runs.csv")
res = vector("list", nrow(cl))
for (i in 1:nrow(cl)) {
  run = cl[i, "run"]
  k = cl[i, "k"]
  group = cl[i, "group"]
  a = read.table(paste0("data/admixture/", group, ".", run, ".miss0.7.MAC2.thinned.", k, ".Q"), sep = " ")
  n = read.table(paste0("data/input_files/", group, ".0.miss0.7.MAC2.thinned.fam"))
  
  a$ind = n$V2
  a = cbind(a, sp[match(a$ind, sp$sample_fix), c("OTU", "cluster", "latitude", "longitude")])

  # remove the likely mis-identified ind
  a = a[a$ind != "LACM136734",]
  
  # summarize data a bit
  a1 = a %>% mutate(cluster=replace(cluster, cluster=="magdalena", "sherbrookei")) %>%
    filter(complete.cases(latitude)) %>%
    dplyr::select(-ind, -OTU) %>% 
    mutate(latitude2 = round(latitude, 2), longitude2 = round(longitude, 2)) %>%
    group_by(cluster, latitude2, longitude2) %>% summarize_all(mean) %>%
    ungroup()
  
  # need to fix the colors
  map1 = a %>% mutate(cluster=replace(cluster, cluster=="magdalena", "sherbrookei")) %>%
    dplyr::select(-latitude, -longitude) %>% 
    gather("pop", "prob", -OTU, -ind, -cluster) %>% 
    group_by(ind) %>% slice_max(prob) %>%
    group_by(cluster, OTU, pop) %>% summarize(n = n())
  map2 = map1 %>% group_by(cluster) %>% slice_max(n)
  cols = cl_cols[map2$cluster]
  map3 = map1 %>% group_by(cluster) %>%
    summarize(OTU = paste(unique(OTU), collapse = ":"))
  names(cols) = map3[match(names(cols), map3$cluster), "OTU"] %>% pull()
  map2$OTU = pull(map3[match(map2$cluster, map3$cluster), "OTU"])
  
  # need to fix the range
  ylim1 = min(a1$latitude2) - abs(max(a1$latitude2) -  min(a1$latitude2)) * 0.1
  ylim2 = max(a1$latitude2) + abs(max(a1$latitude2) -  min(a1$latitude2)) * 0.1
  
  xlim1 = min(a1$longitude2) - abs(max(a1$longitude2) -  min(a1$longitude2)) * 0.1
  xlim2 = max(a1$longitude2) + abs(max(a1$longitude2) -  min(a1$longitude2)) * 0.1
  
  # rename rows
  for (xx in 1:nrow(map2)) {
    names(a1)[ match(map2[xx, "pop"], names(a1)) ] = map2[xx, "OTU"]
  }
  
  sizes = c(3, 0.7, 1)

  map = ggplot() +
    ## Adding baseline map:
    geom_sf(data = world, fill = "white", color = "gray70", size = 0.15) +
    ## Pie chart of the ancestry coefficients:
    geom_scatterpie(data = a1, aes(x = longitude2, y = latitude2),
                    pie_scale = sizes[i], alpha = 1,
                    size = 0.3, 
                    cols = map3$OTU) +
    # zoom into where it matters
    xlim(xlim1, xlim2) + ylim(ylim1, ylim2) +
    xlab("longitude") + ylab("latitude") +
    theme_classic() +
    scale_fill_manual(values = cols) +
    ## Setting some other elements:
    theme(panel.background = element_rect(fill = "aliceblue"),
          panel.border = element_rect(fill = NA, color = "gray30", linewidth = 0.5),
          legend.title = element_blank())
  
  save_plot(paste0("figures/admixture_map.", group, ".png"), map)
  
}

