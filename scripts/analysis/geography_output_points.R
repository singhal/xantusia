rm(list = ls())
library(sf)
library(dplyr)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Xantusia/")

source("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_speciation/scripts/geography_helper.R")

# nuclear data
d = read.csv("metadata/xantusia_samples_v9.csv")
d1 = d[complete.cases(d$OTU), ]
d1 = d1[complete.cases(d1$latitude), ]

xx = d1 %>% dplyr::select(SAMPLE = sample_fix, LAT = latitude, LON = longitude, OTU)

# get rid of some outliers
xx = xx[!xx$SAMPLE %in% c("LACM136734"), ]

world_map <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# make comparision set
taxa = unique(xx$OTU)
combo = combn(taxa, 2)
res = vector("list", ncol(combo))

pts = vector("list", length(taxa))
for (i in 1:length(taxa)) {
  # make list of points
  pts[[i]] = xx %>% dplyr::filter(OTU == taxa[i]) %>% 
    dplyr::select(LON, LAT) %>% 
    st_as_sf(coords = c("LON", "LAT")) %>% 
    st_set_crs(sf::st_crs(world_map))
}
names(pts) = taxa


for (i in 1:ncol(combo)) {
  t1 = combo[1, i]
  t2 = combo[2, i]
  if (i %% 100 == 0) {
    message("doing combo ", i, " of ", ncol(combo), "\n")
  }
  
  # kilometer
  mean_distance = as.numeric(mean(st_distance(pts[[t1]], pts[[t2]]))) / 1000
  # smallest distance
  min_distance = as.numeric(min(st_distance(pts[[t1]], pts[[t2]]))) / 1000
  median_distance = as.numeric(median(st_distance(pts[[t1]], pts[[t2]]))) / 1000
  
  res[[i]] = c(t1, t2, mean_distance, min_distance, median_distance)
}

res2 = data.frame(do.call("rbind", res))
names(res2) =  c("species1", "species2", "mean_pt_distance", "min_pt_distance", "median_pt_distance")
res3 = res2[complete.cases(res2$species1), ]
nums = c("mean_pt_distance", "min_pt_distance", "median_pt_distance")
for (num in nums) {
  res3[, num] = as.numeric(res3[,num])
}
write.csv(res3, paste0("data/geography/point_distances.csv"))
