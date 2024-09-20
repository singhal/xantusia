rm(list = ls())
library(sf)
library(dplyr)
library(alphahull)
library(hull2spatial)
library(ggplot2)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Xantusia/")

source("~/Dropbox (Personal)/research_projects-ACTIVE/Sphenomorphine_speciation/scripts/geography_helper.R")

# nuclear data
d = read.csv("metadata/xantusia_samples_v10.csv")
d1 = d[complete.cases(d$OTU), ]
d1 = d1[complete.cases(d1$latitude), ]

xx = d1 %>% dplyr::select(SAMPLE = sample_fix, LAT = latitude, LON = longitude, OTU)
# get rid of outlier - likely wrong georef
xx = xx[!xx$SAMPLE %in% c("LACM136734"), ]

world_map <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# species
taxa = unique(xx$OTU)

# make alphahull
min_pts = 3

outdir = "data/geography/alphahull.ranges"
if (!dir.exists(outdir)) {
  dir.create(outdir)
}


for (i in 1:length(taxa)) {
  sp = taxa[i]
  
  # combine with our data occurrences
  pts2 = xx[which(xx$OTU == sp), c("LON", "LAT")]
  pts3 = pts2 %>% filter(complete.cases(LON)) %>%
    mutate_all(round, digits = 6) %>% distinct()
  
  if (nrow(pts3) >= min_pts) {
    pts4 = st_as_sf(pts3, coords = c("LON", "LAT")) %>% 
      st_set_crs(sf::st_crs(world_map))
    
    # make shape files
    alphaval = 10
    all_r4 = get_alphahull(pts3, pts4, alphaval)
    # clip to the world
    all_r5 = st_intersection(all_r4, world_map)
    
    st_write(all_r5, paste0(outdir, "/", sp, ".shp"), append = F)
    
  } else if (nrow(pts3) > 0) {
    pts4 = st_as_sf(pts3, coords = c("LON", "LAT")) %>% 
      st_set_crs(sf::st_crs(world_map))
    all_r4 = st_union(st_buffer(pts4, dist = 0.1)) %>% 
      st_set_crs(sf::st_crs(world_map))
    # clip to the world
    all_r5 = st_intersection(all_r4, world_map)
    
    st_write(all_r5, paste0(outdir, "/", sp, ".shp"), append = F)
  }
  
  # plot it
  outfile = paste0(outdir, "/", sp, ".pdf")
  mapplt = ggplot() +
    geom_sf(data = world_map, fill = "transparent", color = "gray30") +
    geom_sf(data = all_r5, fill = "forestgreen", alpha = 0.3,
            color = "transparent") +
    geom_sf(data = pts4, shape = 21, fill = "white", size = 1) +
    coord_sf(xlim = c(-122, -101), ylim = c(20, 38)) +
    theme_void() + labs(title = sp)
  cowplot::save_plot(outfile, mapplt)
}
