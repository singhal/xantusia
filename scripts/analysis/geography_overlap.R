rm(list = ls())
library(sf)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Xantusia/data/geography/")

rr = list.files("alphahull.ranges/", full.names = T,
           pattern = ".shp")
rr1 = lapply(rr, st_read)

allres = vector("list", length(rr1))
world_map <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

for (i in 1:length(rr1)) {
  allres[[i]] = rr1[[i]] %>% 
      st_simplify(preserveTopology = TRUE, dTolerance = 0.1) %>%
      st_buffer(dist=0) %>%
      st_set_crs(sf::st_crs(world_map))
}
  
taxa = gsub(".shp", "", gsub(".*//", "", rr))
names(allres) = taxa
combo = combn(taxa, 2)
res = vector("list", ncol(combo))
  
for (i in 1:ncol(combo)) {
  t1 = combo[1, i]
  t2 = combo[2, i]
  if (i %% 100 == 0) {
    message("doing combo ", i, " of ", ncol(combo), "\n")
  }
    
  if (t1 %in% names(allres)) {
    r1 = allres[[t1]]
    if (length(as.numeric(st_area(r1))) == 0) {
      r1 = NULL
    }
  } else {
    r1 = NULL
  }
  if (t2 %in% names(allres)) {
    r2 = allres[[t2]]
    if (length(as.numeric(st_area(r2))) == 0) {
      r2 = NULL
    }
  } else {
    r2 = NULL
  }
    
  overlap = NA
  distance = NA
  # both are ranges
  sf_use_s2(FALSE)
  if (is.null(r1) | is.null(r2)) {
    res[[i]] = c(t1, t2, NA, NA)
  } else {
    min_area = min(c(sum(sf::st_area(r1)), sum(sf::st_area(r2))))
    intersect = suppressMessages(sf::st_intersection(r1, r2))
    overlap = as.numeric(sum(sf::st_area(intersect)) / min_area)
    if (overlap > 0) {
      distance = 0
    } else {
      distance = min(sf::st_distance(r1, r2))
    }
    res[[i]] = c(t1, t2, overlap, distance)
  }
}

res2 = data.frame(do.call("rbind", res))
names(res2) =  c("species1", "species2", "overlap", "geo_dist")
res3 = res2[complete.cases(res2$species1), ]
nums = c("overlap", "geo_dist")
for (num in nums) {
  res3[, num] = as.numeric(res3[,num])
}

source("../../scripts/gps_colors.R")

t = ape::read.nexus("data/phylogeny/mcc.tre")
dd = as.data.frame(ape::cophenetic.phylo(t)) 
dd$species1 = rownames(dd)
dd1 = dd %>% tidyr::gather(species2, phy_dist, -species1) %>%
  mutate(phy_dist = phy_dist / 2)
res4 = res3 %>% left_join(dd1) %>% left_join(combo)
x = read.csv("data/geography/point_distances.csv", row.names = 1)
x2 = x %>% dplyr::select(taxa1 = species2, taxa2 = species1, mean_pt_distance,
             min_pt_distance, median_pt_distance)
names(x2) = names(x)
xx = rbind(x, x2)
res5 = res4 %>% left_join(xx) %>%
  filter(complete.cases(phy_dist)) %>%
  mutate(overlapT = ifelse(overlap > 0.1, TRUE, FALSE))

# add group identity
get_group <- function(y) {
  return(names(which(unlist(lapply(otus, function(x) y %in% x)) == TRUE)))
}
res5$group1 = unlist(sapply(res5$species1, get_group))
res5$group2 = unlist(sapply(res5$species2, get_group))

res5 = res5[res5$group1 == res5$group2, ]

res5 %>% filter(complete.cases(phy_dist)) %>% nrow()
res5 %>% filter(complete.cases(phy_dist)) %>% filter(overlap > 0.2) %>% 
  dplyr::select(species1, species2, overlap, median_pt_distance, min_pt_distance) %>%
  View()

 
a = ggplot(res5, aes(phy_dist, geo_dist / 1000, fill = type)) + 
   geom_jitter(alpha = 0.8, size = 2, shape = 21) + 
   xlab("phylogenetic distance") + ylab("geographic distance (km)") +
   theme_classic() +
   scale_fill_manual(values = habcols)
b = ggplot(res5, aes(phy_dist, overlap, fill = type)) + 
  geom_jitter(alpha = 0.8, size = 2, shape = 21, width = 0.001, height = 0.01) + 
  xlab("phylogenetic distance") + ylab("range overlap") +
  theme_classic() + theme(legend.title = element_blank()) +
  scale_fill_manual(values = habcols) 

prow <- plot_grid(
  a + theme(legend.position="none"),
  b + theme(legend.position="none"),
  align = 'vh',
  labels = c("A", "B"),
  hjust = -1,
  nrow = 1
)
prow

legend <- get_legend(
  # create some space to the left of the legend
  a + theme(legend.box.margin = margin(0, 0, 0, 0))
)

# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
geoplt = plot_grid(prow, legend, rel_widths = c(3, 0.5))

ab = cowplot::save_plot("~/Dropbox (Personal)/research_projects-ACTIVE/Xantusia/figures/overlap.png", 
                        geoplt, base_width = 8, base_height = 3)
