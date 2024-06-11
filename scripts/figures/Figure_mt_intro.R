rm(list = ls())

library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Xantusia/")

source("./scripts/gps_colors.R")

x = read.csv("./metadata/xantusia_samples_v8.csv")
d = read.csv("./data/mtDNA/mtDNA_v6.csv")
world = ne_countries(scale = "medium", returnclass = "sf")

x$mtDNA_OTU = d[match(x$sample_fix, gsub(" ", "", d$museum)), "OTU"]

x1 = x[complete.cases(x$mtDNA_OTU), ]

x1 %>% filter(complete.cases(OTU)) %>% 
  group_by(OTU, mtDNA_OTU) %>% summarize(n = n()) %>%
  filter(OTU != mtDNA_OTU) 

otus = c("vigilis_EM", "vigilis_KCAV", "vigilis_OV")
xx = x %>% filter(mtDNA_OTU %in% otus | OTU %in% otus)  %>%
  filter(complete.cases(OTU)) %>%
  # remove outlier individual taht i think might be misid
  filter(sample_fix != "LACM136734")

# add arrows for mismatch
xx1 = xx[which(xx$mtDNA_OTU != xx$OTU), ]
xx1$xend = xx1$longitude + 0.2
xx1$xstart = xx1$longitude + 0.45
xx1$yend = xx1$latitude 
xx1$ystart = xx1$latitude 

xx = ggplot() +
  geom_sf(data = world, fill = "white", color = "gray70", size = 0.15) +
  geom_point(data = xx, 
             aes(x = longitude, y = latitude, col = OTU), 
             size = 2.5, shape = 16) +
  geom_point(data = xx %>% filter(complete.cases(mtDNA_OTU)), 
             aes(x = longitude, y = latitude, col = mtDNA_OTU), 
             size = 3, shape = 2, stroke = 2) +
  geom_segment(data = xx1, aes(x = xstart, y = ystart, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.2, "cm"))) +
  xlim(-120, -110) + ylim(32.8, 38) +
  theme_bw() +
  scale_colour_manual(values = otu_cols) +
  theme(panel.background = element_rect(fill = "aliceblue"),
        panel.border = element_rect(fill = NA, color = "gray30", size = 0.5))

save_plot("figures/mtDNA_introgression.png", xx, base_height = 5, 
          base_width = 8)
