rm(list = ls())

library(ape)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggtree)
library(patchwork)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Xantusia/")
source("./scripts/gps_colors.R")

# read in & root tree
t = read.tree("data/mtDNA/treeinference/Xantusia.treefile")
t = root(t, c("EU116520.1", "EU116521.1"), resolve.root = T)

d = read.csv("data/mtDNA/mtDNA_v7.csv", na.strings = c("NA", ""))
d = d[complete.cases(d$OTU), ]
d = d[complete.cases(d$museum), ]

# check that all tips in the tree are in our dataset
t$tip.label = gsub("\\.1", "", t$tip.label)
# no, some are missing
# but not worth it to track down
t1 = drop.tip(t, t$tip.label[ ! t$tip.label %in% d$genbank ])

# get rid of duplicates
tips = gsub(" ", "", d[match(t1$tip.label, d$genbank), "museum"])
t2 = drop.tip(t1, t1$tip.label[duplicated(tips)])

t2 = read.tree(text = write.tree(ladderize(t2)))

# tx = t2
# tx$tip.label = paste(d[match(tx$tip.label, d$genbank), "OTU"], 
#                       d[match(tx$tip.label, d$genbank), "museum"])
# otus1 = unique(d$OTU)
# for (otu in otus1) {
#   pdf(paste0("~/Desktop/mtdna.", otu, ".pdf"), width = 10, height = 50)
#   plot.phylo(tx, show.tip.label = F)
#   otulist = d[match(t2$tip.label, d$genbank), "OTU"]
#   for (i in 1:length(tx$tip.label)) {
#     if (otulist[i] == otu) {
#       ape::tiplabels(tx$tip.label[i], i, frame = "rec", adj = c(0, 0.5), cex = 0.5, bg = "yellow")
#     } else {
#       ape::tiplabels(tx$tip.label[i], i, frame = "none", adj = c(0, 0.5), cex = 0.2)
#     }
#   }
#   dev.off()
# }

d %>% group_by(OTU, species) %>% summarize(n = n()) %>% View()

t3 = t2
t3$tip.label = paste0(d[match(t3$tip.label, d$genbank), "OTU"], " ",
                      d[match(t3$tip.label, d$genbank), "museum"])
write.tree(t3, "~/Desktop/mtDNA_tree.tre")
t2$tip.label = d[match(t2$tip.label, d$genbank), "museum"]


# otus in tree
otus1 = unique(d$OTU)
# what are the new OTUs?
otus1[! (otus1 %in% x$OTU) ]

otu_cols1 = otu_cols[otus1]
names(otu_cols1) = otus1
otu_cols1[!complete.cases(otu_cols1)] = "gray80"
# cols1 = c("black", otu_cols1)
# names(cols1)[1] = "0"
# treeplot = ggtree(t2)
# for (i in 1:length(otus1)) {
#   sp = otus1[i]
#   tmpinds = t2$tip.label[ d[match(t2$tip.label, d$museum), "OTU"] == sp]
#   if (length(tmpinds) > 1) {
#     node = phytools::findMRCA(t2, tmpinds)
#   } else {
#     node = match(tmpinds, t2$tip.label)
#   }
#   
#   if (sp == "bolsonae") {
#     tmpinds = c("LACM138478", "MZFC 4771")
#     node = phytools::findMRCA(t2, tmpinds)
#   }
#   if (sp == "vigilis_YV") {
#     tmpinds = c("DHL-JT6", "CAS224956")
#     node = phytools::findMRCA(t2, tmpinds)
#   } 
#   
#   ypos = (max(which(t2$tip.label %in% tmpinds)) + min(which(t2$tip.label %in% tmpinds))) / 2
#   ypos = Ntip(t2) - ypos + 1
#   
#   # treeplot = treeplot + geom_hilight(node = node, fill = otu_cols1[i], 
#   #                                   type = "gradient", gradient.direction = "tr", 
#   #                                   alpha = 0.8) +
#   treeplot = collapse(treeplot, node, 'min', fill= otu_cols1[i]) +
#     geom_text(x = 0.3, y = ypos, label=sp)
#   }
# treeplot = treeplot + theme(legend.position = "none",
#                             plot.margin = unit(c(0,0,0,0), "cm"))
# cowplot::save_plot("~/Desktop/mtdna.pdf", treeplot, base_height = 25, 
#                    base_width = 10)

# otu order
otuorder = unique(d[match(t2$tip.label, d$museum), "OTU"])


# read in all mapping data
world = ne_countries(scale = "medium", returnclass = "sf")

# get nuclear dta
x = read.csv("metadata/xantusia_samples_v10.csv")
dx = rbind(d %>% dplyr::select(sample, LON, LAT, OTU) %>% mutate(type = "mt"),
           x %>% dplyr::select(sample = sample_fix, LON = longitude,
                        LAT = latitude, OTU) %>% mutate(type = "nuc")) %>%
  filter(complete.cases(OTU))
dx$OTU2 = factor(dx$OTU, levels = otuorder)

mapplot = ggplot() +
  geom_sf(data = world, fill = "white", color = "gray70", size = 0.15) +
  geom_point(data = dx, aes(x = LON, y = LAT, fill = OTU, shape = type), 
             alpha = 0.7, size = 2, color = "black", stroke = 0.2) +
  facet_wrap(~OTU2, ncol = 4) +
  xlim(-121, -102) + ylim(21, 38) +
  theme_void() +
  guides(fill = "none", size = "none", alpha = "none", shape = "none") +
  scale_fill_manual(values = otu_cols1) +  
  scale_shape_manual(values = c(21, 24)) +
  theme(panel.background = element_rect(fill = "aliceblue"),
        panel.border = element_rect(fill = NA, color = "gray30", size = 0.5),
        strip.text = element_text(size = 9, face = "italic", 
                                  margin = margin(t = 0, r = 0, b = 2, l = 0)), 
        strip.background = element_blank(),
        plot.margin = margin(r = 2, l = 2, t = 0, b = 0))
cowplot::save_plot("~/Desktop/maps.png", mapplot, base_height = 10, base_width = 6)
