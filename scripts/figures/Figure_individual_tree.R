rm(list = ls())

library(ape)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggtree)
library(patchwork)
library(rnaturalearth)
library(rnaturalearthdata)

### Setting up the Key Files ----
setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Xantusia/")
source("./scripts/gps_colors.R")

x = read.csv("./metadata/xantusia_samples_v10.csv")

### Preparing Tree ----

plot_tree <- function() {
  # read in tree
  t = read.tree("data/phylogeny/xantusia.low_inds_removed.treefile")
  outs = x[x$genus == "Lepidophyma", "sample_fix"]
  outs = outs[outs %in% t$tip.label]

  t1 = root(t, outs, resolve.root = TRUE)
  t2 = read.tree(text = write.tree(ladderize(keep.tip(t1, inds))))
  
  sh = as.numeric(gsub("/\\S+", "", t2$node.label))
  nodes =  seq(Ntip(t2) + 1, Ntip(t2) + Nnode(t2))
  names(sh) = nodes
  desc = lapply(nodes, function(x) t2$tip.label[ phangorn::Descendants(t2, x, type = "tips")[[1]]])
  num_otu = lapply(desc, function(y) length(unique(x[x$sample_fix %in% y, "OTU"])))
  keep = sh[num_otu > 1]
  keep = keep[keep > 95]
  
  # color branches by clade
  gps = lapply(otus1, function(sp) x[which(x$OTU == sp), "sample_fix"])
  names(gps) = otus1
  t3 = groupOTU(t2, gps)

  cols1 = c("black", otu_cols[otus1])
  names(cols1)[1] = "0"
  treeplot = ggtree(t3, size = 0.75, ladderize = FALSE, aes(color = group)) +
    scale_color_manual(values = cols1) + 
    geom_nodepoint(aes(subset = node %in% as.numeric(names(keep))), 
                   size=1, fill='white', shape = 21)
  for (i in 1:length(otus1)) {
    sp = otus1[i]
    tmpinds = t3$tip.label[ x[match(t3$tip.label, x$sample_fix), "OTU"] == sp]
    treeplot = treeplot + 
      geom_strip(taxa1 = tmpinds[1], taxa2 = tmpinds[length(tmpinds)], 
                 label = sp, offset = 0.0001, face = "italic",
                 offset.text = 0.0003, color = otu_cols[sp])}
  treeplot = treeplot + theme(legend.position = "none",
    plot.margin = unit(c(0,0,0,0), "cm")) + 
    xlim(0, 0.028)
  return(list(t3, treeplot))
}

### Admixture plotting ----

plot_admixture = function(t3) {
  # read in admixture data
  cl = read.csv("data/admixture/best_runs.csv")
  res = vector("list", nrow(cl))
  for (i in 1:nrow(cl)) {
    run = cl[i, "run"]
    k = cl[i, "k"]
    group = cl[i, "group"]
    a = read.table(paste0("data/admixture/", group, ".", run, ".miss0.7.MAC2.thinned.", k, ".Q"), sep = " ")
    n = read.table(paste0("data/input_files/", group, ".0.miss0.7.MAC2.thinned.fam"))
    a$ind = n$V2
    a[] <- lapply(a, function(a) ifelse(a < 0.001, 0.0001, a))
  
    xx = gather(a, "pop", "prob", -ind) %>% filter(ind %in% inds)
    pops = xx %>% group_by(ind) %>% slice_max(prob) %>%
      dplyr::select(ind, pop) %>% ungroup() %>% 
      left_join(x %>% dplyr::select(ind = sample_fix, sp = cluster))
    # get group to species identify
    popsp = pops %>% group_by(sp, pop) %>% summarize(n = n()) %>% slice_max(n)
    res[[i]] = xx %>% left_join(popsp)
  }
  
  xx = do.call("rbind", res)
  order = t3$tip.label
  sps = unique(xx$sp)
  xx$ind = factor(xx$ind, levels = order)
  
  admixplot = ggplot(data = xx, aes(y = ind)) +
    geom_bar(aes(x = prob, fill = sp), 
             color = "white",
             size = 0.01, 
             stat = "identity", position = "fill",
             show.legend = FALSE) +
    scale_fill_manual(values = cl_cols[sps]) + 
    labs(x = "", y = "") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0), breaks = round(c(0, 1), 0)) +
    theme_minimal() +
    theme(
      axis.text.y = element_blank(), 
      axis.text.x = element_blank(), 
      axis.title.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid = element_blank(), 
      axis.title.x = element_blank(),
      plot.margin = margin(r = 15, l = 10, t = 0, b = 0))
    
  return(admixplot)
}

### Plotting the map ----

plot_map = function(t3) {
  ll = x[x$OTU %in% unlist(otus), ]
  
  # drop weird sample
  ll = ll[ll$sample_fix != "LACM136734", ]
  
  # read in all mapping data
  world = ne_countries(scale = "medium", returnclass = "sf")
  
  # sps = unique(ll$OTU)
  # names2 = unlist(lapply(sps, function(x) {paste0(unique(ll[which(ll$cluster == x), "OTU"]), collapse = ", ")}))
  # ll$cluster2 = names2[match(ll$cluster, sps)]
  ll$OTU2 = factor(ll$OTU, levels = rev(unique(ll[match(t3$tip.label, ll$sample_fix), "OTU"])))

  mapplot = ggplot() +
    geom_sf(data = world, fill = "white", color = "gray70", size = 0.15) +
    geom_point(data = ll, aes(x = longitude, y = latitude, fill = OTU), 
               alpha = 0.7, size = 2, color = "black", stroke = 0.2, shape = 21) +
    facet_wrap(~OTU2, ncol = 3) +
    xlim(-121, -102) + ylim(21, 38) +
    theme_void() +
    guides(fill = "none", size = "none", alpha = "none", shape = "none") +
    scale_fill_manual(values = otu_cols[unique(ll$OTU)]) +  
    theme(panel.background = element_rect(fill = "aliceblue"),
          panel.border = element_rect(fill = NA, color = "gray30", size = 0.5),
          strip.text = element_text(size = 9, face = "italic", 
                                    margin = margin(t = 0, r = 0, b = 2, l = 0)), 
          strip.background = element_blank(),
          plot.margin = margin(r = 2, l = 2, t = 0, b = 0))
  return(mapplot)
}

otus1 = unlist(otus)
# drop duplicated inds
x = x[x$duplicate == FALSE, ]
inds = x[x$OTU %in% otus1, "sample_fix"]

trees = plot_tree()
admixplot = plot_admixture(trees[[1]])
mapplot = plot_map(trees[[1]])
 
allplt = cowplot::plot_grid(trees[[2]], admixplot, mapplot, ncol = 3,
                            rel_widths =  c(2, 1, 2))
cowplot::save_plot("figures/all_ind_nuc_tree.png", allplt, 
                   base_height = 8, base_width = 10)

