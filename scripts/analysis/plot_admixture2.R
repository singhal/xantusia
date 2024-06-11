rm(list = ls())
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(zoo)

cols = RColorBrewer::brewer.pal(12, "Paired")
names(cols) = paste0("V", seq(1:12))
group = "vigilis"

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Xantusia/")
d = read.csv(paste0("data/admixture/", group, ".summary.csv"))

# get last k with decreasing CV value
best = d %>% arrange(run, k) %>% group_by(run) %>% 
  mutate(diff = lag(cv) - cv) %>% 
  dplyr::filter(diff < 0) %>%
  dplyr::filter(row_number() == 1) %>% 
  ungroup() %>% dplyr::select(run, k) %>%
  mutate(k = k - 1)
best$best = "yes"
d = left_join(d, best)
d[!complete.cases(d$best), "best"] = "no"

table(best$k)

ggplot(d, aes(k, cv)) + geom_line() + facet_wrap(~run, ncol = 3) +
  theme_classic() + geom_point(aes(fill = best), shape = 21) +
  scale_fill_manual(values = c("black", "red"))

d = read.csv("metadata/xantusia_samples_v8.csv")

kk = best %>% group_by(k) %>% slice_max(run)
pdf(paste0("figures/", group, ".admixture.pdf"), height = 4, width = 10)
for (i in 1:nrow(kk)) {
  file = paste0("data/admixture/", group, ".", kk[i, "run"], ".miss0.7.MAC2.thinned.",  kk[i, "k"], ".Q")
  x = read.table(file)
  n = read.table(paste0("data/input_files/", group, ".", kk[i, "run"], ".miss0.7.MAC2.thinned.fam"))
  x$ind = n$V2
  x$sp = d[match(x$ind, d$sample_fix), "revised_species"]
  x$sp = substr(x$sp, 1, 3)
  
  x[] <- lapply(x, function(x) ifelse(x < 0.001, 0.0001, x))
  
  xx = gather(x, "pop", "prob", -sp, -ind)
  pops = xx %>% group_by(ind) %>% slice_max(prob) %>%
    dplyr::select(ind, pop) %>% ungroup()
  x2 = x %>% left_join(pops)
  
  order = x2 %>% arrange(pop, across(names(x)[1:(ncol(x) - 2)])) %>% 
    ungroup() %>% dplyr::select(ind) %>% pull()
  
  plt = ggplot(xx, aes(factor(ind, levels = order), prob, fill = pop)) +
    geom_col(color = "gray", size = 0.1) +
    facet_grid(~sp, switch = "x", scales = "free", space = "free") +
    theme_minimal() + labs(x = "",  y = "") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = expand_scale(add = 1)) +
    theme(
      panel.spacing.x = unit(0.1, "lines"),
      axis.text.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none"
    ) +
    scale_fill_manual(values = cols)
  print(plt)
}
dev.off()

bestk = names(sort(table(best$k),decreasing=TRUE)[1])
kk2 = best %>% filter(k == as.numeric(bestk)) 
pdf(paste0("preliminary_figures/", group, ".best", bestk, ".admixture.pdf"), height = 4, width = 10)
for (i in 1:nrow(kk2)) {
  file = paste0("data/admixture/", group, ".", kk2[i, "run"], ".miss0.7.MAC2.thinned.",  kk2[i, "k"], ".Q")
  x = read.table(file)
  n = read.table(paste0("data/input_files/", group, ".", kk2[i, "run"], ".miss0.7.MAC2.thinned.fam"))
  x$ind = n$V2
  x$sp = d[match(x$ind, d$sample_fix), "revised_species"]
  x$sp = substr(x$sp, 1, 3)
  
  x[] <- lapply(x, function(x) ifelse(x < 0.001, 0.0001, x))
  
  xx = gather(x, "pop", "prob", -sp, -ind)
  pops = xx %>% group_by(ind) %>% slice_max(prob) %>%
    dplyr::select(ind, pop) %>% ungroup()
  x2 = x %>% left_join(pops)
  
  order = x2 %>% arrange(pop, across(names(x)[1:(ncol(x) - 2)])) %>% 
    ungroup() %>% dplyr::select(ind) %>% pull()
  
  plt = ggplot(xx, aes(factor(ind, levels = order), prob, fill = pop)) +
    geom_col(color = "gray", size = 0.1) +
    facet_grid(~sp, switch = "x", scales = "free", space = "free") +
    theme_minimal() + labs(x = "",  y = "") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = expand_scale(add = 1)) +
    theme(
      panel.spacing.x = unit(0.1, "lines"),
      axis.text.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none"
    ) +
    scale_fill_manual(values = cols)
  print(plt)
}
dev.off()


ele = terra::rast("data/rasters/HYP_HR_SR_W/HYP_HR_SR_W.tif")

pdf(paste0("figures/", group, ".map.pdf"), width = 4, height = 5)
for (i in 1:nrow(kk)) {
  file = paste0("data/admixture/", group, ".", kk[i, "run"], ".miss0.7.MAC2.thinned.",  kk[i, "k"], ".Q")
  x = read.table(file)
  n = read.table(paste0("data/input_files/", group, ".", kk[i, "run"], ".miss0.7.MAC2.thinned.fam"))
  x$ind = n$V2
  x$sp = d[match(x$ind, d$sample_fix), "revised_species"]
  x$sp = substr(x$sp, 1, 3)
  
  x[] <- lapply(x, function(x) ifelse(x < 0.001, 0.0001, x))
  
  xx = gather(x, "pop", "prob", -sp, -ind)
  xx1 = xx %>% group_by(ind) %>% slice_max(prob) %>% ungroup()
  xx2 = left_join(xx1, d, by = c("ind" = "sample_fix"))
  
  ymin = min(xx2$latitude, na.rm = T) - 1
  ymax = max(xx2$latitude, na.rm = T) + 1
  xmin = min(xx2$longitude, na.rm = T) - 1
  xmax = max(xx2$longitude, na.rm = T) + 1
  
  ele2 = raster::crop(ele, raster::extent(xmin, xmax, ymin, ymax))

  par(mar = c(0, 0, 0, 0))
  plot(NA, xlim = c(xmin, xmax), ylim = c(ymin, ymax), 
       axes = F, xlab = "", ylab = "")
  terra::plotRGB(ele2, add = T,maxpixels = raster::ncell(ele2))
  points(xx2$longitude, xx2$latitude, bg = cols[ xx2$pop ],
         pch = 21, cex = 1)
}
dev.off()

# world = rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
# ggplot() +
#     geom_sf(data = world) +
#     geom_point(data = xx2, 
#                 aes(x = longitude, y = latitude, fill = pop), 
#                 shape = 21, alpha = 0.99,
#                 size = 2, color = "black", stroke = 0.2) +
#     coord_sf(xlim = c(-118, -115), ylim = c(30, 35)) +
#     theme_void() +
#     scale_fill_manual(values= cols) +
#     guides(fill = "none", size = "none", alpha = "none", shape = "none")

