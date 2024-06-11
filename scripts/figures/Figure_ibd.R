rm(list = ls())

library(sf)
library(tidyverse)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Xantusia/")
source("scripts/gps_colors.R")

x = read.csv("metadata/xantusia_samples_v9.csv")

pts = x %>% filter(complete.cases(longitude)) %>% 
  st_as_sf(coords = c("longitude", "latitude")) %>% 
  st_set_crs(4979)
d = as.matrix(pts %>% st_distance())
dimnames(d)[[1]] = pts$sample_fix
dimnames(d)[[2]] = pts$sample_fix

d1 = data.frame(ind1 =colnames(d)[col(d)], ind2=rownames(d)[row(d)], dist=c(d))
d1$sp1 = x[match(d1$ind1, x$sample_fix), "OTU"]
d1$sp2 = x[match(d1$ind2, x$sample_fix), "OTU"]

d2 = read.csv("data/divergence_between_inds.16Jan23.csv") 
d3 = d2 %>% left_join(d1)

# anyone missing?
d3 %>% filter(!complete.cases(dist)) %>% nrow()

# filter out low values
d4 = d3 %>% filter(fst_denom >= 100) %>% 
  mutate(dist = as.numeric(dist), log_dist = log(dist),
         inv_fst = fst / (1 - fst))

# who am i losing?
# this is pretty much just dropping XAHE14
d3 %>% filter(fst_denom < 100) %>% View()

keepsp = d4 %>% filter(sp1 == sp2) %>% 
  group_by(sp1) %>% 
  summarize(n = length(unique(c(ind1, ind2)))) %>%
  filter(n >= 4) %>% pull(sp1)
# drop bolsonae b/c sampling too linear
keepsp = keepsp[keepsp != "bolsonae"]

lms = vector(mode = "list", length = length(keepsp))
for (i in 1:length(keepsp)) {
  dx = d4 %>% filter(sp1 == sp2, sp1 %in% keepsp[i]) %>%
    filter(is.finite(log_dist))
  lms[[i]] = lm(dx$inv_fst ~ dx$log_dist)
}

# pvals 
pvals = unlist(lapply(lms, function(x) summary(x)$coefficients[2, 4]))
names(pvals) = keepsp


# rsq
rsq = unlist(lapply(lms, function(x) summary(x)$adj.r.squared))
names(rsq) = keepsp

# change title, change color
d5 = d4 %>% filter(sp1 == sp2, sp1 %in% keepsp) %>%
  filter(is.finite(log_dist)) %>%
  mutate(species = gsub("_", " ", sp1))
pltibd = ggplot(d5,
       aes(log_dist, inv_fst)) +
  geom_point(alpha = 0.5, shape = 16, color = habcols[1]) + 
  geom_smooth(data = subset(d5, sp1 %in% names(pvals[pvals < 0.05])), 
              method = "lm", color = "black", size = 0.5) +
  facet_wrap(~species, ncol = 4, scales = "free") +
  xlab("log distance (km)") +
  ylab(expression("inverse" ~ F[ST])) +
  theme_classic() +
  theme(
    strip.background = element_rect(
      color="white", fill="white", size=0, linetype="solid"
    ),
    strip.text.x = element_text(
      size = 10, face = "italic", hjust = 0
    )
  )
cowplot::save_plot("figures/within_IBD.png", pltibd, base_height = 5, 
                   base_width = 8)




cc = list(
  c("henshawi_1", "henshawi_2"),
  c("henshawi_2", "henshawi_3"),
  c("bezyi", "wigginsi"),
  c("vigilis_SJ", "vigilis_YV"),
  c("sierra", "vigilis_OV"),
  c("vigilis_EM", "vigilis_KCAV"))

res = vector("list", length(cc))
  
for (j in 1:length(cc)) {
    a = cc[[j]][1]
    b = cc[[j]][2]
    dd = d4 %>% filter((sp1 == a & sp2 == a) | 
                         (sp1 == b & sp2 == b) |
                         (sp1 == a & sp2 == b) |
                         (sp1 == b & sp2 == a))
    dd$type = NA
    dd[dd$sp1 == a & dd$sp2 == a, "type"] = paste0("within")
    dd[dd$sp1 == b & dd$sp2 == b, "type"] = paste0("within")
    dd[dd$sp1 == a & dd$sp2 == b, "type"] = paste0("between")
    dd[dd$sp1 == b & dd$sp2 == a, "type"] = paste0("between")
    res[[j]] = 
      ggplot(dd, aes(log_dist, inv_fst)) + geom_point(aes(fill = type), pch = 21, stroke = 0.3) +
      xlab("log distance (km)") +
      ylab(expression("inverse" ~ F[ST])) +
      ggtitle(paste0(a, " vs. ", b)) +
      scale_fill_manual(values = c("gray20", "white")) +
      theme_classic() +
      theme(legend.position = "none")
}
  

xx = cowplot::plot_grid(plotlist = res, ncol = 3)
cowplot::save_plot("figures/IBD.png", xx, base_height = 5, 
                   base_width = 10)
