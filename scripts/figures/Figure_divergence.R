rm(list = ls())

library(ape)
library(tidyverse)
library(ggplot2)
library(cowplot)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Xantusia/")
source("scripts/gps_colors.R")

dd = list.files("data/divergence/", full.names = T)
dd1 = lapply(dd, read.csv)
dd2 = do.call("rbind", dd1)

fst = dd2 %>% filter(metric == "fst", submetric == "fst") %>%
  spread(submetric, value) %>% dplyr::select(species1, species2, fst)
# check d_a
dxy = dd2 %>% filter(metric == "dxy") %>%
  spread(submetric, value) %>% 
  mutate(d_xy = diff / denom, d_a = (diff - (pi1 + pi2) / 2) / denom) %>%
  dplyr::select(species1, species2, d_xy, d_a)
mtdxy = dd2 %>% filter(metric == "mtdxy", submetric == "mt_btn") %>%
  spread(submetric, value) %>% dplyr::select(species1, species2, mt_dxy = mt_btn)

div = left_join(dxy, fst) %>% left_join(mtdxy) %>%
  left_join(combo)

# add group identity
get_group <- function(y) {
  return(names(which(unlist(lapply(otus, function(x) y %in% x)) == TRUE)))
}
div$group1 = sapply(div$species1, get_group)
div$group2 = sapply(div$species2, get_group)

div2 = div[div$group1 == div$group2, ]


t = ape::read.nexus("data/phylogeny/mcc.tre")
td = cophenetic.phylo(t) %>% as.data.frame() 
td$species1 = rownames(td)
td = td %>% gather(species2, dist, -species1) %>% mutate(dist = dist / 2)

div3 = div2 %>% left_join(td)


a = ggplot(div3, aes(fst, d_xy)) + 
  geom_point(aes(fill = type), shape = 21, size = 2.5) + 
  ylab(expression(d[xy])) + xlab(expression(F[ST])) +
  theme_classic() + 
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = habcols)
b = ggplot(div3, aes(fst, dist)) + 
  geom_point(aes(fill = type), shape = 21, size = 2.5) + 
  xlab(expression(F[ST]))  + ylab("relative div. time") +
  theme_classic() + 
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = habcols)
c = ggplot(div3, aes(fst, mt_dxy)) + 
  geom_point(aes(fill = type), shape = 21, size = 2.5) + 
  xlab(expression(F[ST]))  + ylab(expression("mtDNA" ~ d[xy])) +
  theme_classic() + 
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = habcols)

prow <- plot_grid(
  a + theme(legend.position="none"),
  b + theme(legend.position="none"),
  c + theme(legend.position="none"),
  align = 'vh',
  labels = c("A", "B", "C"),
  hjust = -1,
  nrow = 1
)

legend <- get_legend(
  # create some space to the left of the legend
  a + theme(legend.box.margin = margin(0, 0, 0, -10))
)

abc = plot_grid(prow, legend, rel_widths = c(3, .5))
save_plot("figures/divergence.png", abc, base_height = 2, base_width = 8)




d = read.csv("data/demographic_results.16Jan24.csv")
div3 = div2 %>% left_join(d)

div3 = div3 %>% arrange(fst) %>%
  mutate(id = 1:nrow(div3))

ggplot(div3, aes(id, fst)) + 
  geom_point(aes(fill = type), shape = 21, size = 2.5) +
  theme_classic() + 
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = habcols)

box = data.frame(x = c(0.005, 0.02, 0.02, 0.005),
                 y = c(0, 0, 1, 1))

migr = ggplot(div3, aes(d_a, p_migration_no_sfs)) +
  geom_polygon(data = box, aes(x, y), alpha = 0.5, fill = "gray") +
  geom_point(aes(fill = type), shape = 21, size = 2.5) +
  scale_x_log10() +
  xlab(expression(d[A])) + ylab("p(migration)") +
  scale_fill_manual(values = habcols) + 
  theme_classic() + 
  theme(legend.title = element_blank(),
        legend.position = c(.85, .85)) +
  ylim(0, 1)
save_plot("figures/migration.png", migr, base_height = 3, base_width = 5)



