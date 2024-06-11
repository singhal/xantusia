rm(list = ls())
library(adegenet)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

source("scripts/gps_colors.R")

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Xantusia/")
sampfile = "metadata/xantusia_samples_v10.csv"
indfile =  "data/archived/xantusia.1.miss0.7.MAC2.thinned.fam"
snpfile =  "data/archived/xantusia.1.miss0.7.MAC2.thinned.snp"

d = read.csv(sampfile, stringsAsFactors = F, na.string = c(""))
inds = read.table(indfile, stringsAsFactors = F, header = F, sep = " ")
inds = data.frame(inds = inds$V2, 
                  sp = d[match(inds$V2, d$sample_fix), "OTU"], 
                  stringsAsFactors = F)
rownames(inds) = inds$inds

##############
# PCA results
##############

t = read.snp(snpfile)
pcasnp = glPca(t, nf=4)

eig.perc <- pcasnp$eig/sum(pcasnp$eig)
(eig.perc[1:10])
pcal = data.frame(pcasnp$scores)
pcal$sp = inds$sp

pcal1 = pcal[pcal$sp != "NA", ]
pcal1$group = NA
pcal1$group = ifelse(pcal1$sp %in% otus[["henshawi"]], "group 1", pcal1$group)
pcal1$group = ifelse(pcal1$sp %in% otus[["extorris"]], "group 2", pcal1$group)
pcal1$group = ifelse(pcal1$sp %in% otus[["vigilis"]], "group 3", pcal1$group)

order = pcal1 %>% arrange(group) %>% dplyr::select(sp) %>% unique() %>% pull(sp)
pcal1$sp2 = factor(pcal1$sp, levels = order)

xx = ggplot(pcal1, aes(PC1, PC2, fill = sp2)) + 
  ggforce::geom_mark_ellipse(aes(x = PC1, y = PC2, fill = group), 
                             expand = unit(2, "mm"), size = 0.1, alpha = 0.3) +
  geom_jitter(shape = 21, size = 2.5) +
  theme_cowplot() +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = otu_cols) +
  xlab(paste0("PC1 (", round(eig.perc[1] * 100, 1), "%)")) +
  ylab(paste0("PC1 (", round(eig.perc[2] * 100, 1), "%)"))
cowplot::save_plot("figures/PCA.pdf", xx, base_height = 5, base_width = 7)

