rm(list = ls())

library(ape)
library(tidyverse)
library(ggplot2)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Xantusia/")
source("scripts/gps_colors.R")

# plot it!!!
t = treeio::read.beast("data/phylogeny/mcc.tre")
nodelabs = t@data$posterior
names(nodelabs) = t@data$node
t = t@phylo
nodelabs = nodelabs[ as.character(seq(Ntip(t) + 1, Ntip(t) + Nnode(t))) ]
t$node.label = nodelabs
t = read.tree(text = write.tree(ladderize(t)))
t1 = t
t1$tip.label = paste("X.", t1$tip.label)

pdf("figures/species_treephylogeny.pdf", width = 3.5, height = 4)
par(mar = c(0, 0, 0, 6), xpd = T)
plot.phylo(t1, edge.width = 2, show.tip.label = F)
pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)

for (i in 1:length(nodelabs)) {
  if (nodelabs[i] >= 0.95) {
    nodelabels("", i + Ntip(t), pch = 16, col = "black", size = 1,
               frame = "none")
  }
}

boxlabel<-function(x,y,text,cex=1,bg="transparent",offset = 0){
  w<-strwidth(text)*cex*1.1
  h<-strheight(text)*cex*1.4
  os<-offset*strwidth("W")*cex
  rect(x+os,y-0.5*h,x+w+os,y+0.5*h,col=bg,border=0)
  text(x,y,text,pos=4,offset=offset,font=3)
}

c = combo %>% dplyr::select(species1, habitat1) %>% distinct()
tipbg = rep(NA, Ntip(t1))
tipbg[which(t$tip.label %in% c[c$habitat1 == "rocks", "species1"])] = habcols["rock-rock"]
tipbg[which(t$tip.label %in% c[c$habitat1 == "plants", "species1"])] = habcols["plant-plant"]
# tiplabels(t1$tip.label, seq(1, Ntip(t1)), , font = 3,
#          bg = scales::alpha(tipbg, 0.5), 
#          adj = c(0, 0.5), col = NA)

for(i in 1:Ntip(t1)) boxlabel(pp$xx[i],pp$yy[i],t1$tip.label[i], 
                              bg= alpha(tipbg[i], 0.5), offset = 0.5)
dev.off()
