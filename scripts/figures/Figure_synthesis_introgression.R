rm(list = ls())

library(ape)
library(tidyverse)
library(ggplot2)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Xantusia/")
source("scripts/gps_colors.R")

###############
# DILS results
###############

# add group identity
get_group <- function(y) {
  return(names(which(unlist(lapply(otus, function(x) y %in% x)) == TRUE)))
}

d = read.csv("data/demographic_results.16Jan24.csv")
d[d$species1 == "magdalena", "species1"] = "sherbrookei"
d[d$species2 == "magdalena", "species2"] = "sherbrookei"
d$group1 = sapply(d$species1, get_group)
d$group2 = sapply(d$species2, get_group)
d = d %>% filter(group1 == group2) %>%
  mutate(analysis = "DILS") %>%
  filter(p_migration_no_sfs > 0.5) %>%
  dplyr::select(analysis, species1, species2)

###############
# mtDNA & admixture
###############

x = read.csv("data/introgression.csv") %>% rbind(d)

###############
# dstat
###############

res = vector("list", length(clusters))
for (i in 1:length(clusters)) {
  s = paste0("data/dsuite/", names(clusters)[i], "_tree.txt")
  res[[i]] = read.csv(s, sep = "\t")
}

# In all three types of output, we order P1 and P2 so that nABBA>= nBABA.
# As a result, the D statistic is always positive and all the results, 
# including the f4-ratio and other statistics reflect
# evidence of excess allele sharing between P3 and P2 for each trio.
res = do.call("rbind", res)
res1 = res %>% filter(p.value < 0.05 / nrow(res)) %>% 
  dplyr::select(species1 = P2, species2 = P3) %>% distinct()

res = vector("list", length(clusters))
for (i in 1:length(clusters)) {
  s = paste0("data/dsuite/", names(clusters)[i], ".fbranch")
  s = read.csv(s, sep = "\t")
  
  # don't map internal branches
  s1 = s[grep(",", s$branch_descendants, invert = T), ]
  # keeps outgroup, but this is okay b/c all NAs
  res[[i]] = s1 %>% dplyr::select(-branch) %>% 
    gather(key, value, -branch_descendants) 
}
res1 = do.call("rbind", res) %>% 
  filter(complete.cases(value))
hist(res1$value)
res2 = res1 %>% filter(value > 0.05) %>%
  mutate(analysis = "dstat") %>%
  dplyr::select(analysis, species1 = branch_descendants, species2 = key)

x = x %>% rbind(res2)
x$sp1 = ifelse(x$species1 > x$species2, x$species1, x$species2)
x$sp2 = ifelse(x$species1 > x$species2, x$species2, x$species1)
x2 = x %>% dplyr::select(analysis, sp1, sp2) %>% group_by(analysis) %>%
  distinct() %>% ungroup()

# plot it!!!
t = treeio::read.beast("data/phylogeny/mcc.tre")
nodelabs = t@data$posterior
names(nodelabs) = t@data$node
t = t@phylo
nodelabs = nodelabs[ as.character(seq(Ntip(t) + 1, Ntip(t) + Nnode(t))) ]
t$node.label = nodelabs
t = read.tree(text = write.tree(ladderize(t)))

cols = cl_cols[1:4]
names(cols) = c("mtDNA", "admixture", "dstat", "DILS")
lbltxt = c("mtDNA introgression", "admixture", "d-statistic", "demographic modeling")

anc = phangorn::Ancestors(t, 1:Ntip(t), type = "parent")
hts = unlist(lapply(anc, function(x) phytools::nodeheight(t, x)))
names(hts) = t$tip.label

x2$ht1 = hts[ x2$sp1 ]
x2$ht2 = hts[ x2$sp2 ]
x2$ht = apply(x2, 1, function(x) { max(x[4:5])})

x2 = x2 %>% group_by(ht) %>% arrange(sp1, sp2) %>% 
  mutate(id = row_number()) %>% ungroup()

t1 = t
t1$tip.label = paste("X.", t1$tip.label)
t1$tip.label[which(t1$tip.label == "X. magdalena")] = "X. sherbrookei"

pdf("figures/introgression_phylogeny.pdf", width = 8, height = 4)
par(mar = c(0, 0, 0, 6), xpd = T)
plot.phylo(t1, edge.width = 2, show.tip.label = F)
pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)

boxlabel<-function(x,y,text,cex=1,bg="transparent",offset = 0){
  w<-strwidth(text)*cex*1.1
  h<-strheight(text)*cex*1.4
  os<-offset*strwidth("W")*cex
  rect(x+os,y-0.5*h,x+w+os,y+0.5*h,col=bg,border=0)
  text(x,y,text,pos=4,offset=offset,font=3)
}

c = combo %>% dplyr::select(species1, habitat1) %>% distinct()
tipbg = rep(NA, Ntip(t1))
t$tip.label[which(t$tip.label == "magdalena")] = "sherbrookei"
tipbg[which(t$tip.label %in% c[c$habitat1 == "rocks", "species1"])] = habcols["rock-rock"]
tipbg[which(t$tip.label %in% c[c$habitat1 == "plants", "species1"])] = habcols["plant-plant"]
# tiplabels(t1$tip.label, seq(1, Ntip(t1)), , font = 3,
#          bg = scales::alpha(tipbg, 0.5), 
#          adj = c(0, 0.5), col = NA)

for(i in 1:Ntip(t1)) boxlabel(pp$xx[i],pp$yy[i],t1$tip.label[i], 
                              bg= alpha(tipbg[i], 0.5), offset = 0.5)


xlim = max(pp$xx)
for (i in 1:length(t$node.label)) {
  if (as.numeric(t$node.label[i]) > 0.95) {
    nodenum = i + Ntip(t)
    nodelabels("", nodenum, pch = 16, col = "black", frame = "none")
  }
}
x2[x2$sp1 == "magdalena", "sp1"] = "sherbrookei"
x2[x2$sp2 == "magdalena", "sp2"] = "sherbrookei"
for (i in 1:nrow(x2)) {
  sp1 = pull(x2[i, "sp1"])
  sp2 = pull(x2[i, "sp2"])
  ystart = match(sp1, t$tip.label)
  yend = match(sp2, t$tip.label)
  
  ys = min(c(ystart, yend))
  ye = max(c(ystart, yend))
  
  ht = as.numeric(pull(x2[i, "ht"]))
  
  n = pull(x2[i, "id"])
  diff = xlim - ht
  xloc = ht + n * 0.15 * diff
  
  lines(x = c(xloc, xloc),
        y = c(ys + 0.05, ye - 0.05), col = cols[pull(x2[i, "analysis"])])
}
legend(0, 5, lbltxt, text.col = cols, bty = "n", xjust = 0)
dev.off()


###########################################
# compare patterns of introgression to null
###########################################

# combine x2 with phy dist & group data & hab type
x3 = x2 %>% select(species1 = sp1, species2 = sp2) %>% 
  unique() %>% left_join(combo)
x3$group = names(unlist(lapply(x3$species1, function(x) {which(unlist(lapply(otus, function(y) {x %in% y})) == TRUE)})))

t = ape::read.nexus("data/phylogeny/mcc.tre")
dd = as.data.frame(ape::cophenetic.phylo(t)) 
dd$species1 = rownames(dd)
dd1 = dd %>% tidyr::gather(species2, phy_dist, -species1)

x4 = x3 %>% left_join(dd1)

### do 1000 simulations
nsim = 1000
res = vector("list", nsim)
g = table(x4$group)
for (i in 1:nsim) {
  gx = vector("list", length(g))
  for (j in 1:length(g)) {
    gxx = lapply(1:g[j], function(x) sample(otus[[names(g)[j]]], 2))
    gx[[j]] = as.data.frame(do.call(rbind, gxx))  
  }
  gxx = do.call(rbind, gx)
  names(gxx) = c("species1", "species2")
  res[[i]] = gxx %>% left_join(dd1)  %>% left_join(combo)
}

# get average phy dist relative to ours
mean_phy_dist = unlist(lapply(res, function(x) {mean(x$phy_dist)}))
a = ggplot(data.frame(mean_phy_dist), aes(mean_phy_dist)) +
  geom_histogram() + theme_classic() +
  geom_vline(xintercept = mean(x4$phy_dist), col = "red") +
  xlab("mean phylogenetic distance")

# get average % btn eco relative to ours
mean_btn = unlist(lapply(res, function(x) {sum(x$type == "rock-plant") / nrow(x)}))
our_btn = sum(x4$type == "rock-plant") / nrow(x4)

b = ggplot(data.frame(mean_btn), aes(mean_btn)) +
  geom_histogram(bins = 10) + theme_classic() +
  geom_vline(xintercept = our_btn, col = "red") +
  xlab("% of gene flow events between morphs")

ab = plot_grid(a,b, labels = c("A", "B"))
save_plot("figures/introgression_simulations.png", ab,
          base_height = 3, base_width = 8)
