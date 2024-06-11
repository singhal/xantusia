library(ape)

source("~/Dropbox (Personal)/research_projects-ACTIVE/Xantusia/scripts/gps_colors.R")

t = read.tree("~/Dropbox (Personal)/research_projects-ACTIVE/Xantusia/data/phylogeny/xantusia.low_inds_removed.treefile")
d = read.csv("~/Dropbox (Personal)/research_projects-ACTIVE/Xantusia/metadata/xantusia_samples_v7.csv")

outs = d[d$genus == "Lepidophyma", "sample_fix"]
outs = outs[outs %in% t$tip.label]
t1 = root(t, outs, resolve.root = TRUE)
t2 = drop.tip(t1, outs)

d1 = d[complete.cases(d$OTU), ]
sps = unique(d1$OTU)
keep = rep(NA, length(sps))
for (i in 1:length(sps)) {
  sp = sps[i]
  inds = d[which(d$OTU == sp), "sample_fix"]
  inds = inds[inds %in% t$tip.label]
  keep[i] = sample(inds, 1)
}

t3 = keep.tip(t2, keep)
t3$tip.label = sps[match(t3$tip.label, keep)]
t3$node.label = NULL

d2 = d[d$genus == "Xantusia", ]
d2[!complete.cases(d2$OTU), "OTU"] = "xxx"
d2[!(d2$sample_fix %in% t$tip.label), "OTU"] = "xxx"

outs = c("vigilis_KCAV", "vigilis_KCAV", "henshawi_2")
for (i in 1:length(groups)) {
  # identify keep
  keep = c(outs[i], groups[[i]])
  
  # make species tree
  t4 = keep.tip(t3, keep)
  
  # modify sets
  d3 = d2
  d3[!d3$OTU %in% keep, "OTU"] = "xxx"
  d3[d3$OTU == outs[i], "OTU"] = "Outgroup"
  
  write.tree(t4, paste0("~/Desktop/species_tree.", names(groups)[i], ".tre"))
  write.table(d3[,c("sample_fix", "OTU")], paste0("~/Desktop/SETS.", names(groups)[i], ".txt"),
              quote = F, row.names = F, sep = "\t", col.names = F)
}
