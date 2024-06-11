rm(list = ls())

library(ape)
library(phytools)
library(phangorn)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Xantusia/")
t1 = read.nexus("data/phylogeny/mcc.tre")
t = read.nexus("~/Dropbox (Personal)/research_projects-ACTIVE/NSF_NALDP/Species_Specific/Xantusia/for_Adam_Hayden/species_tree_snapp/combined_runs/combined.trees")
h = read.csv("data/BayesTraits/Xantusia.habitat.txt", header = F, sep = "\t")
h1 = h$V2
names(h1) = h$V1
h1[which(h1 == 0)] = "rocks"
h1[which(h1 == 0)] = "plants"

er = fitMk(t1, h1, model = "ER")
ard = fitMk(t1,  h1, model = "ARD")
fit01 = fitMk(t1,  h1, model = matrix(c(0, 1, 0, 0), 2, 2, byrow = T))
fit10 = fitMk(t1,  h1, model = matrix(c(0, 0, 1, 0), 2, 2, byrow = T))

# no reason to not use ER, which has greatest log likelihood
unlist(lapply(list(er, ard, fit01, fit10), function(x) logLik(x)))

mtrees = make.simmap(t[1:100], h1, model = "ER", nsim = 10,
            Q = "mcmc", vQ = 0.01,
            prior = list(use.empirical = TRUE), samplefreq = 10)
summary(mtrees)

mtrees1 = make.simmap(t1, h1, model = "ER", nsim = 1000,
                     Q = "mcmc", vQ = 0.01,
                     prior = list(use.empirical = TRUE), samplefreq = 10)
source("./scripts/gps_colors.R")
cols = habcols[c("rock-rock", "plant-plant")]
mtree1 = densityMap(mtrees1)

pd = summary(mtrees1)
mtree1 = setMap(mtree1, rev(cols))
png("figures/habitat_map.png", res = 300, height = 3, width = 6, units = "in")
mtree1[[1]]$tip.label[which(mtree1[[1]]$tip.label == "magdalena")] = "sherbrookei"
plot(mtree1,cols,fsize=0.8,ftype="i")
for (i in 1:Nnode(mtree1[[1]])) {
  nnode = i + Ntip(mtree1[[1]])
  nodelabels("", nnode, frame = "none", cex = 0.7, pie = pd$ace[as.character(nnode), ], piecol = rev(cols))
}
dev.off()
