rm(list = ls())

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Xantusia/")
sampfile = "metadata/xantusia_samples_v8.csv"
indfile =  "data/archived/xantusia.1.miss0.7.MAC2.thinned.fam"
snpfile =  "data/archived/xantusia.1.miss0.7.MAC2.thinned.snp"

t = read.snp(snpfile)

d = read.csv(sampfile, stringsAsFactors = F, na.string = c(""))
inds = read.table(indfile, stringsAsFactors = F, header = F, sep = " ")
inds = data.frame(inds = inds$V2, 
                  sp = d[match(inds$V2, d$sample_fix), "OTU"], 
                  miss = unlist(lapply(t@gen, function(x) length(x@NA.posi))) / t@n.loc)

dups = d[d$duplicate == TRUE, "sample_fix"]
inds = inds[inds$sp != "NA", ]
inds = inds[!inds$inds %in% dups, ]

otus = unique(inds$sp)
sinds = vector("list", length(otus))
for (i in 1:length(otus)) {
  dd = inds[which(inds$sp == otus[i]), ]
  dd1 = dd[dd$miss < 0.2, ]
  if (nrow(dd1) > 4) {
    sinds[[i]] = sample(dd1$inds, 4)
  } else if (nrow(dd) > 4) {
    sinds[[i]] = sample(dd$inds, 4)
  } else {
    sinds[[i]] = dd$inds
  }
}

xx = data.frame(ind = unlist(sinds))
write.csv(xx, "~/Desktop/subsampled.csv", row.names = F, quote = F)

View(cbind(xx, inds[ match(xx$ind, inds$inds), "miss"]))

library(ape)
a = read.dna("~/Desktop/subsampled.usnps", format = "sequential")
dimnames(a)[[1]] = paste0(dimnames(a)[[1]], "_" , d[match(dimnames(a)[[1]], d$sample_fix), "OTU"])
write.dna(a, "~/Desktop/subsampled.renamed.usnps", format = "sequential")
