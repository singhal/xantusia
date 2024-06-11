library(RColorBrewer)

setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Xantusia/")

otus = list(c("henshawi_1", "henshawi_2", "henshawi_3", "gracilis"),
              c("gilberti", "sherbrookei", "bolsonae", "sanchezi", "extorris"),
              c("bezyi", "wigginsi", "vigilis_SJ", "vigilis_YV",
                "arizonae", "sierra", "vigilis_OV", "vigilis_KCAV",
                "vigilis_EM"))
clusters = list(c("henshawi_1", "henshawi_2", "henshawi_3"),
                c("gilberti", "sherbrookei", "bolsonae", "sanchezi", "extorris"),
                c("bez_wig", "vigilis_1", "arizonae", "vigilis_4", "vigilis_3",
                  "vigilis_2"))
names(otus) = c("henshawi", "extorris", "vigilis")
names(clusters) = c("henshawi", "extorris", "vigilis")

cl_cols = c("#ce437c","#9ab443", "#6092d1",
           "#bc50bb","#57b752", "#7065cc", "#c5733c","#d1453a",
           "#49803d","#bc79b9", "#53b798", "#c76b71","#d5a339", "#92833d")
names(cl_cols) = unlist(clusters)

otu_cols = c("#ce437c","#9ab443", "#6092d1", "#336bb1",
           "#bc50bb","#57b752", "#7065cc", "#c5733c","#d1453a",
           "#49803d","#326826", "#bc79b9", "#982b97", "#53b798", "#ad1d65", "#c76b71","#d5a339", "#92833d")
names(otu_cols) = unlist(otus)

x = read.csv("metadata/xantusia_samples_v10.csv")
x = x[complete.cases(x$OTU), ]
sps = unique(x$OTU)

combo = combn(sps, 2)
combo = data.frame(species1 = c(combo[1, ], combo[2, ]),
                   species2 = c(combo[2, ], combo[1, ]))
combo$habitat1 = x[match(combo$species1, x$OTU), "habitat"]
combo$habitat2 = x[match(combo$species2, x$OTU), "habitat"]
combo$type = NA
combo[(combo$habitat1 == "rocks" & combo$habitat2 == "rocks"), "type"] = "rock-rock"
combo[(combo$habitat1 == "plants" & combo$habitat2 == "plants"), "type"] = "plant-plant"
combo[(combo$habitat1 == "plants" & combo$habitat2 == "rocks"), "type"] = "rock-plant"
combo[(combo$habitat1 == "rocks" & combo$habitat2 == "plants"), "type"] = "rock-plant"


habcols = c("#4daf4a", "#a65628", "#ff7f00")
names(habcols) = c("plant-plant", "rock-rock", "rock-plant")
