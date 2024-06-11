rm(list = ls())

library(tidyverse)

# get directory
setwd("~/Dropbox (Personal)/research_projects-ACTIVE/Xantusia/")

d = read.csv("data/moments/all_moments.10January24.csv")
snp = read.csv("data/moments/sfs_snp_counts.10Jan24.csv")
d = left_join(d, snp)
d = d %>% mutate(AIC = as.numeric(d$AIC)) %>% 
  filter(complete.cases(AIC))

# attempt to get at AIC weights
# this is okay b/c i removed linked sites
wt = d %>% group_by(species1, species2, Model) %>% 
  slice_min(AIC, with_ties = FALSE) %>% ungroup() %>%
  group_by(species1, species2) %>% mutate(rel_like = exp((min(AIC) - AIC) * 0.5)) %>%
  mutate(weight = rel_like / sum(rel_like)) %>% ungroup()
wt$migration = ifelse(wt$Model %in% c("no_mig", "no_mig_size"), FALSE, TRUE)
wt2 = wt %>% group_by(species1, species2, migration) %>% 
  summarize(sumwt = sum(weight)) %>% ungroup() %>% 
  filter(migration == TRUE) %>% left_join(d %>% dplyr::select(species1, species2, snp_count))
# snp count does not appear to be affecting modelw eight
ggplot(wt2, aes(snp_count, sumwt)) + geom_point() + 
   theme_classic()

# calculate the difference in AIC
# between a model with migration
# and one without
nomig = d %>% group_by(species1, species2) %>% 
  filter(Model %in% c("no_mig", "no_mig_size")) %>%
  summarize(minAIC_no = min(AIC)) %>% ungroup()
mig = d %>% group_by(species1, species2) %>% 
  filter(Model %in% c("asym_mig", "asym_mig_size", "sym_mig", "sym_mig_size")) %>%
  summarize(minAIC_mig = min(AIC)) %>% ungroup()
d1 = d %>% left_join(nomig %>% dplyr::select(species1, species2, minAIC_no)) %>%
  left_join(mig %>% dplyr::select(species1, species2, minAIC_mig))

# include AIC diff
# include wt
d2 = d1 %>% group_by(species1, species2) %>%
  slice_min(AIC, with_ties = FALSE) %>% ungroup() %>%
  mutate(AIC_diff = abs(minAIC_mig - minAIC_no)) %>% 
  left_join(wt2 %>% dplyr::select(species1, species2, sumwt) %>% distinct())
d2$migration = ifelse(d2$Model %in% c("no_mig", "no_mig_size"), FALSE, TRUE)
table(d2$migration)
View(d2[d2$AIC_diff < 2, ])

params = list(
  no_mig = c("n1", "n2", "T"),
  no_mig_size = c('n1_split', 'n2_split', 'n1', 'n2', 'T', 'T_size'),
  asym_mig = c('n1', 'n2', 'm12', 'm21', 'T'),
  sym_mig = c('n1', 'n2', 'm', 'T'),
  asym_mig_size  = c('n1_split', 'n2_split', 'n1', 'n2', 'm12', 'm21', 'T', 'T_size'),
  sym_mig_size  = c('n1_split', 'n2_split', 'n1', 'n2', 'm', 'T', 'T_size')
)

get = c("n1", "n2", "T", "m12", "m21", "m")
res = vector("list", length = nrow(d2))
for (i in 1:nrow(d2)) {
  pp = as.numeric(strsplit(pull(d2[i, "param"]), ",")[[1]])
  names(pp) = params[[pull(d2[i, "Model"])]]
  pp = pp[get]
  names(pp) = get
  pp[!complete.cases(pp)] = 0
  vec1 = d2[i, c("species1", "species2", "Model", "log.likelihood", "AIC", "theta")] %>% unlist(., use.names=FALSE)
  vec2 = c(vec1, pp)
  res[[i]] = vec2
}

res2 = data.frame(do.call("rbind", res))
names(res2) = c("species1", "species2", "model",
                "log.likelihood", "AIC", "theta",
                "n1", "n2", "T", "M12", "M21", "M")
num = c("log.likelihood", "AIC", "theta",
        "n1", "n2", "T", "M12", "M21", "M")
for (i in 1:length(num)) {
  res2[, num[i]] = as.numeric( res2[, num[i]] )
}

res2[res2$model %in% c("sym_mig_size", "sym_mig"), "M12"] = 
  res2[res2$model %in% c("sym_mig_size", "sym_mig"), "M"]
res2[res2$model %in% c("sym_mig_size", "sym_mig"), "M21"] = 
  res2[res2$model %in% c("sym_mig_size", "sym_mig"), "M"]

d3 = res2 %>% left_join(d2 %>% dplyr::select(species1, species2, sumwt, migration, snp_count))

ggplot(d3 %>% filter(snp_count > 100), aes(M12, M21)) +
  geom_point(shape = 16, col = "black", alpha  = 0.3) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_vline(xintercept = 1, linetype = "dotted") +
  geom_point(data = d3 %>% filter(migration == TRUE), fill = "red", shape = 21) +
  theme_classic() + xlab(expression(2 * Nm[12])) +
  ylab(expression(2 * Nm[21]))

x = read.csv("data/abc/all_Xantusia.16Jan24.csv")
xx3 = x %>% filter(metric %in% c("probability", "best_model")) %>% 
   spread(metric, value) %>% 
   mutate_at(c('probability'), as.numeric)  %>%
   rename(best_abc_model_no_sfs = best_model,
          probability_no_sfs = probability)

 
x2 = d3 %>% full_join(xx3) %>% 
   mutate(best_moment_model = ifelse(model %in% c("sym_mig_size", "sym_mig", "asym_mig_size", "asym_mig"),
                                     "migration", "isolation")) %>%
   mutate(p_migration_no_sfs = ifelse(best_abc_model_no_sfs == "migration", probability_no_sfs, 1 - probability_no_sfs)) %>%
   select(-M, -probability_no_sfs) 
 
x2 %>% group_by(best_moment_model, best_abc_model_no_sfs) %>% summarize(n = n())
write.csv(x2, "data/demographic_results.16Jan24.csv",
           row.names = F)
ggplot(x2 %>% filter(snp_count > 200), aes(M12, M21)) +
  geom_point(shape = 16, col = "black", alpha  = 0.3) +  
   geom_hline(yintercept = 1, linetype = "dotted") +
   geom_vline(xintercept = 1, linetype = "dotted") +
   geom_point(data = x2 %>% filter(best_abc_model_no_sfs == "migration"), fill = "red", shape = 21) + 
   theme_classic() + xlab(expression(2 * Nm[12])) +
   ylab(expression(2 * Nm[21]))

x2 %>% filter(complete.cases(M12), complete.cases(M21)) %>% 
   filter(best_moment_model == "migration") %>% 
   summarize(quantile(M21, c(0.01, 0.1, 0.5, 0.9, 0.99)),
             quantile(M12, c(0.01, 0.1, 0.5, 0.9, 0.99)))
x2 %>% filter(M12 > 1 | M21 > 1) %>% group_by(best_abc_model_no_sfs) %>%
  summarize(n = n())

ggplot(x2 %>% filter(snp_count > 100)) +
  geom_histogram(aes(p_migration_no_sfs)) +
  theme_classic()
x2 %>% arrange(desc(p_migration_no_sfs)) %>% View()
