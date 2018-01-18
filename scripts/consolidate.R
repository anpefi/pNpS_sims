# Consolidate samples

library(tidyverse, quietly = T, warn.conflicts = FALSE)

samples03 <- do.call("rbind", lapply(list.files("results/", pattern = "sel_03_McFL_[0-9]*.samples.rds", full.names = TRUE),readRDS))

samples03 <- samples03 %>% filter(d==1e-09)
r <- samples03 %>% filter(t==0, d==1e-09) %>% count %>% as.numeric
times <- samples03 %>% filter(d==1e-09) %>% count %>% as.numeric / r


saveRDS(samples, "results/sel_03_McFL.data.rds")

samples <- do.call("rbind", lapply(list.files("results/", pattern = "sel_01_McFl_[5-9]*.samples.rds", full.names = TRUE),readRDS))
