# Consolidate samples

library(tidyverse, quietly = T, warn.conflicts = FALSE)

samples <- do.call("rbind", lapply(list.files("results/", pattern= "neu_McFL_[0-9]*.samples.1e-09.rds", full.names = TRUE), readRDS)) %>% filter(!is.na(t)) %>% mutate(threshold=1e-09, model="McFL")

samples <- rbind(samples, do.call("rbind", lapply(list.files("results/", pattern= "neu_McFL_[0-9]*.samples.0.001.rds", full.names = TRUE), readRDS)) %>% filter(!is.na(t)) %>% mutate(threshold=0.001, model="McFL"))

samples <- rbind(samples, do.call("rbind", lapply(list.files("results/", pattern= "neu_McFL_[0-9]*.samples.0.01.rds", full.names = TRUE), readRDS)) %>% filter(!is.na(t)) %>% mutate(threshold=0.01, model="McFL"))

samples <- rbind(samples, do.call("rbind", lapply(list.files("results/", pattern= "neu_McFL_[0-9]*.samples.0.05.rds", full.names = TRUE), readRDS)) %>% filter(!is.na(t)) %>% mutate(threshold=0.05, model="McFL"))

r <- samples %>% filter(t==0, threshold==1e-09) %>% count %>% as.numeric
times <- samples %>% filter(threshold==1e-09) %>% count %>% as.numeric / r

samples$rep <- rep(rep(1:r, each=times),4)

saveRDS(samples, "results/neu_McFL.data.rds")
