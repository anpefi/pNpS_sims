# Consolidate samples

library(tidyverse, quietly = T, warn.conflicts = FALSE)

samples <- do.call("rbind",
  lapply(list.files(here::here("results"),
  pattern = "sel_0[0-5]_Exp_[0-9]*.samples.rds",
  full.names = TRUE),
  function(X){readRDS(X) %>%
              mutate(case=X %>% stringr::str_split("_") %>%
              unlist %>% .[3] %>% as.factor,
              rep = X %>%  stringr::str_extract("Exp.*samples") %>%
                stringr::str_sub(5,-9) %>% as.numeric )}))
levels(samples$case) <- c("neutral", "symmetric", "10x_negative", "mostly_negative", "all_positive", "all_negative")

saveRDS(samples, "results/sel_All_Exp.data.rds")
