library(tidyverse)

data_path <- "results"


files <- dir(data_path, pattern = "sel_21_.*result" )
data <- files %>%
  map(~ read_csv(file.path(data_path, .))) %>%
  reduce(rbind) %>%
  mutate(s_pos = 0.1, s_neg= -0.05)

files <- dir(data_path, pattern = "sel_22_.*result" )
data <- files %>%
  map(~ read_csv(file.path(data_path, .))) %>%
  reduce(rbind) %>%
  mutate(s_pos = 0.1, s_neg= -0.1) %>%
  rbind(data)

files <- dir(data_path, pattern = "sel_23_.*result" )
data <- files %>%
  map(~ read_csv(file.path(data_path, .))) %>%
  reduce(rbind) %>%
  mutate(s_pos = 0.1, s_neg= -0.2) %>%
  rbind(data)

files <- dir(data_path, pattern = "sel_31_.*result" )
data <- files %>%
  map(~ read_csv(file.path(data_path, .))) %>%
  reduce(rbind) %>%
  mutate(s_pos = 0.2, s_neg= -0.05) %>%
  rbind(data)

files <- dir(data_path, pattern = "sel_32_.*result" )
data <- files %>%
  map(~ read_csv(file.path(data_path, .))) %>%
  reduce(rbind) %>%
  mutate(s_pos = 0.2, s_neg= -0.1) %>%
  rbind(data)

files <- dir(data_path, pattern = "sel_33_.*result" )
data <- files %>%
  map(~ read_csv(file.path(data_path, .))) %>%
  reduce(rbind) %>%
  mutate(s_pos = 0.2, s_neg= -0.2) %>%
  rbind(data)

files <- dir(data_path, pattern = "sel_11_.*result" )
data <- files %>%
  map(~ read_csv(file.path(data_path, .))) %>%
  reduce(rbind) %>%
  mutate(s_pos = 0.05, s_neg= -0.05) %>%
  rbind(data)

files <- dir(data_path, pattern = "sel_12_.*result" )
data <- files %>%
  map(~ read_csv(file.path(data_path, .))) %>%
  reduce(rbind) %>%
  mutate(s_pos = 0.05, s_neg= -0.1) %>%
  rbind(data)

files <- dir(data_path, pattern = "sel_13_.*result" )
data <- files %>%
  map(~ read_csv(file.path(data_path, .))) %>%
  reduce(rbind) %>%
  mutate(s_pos = 0.05, s_neg= -0.2) %>%
  rbind(data)

files <- dir(data_path, pattern = "sel_41_.*result" )
data <- files %>%
  map(~ read_csv(file.path(data_path, .))) %>%
  reduce(rbind) %>%
  mutate(s_pos = 0.05, s_neg= -0.05) %>%
  rbind(data)

files <- dir(data_path, pattern = "sel_42_.*result" )
data <- files %>%
  map(~ read_csv(file.path(data_path, .))) %>%
  reduce(rbind) %>%
  mutate(s_pos = 0.05, s_neg= -0.1) %>%
  rbind(data)

files <- dir(data_path, pattern = "sel_43_.*result" )
data <- files %>%
  map(~ read_csv(file.path(data_path, .))) %>%
  reduce(rbind) %>%
  mutate(s_pos = 0.05, s_neg= -0.2) %>%
  rbind(data)

write_csv(data, "results/table.csv")


