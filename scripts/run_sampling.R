#!/usr/bin/env Rscript

### LOADING REQUIRED LIBRARIES ###
suppressPackageStartupMessages(library(OncoSimulR))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(require(optparse))

### PARSING COMMAND LINE ARGUMENTS ###
option_list = list(
    make_option(c("-i", "--id"), action="store", default="run", type='character',
                help="ID name for the sim to be read"),
    make_option(c("--dir"), action="store", default="scratch", type='character',
                help="directory for input/output"),

    make_option(c("-d", "--det.limit"), action="store", default=1e-9, type='double',
                help="Frequency threshold to detect a variant")

)
opt <- parse_args(OptionParser(option_list=option_list))
opt_sim <- readRDS(file = here(opt$dir,paste0(opt$id,".opt.rds")))
# In case they are argument shared, this call will update the simulation ones
opt <- c( opt_sim[!is.element(names(opt_sim), names(opt))] , opt )
opt$det.limit <- c(1e-09, 0.001, 0.05)

### LOADING SIMULATION DATA
pp <- readRDS(file = here(opt$dir,paste0(opt$id,".pop.rds")))
devtools::load_all(here(),quiet = TRUE)
samples <- do.call(rbind, Map(extract_samples, pp, list(opt)))

saveRDS(as.tibble(samples), file=here(opt$dir,paste0(opt$id,".samples.rds")), compress = FALSE)

