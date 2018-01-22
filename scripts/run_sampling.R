#!/usr/bin/env Rscript

### PARSING COMMAND LINE ARGUMENTS ###
option_list = list(
    optparse::make_option(c("-i", "--id"), action="store", default="run", type='character',
                help="ID name for the sim to be read"),
    optparse::make_option(c("--dir"), action="store", default="scratch", type='character',
                help="directory for input/output"),

    optparse::make_option(c("-d", "--det.limit"), action="store", default=1e-9, type='double',
                help="Frequency threshold to detect a variant")

)
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))
opt_sim <- readRDS(file = here::here(opt$dir,paste0(opt$id,".opt.rds")))
# In case they are argument shared, this call will update the simulation ones
opt <- c( opt_sim[!is.element(names(opt_sim), names(opt))] , opt )
opt$det.limit <- c(1e-09, 0.05)

### LOADING SIMULATION DATA
pp <- readRDS(file = here::here(opt$dir,paste0(opt$id,".pop.rds")))
devtools::load_all(here::here(),quiet = TRUE)
samples <- do.call(rbind, Map(extract_samples, pp, list(opt)))

saveRDS(tibble::as.tibble(samples), file = here::here(opt$dir,paste0(opt$id,".samples.rds")), compress = FALSE)

