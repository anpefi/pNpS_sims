#!/usr/bin/env Rscript

## Version for three genes (driver, deletereous and neutral), with the same number of NS and S sites

### PARSING COMMAND LINE ARGUMENTS ###
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(tidyverse))
devtools::load_all(here::here(),quiet = TRUE)

option_list = list(
  optparse::make_option(c("-i", "--id"), action="store", default="run", type='character',
                help="ID name for this run"),
  optparse::make_option(c("--model"), action="store", default="Exp", type='character',
                help="Model used in this simulation (Exp/McFL)"),
  optparse::make_option(c("-r", "--reps"), action="store", default=100, type='integer',
                        help="This is the number of replicates of cohorts to run"),
  optparse::make_option(c("--inds"), action="store", default=100, type='integer',
                        help="This is the number of individuals in the cohort"),
  optparse::make_option(c("-n", "--initSize"), action="store", default=10000, type='integer',
                help="Initial size for the cell population"),
  optparse::make_option(c("-t", "--finalTime"), action="store", default=3000, type='integer',
                help="Max time units to evolve the cell population"),
  optparse::make_option(c("--detectionSize"), action="store", default=1e8, type='integer',
                help="Threshold population size to stop the simulation"),
  optparse::make_option(c("--sampleEvery"), action="store", default=0.1, type='double',
                help="Time interval to the simulation check the status"),
  optparse::make_option(c("--keepEvery"), action="store", default=3000, type='double',
                help="#Time units to keep samples for further analysis"),
  optparse::make_option(c("--wall.time"), action="store", default=600, type='integer',
                help="maximum wall time"),


  optparse::make_option(c("-m", "--mu"), action="store", default=2.66e-9, type='double',
                help="Mutation rate (per site and cell division)"), #Following Milholland et al 2017
    # Number of loci. N/S ~ 2.76 following Mulholland et al 2017


  optparse::make_option(c("--nNS_gene"), action="store", default=552, type='integer',
                        help="Number of Non-Synonymous mutations in each gene"),
  optparse::make_option(c("--nS_gene"), action="store", default=200, type='integer',
                        help="Number of Synonymous mutations in each gene"),


  optparse::make_option(c("--s_pos"), action="store", default=0.1, type='double',
                help="Selection coefficient of driver mutations"),
  optparse::make_option(c("--s_neg"), action="store", default=-0.1, type='double',
                help="Selection coefficient of driver mutations"),
  optparse::make_option(c("--mutationPropGrowth"), action="store_true", default=TRUE, type='logical',
                help="Is the mutation rate proportional to the pop growth rate?"),
  optparse::make_option(c("--onlyCancer"), action="store_true", default=TRUE, type='logical',
                help="Repeat if cancer not reached"),
  optparse::make_option(c("-c","--mc.cores"), action="store", default=1, type='integer',
                help="number of cores"),
  optparse::make_option(c("--seed"), action="store", default=0, type='integer',
                help="seed for random number generator (0 for time)"),
  optparse::make_option(c("--debug"), action="store_true", default=FALSE, type='logical',
                help="Run with debugging options")
)
opt = optparse::parse_args(optparse::OptionParser(option_list=option_list))
saveRDS(opt, file= paste0(opt$id,".opt.rds"), compress = FALSE)

### SOME SETUPS ###
attach(opt)
if(!debug) options(warn = -1)

# Initalizing some variables
## gene driver
gene_drv <- c( rep(s_pos,nNS_gene), rep(0,nS_gene) )
gene_neg <- c( rep(s_neg,nNS_gene), rep(0,nS_gene) )
gene_neu <- c( rep(0,nNS_gene), rep(0,nS_gene) )
gene_size <- nNS_gene + nS_gene


loci <- c(gene_drv,gene_neg,gene_neu)
names(loci) <- 1:length(loci)
drvNames = names(loci[1:nNS_gene])

# Setting the Fitness effect on the simulator
fe <- OncoSimulR::allFitnessEffects(noIntGenes = loci, drvNames = drvNames )

### SIMULATION
#Reproducibility setting
RNGkind("Mersenne-Twister")
if(seed>0) set.seed(seed)

#Run
result <- NULL

for (r in 1:reps){ ## Change
  pp <- OncoSimulR::oncoSimulPop(inds, fe, model=model, mu= mu, onlyCancer = onlyCancer,
                                 detectionSize = detectionSize, detectionDrivers = NA,
                                 detectionProb = NA, sampleEvery=sampleEvery,
                                 initSize = initSize, finalTime = finalTime,
                                 keepEvery=keepEvery, mutationPropGrowth = mutationPropGrowth,
                                 mc.cores = mc.cores, max.wall.time = wall.time,
                                 errorHitWallTime=FALSE)

  fitness <- sapply(pp,extract_fitness_genes, loci, opt) %>% t %>% colMeans %>%
    t %>% as.data.frame
  colnames(fitness) <- c("w_total", "w_driver","w_deletereous", "w_neutral")

  samples <- OncoSimulR::samplePop(pp,thresholdWhole = 1e-9 )

  genes_cohorts <- as.data.frame(t(samples)) %>%
    split(rep(1:3, each=(nNS_gene+nS_gene)))
  names(genes_cohorts) <- c("driver","deletereous","neutral")
  data <- lapply(genes_cohorts, function(X){calc_pNpS_cohort(t(X),nNS_gene,nS_gene)}) %>%
    as.data.frame
  FinalSize <- do.call("rbind", lapply(pp, "[[", "TotalPopSize")) %>% mean
  FinalTime <- do.call("rbind", lapply(pp, "[[", "FinalTime")) %>% mean
  result <- rbind(result, cbind(r, data, FinalSize, FinalTime, fitness))
}
write_csv(result,paste0(id,".result.csv"))

detach(opt)
