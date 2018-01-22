#!/usr/bin/env Rscript

### PARSING COMMAND LINE ARGUMENTS ###
suppressPackageStartupMessages(require(optparse))
option_list = list(
  optparse::make_option(c("-i", "--id"), action="store", default="run", type='character',
                help="ID name for this run"),
  optparse::make_option(c("--model"), action="store", default="Exp", type='character',
                help="Model used in this simulation (Exp/McFL)"),
  optparse::make_option(c("-r", "--reps"), action="store", default=10, type='integer',
                help="This is the number of individuals/replicates to tub"),
  optparse::make_option(c("-n", "--initSize"), action="store", default=10000, type='integer',
                help="Initial size for the cell population"),
  optparse::make_option(c("-t", "--finalTime"), action="store", default=10000, type='integer',
                help="Max time units to evolve the cell population"),
  optparse::make_option(c("--detectionSize"), action="store", default=1e9, type='integer',
                help="Threshold population size to stop the simulation"),
  optparse::make_option(c("--sampleEvery"), action="store", default=100, type='double',
                help="Time interval to the simulation check the status"),
  optparse::make_option(c("--keepEvery"), action="store", default=100, type='double',
                help="#Time units to keep samples for further analysis"),
  optparse::make_option(c("--wall.time"), action="store", default=600, type='integer',
                help="maximum wall time"),


  optparse::make_option(c("-m", "--mu"), action="store", default=2.66e-9, type='double',
                help="Mutation rate (per site and cell division)"), #Following Milholland et al 2017
    # Number of loci. N/S ~ 2.76 following Mulholland et al 2017
  optparse::make_option(c("--nNS_pos"), action="store", default=10, type='integer',
                help="Number of Non-Synonymous mutations with positive effect (Drivers)"),
  optparse::make_option(c("--nNS_neg"), action="store", default=10, type='integer',
                help="Number of Non-Synonymous mutations with negative effect"),
  optparse::make_option(c("--nNS_neu"), action="store", default=27580, type='integer',
                help="Number of Non-Synonymous mutations with no effect"),
  optparse::make_option(c("--nS"), action="store", default=10000, type='integer',
                help="Number of Synonymous mutations"),
  optparse::make_option(c("--s_pos"), action="store", default=0.1, type='double',
                help="Selection coefficient of driver mutations"),
  optparse::make_option(c("--s_neg"), action="store", default=-0.1, type='double',
                help="Selection coefficient of driver mutations"),
  optparse::make_option(c("--mutationPropGrowth"), action="store_true", default=TRUE, type='logical',
                help="Is the mutation rate proportional to the pop growth rate?"),
  optparse::make_option(c("--onlyCancer"), action="store_true", default=FALSE, type='logical',
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
nNS <- nNS_pos + nNS_neg + nNS_neu
sNS_pos <- rep(s_pos,nNS_pos)
sNS_neg <- rep(s_neg,nNS_neg)
sNS_neu <- rep(0,nNS_neu)
sS <- rep(0,nS)
loci <- c(sNS_pos,sNS_neg,sNS_neu, sS)
names(loci) <- 1:length(loci)
drvNames = if (nNS_pos > 0) names(loci[1:nNS_pos]) else NULL

# Setting the Fitness effect on the simulator
fe <- OncoSimulR::allFitnessEffects(noIntGenes = loci, drvNames = drvNames )

### SIMULATION
#Reproducibility setting
RNGkind("Mersenne-Twister")
if(seed>0) set.seed(seed)

#Run
pp <- OncoSimulR::oncoSimulPop(reps, fe, model=model, mu= mu, onlyCancer = onlyCancer, detectionSize = detectionSize,
                   detectionDrivers = NA, detectionProb = NA, sampleEvery=sampleEvery,
                   initSize = initSize, finalTime = finalTime, keepEvery=keepEvery,
                   mutationPropGrowth = mutationPropGrowth, mc.cores = mc.cores, max.wall.time = wall.time,
                   errorHitWallTime=FALSE)

saveRDS(pp,file = paste0(id,".pop.rds"), compress = "gzip")

detach(opt)
