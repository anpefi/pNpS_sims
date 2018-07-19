#!/usr/bin/env Rscript

## Version for three genes (driver, deletereous and neutral), with the same number of NS and S sites

### PARSING COMMAND LINE ARGUMENTS ###
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(tidyverse))
devtools::load_all(here::here(),quiet = TRUE)

option_list = list(
  optparse::make_option(c("-i", "--id"), action="store", default="run",
                        type='character',
                help="ID name for this run"),
  optparse::make_option(c("--model"), action="store", default="McFL",
                        type='character',
                help="Model used in this simulation (Exp/McFL)"),
  optparse::make_option(c("-r", "--reps"), action="store", default=1,
                        type='integer',
                        help="This is the number of replicates of cohorts to run"),
  optparse::make_option(c("-n", "--initSize"), action="store", default=20,
                        type='integer',
                help="Initial size for the cell population"),
  optparse::make_option(c("-t", "--finalTime"), action="store", default=50000,
                        type='integer',
                help="Max time units to evolve the cell population"),
  optparse::make_option(c("--detectionSize"), action="store", default=1e7,
                        type='integer',
                help="Threshold population size to stop the simulation"),
  optparse::make_option(c("--sampleEvery"), action="store", default=0.01,
                        type='double',
                help="Time interval to the simulation check the status"),
  optparse::make_option(c("--keepEvery"), action="store", default=50000,
                        type='double',
                        help="#Time units to keep samples for further analysis"),
  optparse::make_option(c("--wall.time"), action="store", default=2000,
                        type='integer',
                help="maximum wall time"),


  optparse::make_option(c("-m", "--mu"), action="store", default=1e-8,
                        type='double',    help="Mutation rate (per site and cell division)"), #Following Milholland et al 2017 2.66e-9
    # Number of loci. N/S ~ 2.76 following Mulholland et al 2017


  optparse::make_option(c("--nNS_gene"), action="store", default=500, type='integer',
                        help="Number of Non-Synonymous mutations in each gene"),
  optparse::make_option(c("--nS_gene"), action="store", default=500, type='integer',
                        help="Number of Synonymous mutations in each gene"),
  optparse::make_option(c("--prop_eff"), action="store", default=0.2,
                        type='double',
                        help="Proportion of NS mutations with effect"),

  optparse::make_option(c("--n_pos"), action="store", default=1, type='integer',
                        help="Number of driver (with s_pos) genes"),
  optparse::make_option(c("--n_neu"), action="store", default=1, type='integer',
                        help="Number of neutral genes"),
  optparse::make_option(c("--n_neg"), action="store", default=1, type='integer',
                        help="Number of deletereous (with s_neg) genes"),

  optparse::make_option(c("--s_pos"), action="store", default=0.5, type='double',
                help="Selection coefficient of driver mutations"),
  optparse::make_option(c("--s_neg"), action="store", default=0, type='double',
                help="Selection coefficient of driver mutations"),


  optparse::make_option(c("--mutationPropGrowth"), action="store_true",
                        default=TRUE, type='logical',
                help="Is the mutation rate proportional to the pop growth rate?"),
  optparse::make_option(c("--onlyCancer"), action="store_true", default=TRUE,
                        type='logical',
                help="Repeat if cancer not reached"),
  optparse::make_option(c("--seed"), action="store", default=0, type='integer',
                help="seed for random number generator (0 for time)"),
  optparse::make_option(c("--debug"), action="store_true", default=FALSE,
                        type='logical',
                help="Run with debugging options")
)
opt = optparse::parse_args(optparse::OptionParser(option_list=option_list))
saveRDS(opt, file= paste0(opt$id,".opt.rds"), compress = FALSE)

### SOME SETUPS ###
attach(opt)
if(!debug) options(warn = -1)

# Initalizing some variables
## gene driver
eff_pos <- floor(nNS_gene*prop_eff)
gene_drv <-
  rep(c( rep(s_pos,floor(nNS_gene*prop_eff)),
         rep(0,nNS_gene-floor(nNS_gene*prop_eff)),
         rep(0,nS_gene) ), n_pos)
gene_neg <-
  rep(c( rep(s_neg,floor(nNS_gene*prop_eff)),
         rep(0,nNS_gene-floor(nNS_gene*prop_eff)),
         rep(0,nS_gene) ), n_neg)
gene_neu <-
  rep(c( rep(0,nNS_gene), rep(0,nS_gene) ), n_neu)
gene_size <- nNS_gene + nS_gene


loci <- c(gene_drv,gene_neg,gene_neu)
names(loci) <- 1:length(loci)
drvNames = names(loci[1:nNS_gene*prop_eff])

# Setting the Fitness effect on the simulator
fe <- OncoSimulR::allFitnessEffects(noIntGenes = loci, drvNames = drvNames )

### SIMULATION
#Reproducibility setting
RNGkind("Mersenne-Twister")
if(seed>0) set.seed(seed)

#Run
headers <- c("N_driver","S_driver","N_neutral","S_neutral","N_deletereous","S_deletereous",
             "sp_N_driver","sp_S_driver","sp_N_neutral","sp_S_neutral","sp_N_deletereous",
             "sp_S_deletereous","cl_N_driver","cl_S_driver","cl_N_neutral","cl_S_neutral",
             "cl_N_deletereous","cl_S_deletereous","w_total","w_driver","w_deletereous",
             "w_neutral","tumor_size","healthy_size","Time","id","rep","pNpS_driver",
             "pNpS_neutral","pNpS_deletereous","pNpS_sp_driver","pNpS_sp_neutral",
             "pNpS_sp_deletereous","dNdS_driver","dNdS_neutral","dNdS_deletereous",
             "s_pos","s_neg")
write_csv(as.data.frame(t(headers)), paste0(id,".result.csv"), col_names = FALSE)

for (r in 1:reps){ ## Change
    pp <- OncoSimulR::oncoSimulIndiv(fe, model=model, mu= mu,
                                     onlyCancer = onlyCancer,
                                     detectionSize = detectionSize,
                                     detectionDrivers = NA,
                                     detectionProb = NA, sampleEvery=sampleEvery,
                                     initSize = initSize, finalTime = finalTime,
                                     keepEvery=keepEvery,
                                     mutationPropGrowth = mutationPropGrowth,
                                     max.wall.time = wall.time,
                                     errorHitWallTime=FALSE)
    if(is.null(pp$TotalPopSize)){
      cat("Failed replicate ",r, "\n")

    } else {
      result <- extract_results(pp,loci, opt)
      colnames(result) <- headers[1:25]

      result <- result %>% as.data.frame %>% mutate(
        id,
        r,
        ((N_driver/nNS_gene)/(S_driver/nS_gene)),
        ((N_neutral/nNS_gene)/(S_neutral/nS_gene)),
        ((N_deletereous/nNS_gene)/(S_deletereous/nS_gene)),
        ((sp_N_driver/nNS_gene)/(sp_S_driver/nS_gene)),
        ((sp_N_neutral/nNS_gene)/(sp_S_neutral/nS_gene)),
        ((sp_N_deletereous/nNS_gene)/(sp_S_deletereous/nS_gene)),
        ((cl_N_driver/nNS_gene)/(cl_S_driver/nS_gene)),
        ((cl_N_neutral/nNS_gene)/(cl_S_neutral/nS_gene)),
        ((cl_N_deletereous/nNS_gene)/(cl_S_deletereous/nS_gene)),
        s_pos, s_neg
      )

      write_csv(result, paste0(id,".result.csv"), col_names = FALSE, append = TRUE)
    }
}


detach(opt)
