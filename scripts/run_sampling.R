#!/usr/bin/env Rscript

### PARSING COMMAND LINE ARGUMENTS ###
suppressPackageStartupMessages(require(optparse))
option_list = list(
    make_option(c("-i", "--id"), action="store", default="run", type='character',
                help="ID name for the sim to be read"),

    make_option(c("-d", "--det.limit"), action="store", default=1e-9, type='double',
                help="Frequency threshold to detect a variant")

)
opt <- parse_args(OptionParser(option_list=option_list))
opt_sim <- readRDS(file = paste0(opt$id,".opt.rds"))
# In case they are argument shared, this call will update the simulation ones
opt <- c( opt_sim[!is.element(names(opt_sim), names(opt))] , opt )

### LOADING REQUIRED LIBRARIES ###
suppressPackageStartupMessages(library(OncoSimulR))
suppressPackageStartupMessages(library(tidyverse))

### LOADING SIMULATION DATA
attach(opt)
pp <- readRDS(file = paste0(id,".pop.rds"))

# Initalizing some variables
nNS <- nNS_pos + nNS_neg + nNS_neu

### SAMPLING
# Function
NSsample <- function(x, det.limit=0.05, nNS, nS, finalTime, keepEvery){
    # adapted form get.mut.vec in OncoSimulR
    if(is.null(x$TotalPopSize)) {
        return(rep(NA, length(x$geneNames)))
    }
    tsamples <- x$pops.by.time[,1]
    n <- length(tsamples)
    result <- tibble(t=tsamples,
                     TS=rep(0,n),
                     N=rep(NA,n),S=rep(NA,n), nDrivers=rep(0,n),
                     pN_pS=rep(NA,n))
    for (i in 1:n){
        #if (!is.element(x$pops.by.time[i,1], tsamples)) break

        pop <- x$pops.by.time[i, -1]

        if(all(pop == 0)) {
            #result <- rbind(result, c(i,0,NA,NA,NA))
            break
        }

        popSize <- x$PerSampleStats[i, 1]
        variants <- as.numeric((tcrossprod(pop, x$Genotypes)/popSize) >= det.limit)
        ndrivers <- sum(variants[1:nNS_pos])
        N <-sum(variants[1:nNS])
        S <-sum(variants[nNS+1:nS])
        if (S>0) pNpS <- (N / nNS) / (S / nS)
        else pNpS <- NA
        result[i,] <- c(tsamples[i],popSize,N,S,ndrivers,pNpS)
    }

    return(result)

}

samples <- do.call(rbind, Map(NSsample, pp, det.limit=det.limit, nNS=nNS, nS=nS, finalTime=finalTime, keepEvery=keepEvery))

saveRDS(as.tibble(samples), file=paste0(id,".samples.",det.limit,".rds"), compress = FALSE)
detach(opt)