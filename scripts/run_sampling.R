#!/usr/bin/env Rscript

### PARSING COMMAND LINE ARGUMENTS ###
suppressPackageStartupMessages(require(optparse))
option_list = list(
    make_option(c("-i", "--id"), action="store", default="run", type='character',
                help="ID name for this run"),

    make_option(c("-t", "--finalTime"), action="store", default=10000, type='integer',
                help="Max time units to evolve the cell population"),

    make_option(c("--keepEvery"), action="store", default=100, type='double',
                help="#Time units to keep samples for further analysis"),


    make_option(c("-d", "--det.limit"), action="store", default=1e-9, type='double',
                help="Frequency threshold to detect a variant"),
    # Number of loci. N/S ~ 2.76 following Mulholland et al 2017
    make_option(c("--nNS_pos"), action="store", default=10, type='integer',
                help="Number of Non-Synonymous mutations with positive effect (Drivers)"),
    make_option(c("--nNS_neg"), action="store", default=10, type='integer',
                help="Number of Non-Synonymous mutations with negative effect"),
    make_option(c("--nNS_neu"), action="store", default=27580, type='integer',
                help="Number of Non-Synonymous mutations with no effect"),
    make_option(c("--nS"), action="store", default=10000, type='integer',
                help="Number of Synonymous mutations"),

    make_option(c("--debug"), action="store_true", default=FALSE, type='logical',
                help="Run with debugging options")
)
opt = parse_args(OptionParser(option_list=option_list))

### LOADING REQUIRED LIBRARIES ###
suppressPackageStartupMessages(library(OncoSimulR))
suppressPackageStartupMessages(library(tidyverse))

### LOADING SIMULATION DATA
attach(opt)
pp <- readRDS(file = paste0(id,".pop.gz"))

# Initalizing some variables
nNS <- nNS_pos + nNS_neg + nNS_neu

### SAMPLING
# Function
NSsample <- function(x, det.limit=0.05, nNS, nS, finalTime, keepEvery){
    # adapted form get.mut.vec in OncoSimulR
    if(is.null(x$TotalPopSize)) {
        return(rep(NA, length(x$geneNames)))
    }
    tsamples <- seq(0, finalTime, keepEvery)
    n <- length(tsamples)
    result <- tibble(t=tsamples,
                     TS=rep(0,n),
                     N=rep(NA,n),S=rep(NA,n), nDrivers=rep(0,n),
                     pN_pS=rep(NA,n))
    for (i in 1:n){
        if (!is.element(x$pops.by.time[i,1], tsamples)) break

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

saveRDS(as.tibble(samples), file=paste0(id,".samples.gz"))
detach(opt)