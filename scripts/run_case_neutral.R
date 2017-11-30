#!/usr/bin/env Rscript
library(OncoSimulR)
library(tidyverse)

## Set the parameters

# Number of loci. N/S ~ 2.76 following Mulholland et al 2017
nNS_pos <- 10 
nNS_neg <- 10
nNS_neu <- 27580
nNS <- nNS_pos + nNS_neg + nNS_neu
nS <- 10000
#Effect of NS mutations. All Neutral, independently of their label
sNS_pos <- rep(0,nNS_pos)
sNS_neg <- rep(0,nNS_neg)
sNS_neu <- rep(0,nNS_neu)
#Effect of S mutations. All neutral
sS <- rep(0,nS)
#Mutation rates
mu <- 2.66e-9  #Following Milholland et al 2017
#Set the loci
loci <- c(sNS_pos,sNS_neg,sNS_neu, sS)
names(loci) <- 1:length(loci)
# Setting the Fitnes effect on the simulator
fe <- allFitnessEffects(noIntGenes = loci, drvNames = names(loci[1:(nNS_pos+1)]) )

#Simulation settings
simName <- "neutral"
nreps <- 100 #Number of individuals (replicates)
initSize <- 10000 #Starting population size
finalTime <- 10000 #Maximum time
detectionSize <- 1e9 #Population size threshold to stop the simulation
sampleEvery <- 10 #Time unit to check status of the population
keepEvery <- 100 #Time units to keep samples for further analysis
mutationPropGroth <- FALSE  #Is the mutation rate proportional to the pop growth rate?
mc.cores <- 5 #cores to use

### SIMULATION
#Reproducibility setting
RNGkind("Mersenne-Twister")
set.seed(2134)

#Run 
pp <- oncoSimulPop(nreps, fe, mu= mu, onlyCancer = F, detectionSize = detectionSize, 
                   detectionDrivers = NA, detectionProb = NA, sampleEvery=sampleEvery,
                   initSize = initSize, finalTime = finalTime, keepEvery=keepEvery,
                   mutationPropGrowth = FALSE, mc.cores = mc.cores)

saveRDS(pp,file = paste0(simName,".pop.gz"), compress = "gzip")

### SAMPLING
# Function
NSsample <- function(x, det.limit=0.05, nNS, nS, finalTime, keepEvery){
  # adapted form get.mut.vec in OncoSimulR
  if(is.null(x$TotalPopSize)) {
    return(rep(NA, length(x$geneNames)))
  }
  tsamples <- seq(0, finalTime, keepEvery)
  n <- length(tsamples)
  result <- data.frame(t=tsamples,
                       TS=rep(0,n),
                       N=rep(NA,n),S=rep(NA,n),pN_pS=rep(NA,n))
   for (i in 1:n){
    if (!is.element(x$pops.by.time[i,1], tsamples)) break
  
    pop <- x$pops.by.time[i, -1]
    
    if(all(pop == 0)) {
      #result <- rbind(result, c(i,0,NA,NA,NA))
      break
    }
    
    popSize <- x$PerSampleStats[i, 1]
    variants <- as.numeric((tcrossprod(pop, x$Genotypes)/popSize) >= det.limit)
    N <-sum(variants[1:nNS])
    S <-sum(variants[nNS+1:nS])
    if (S>0) pNpS <- (N / nNS) / (S / nS)
    else pNpS <- NA
    result[i,] <- c(tsamples[i],popSize,N,S,pNpS)
  }
  
  return(result)
  
}

system.time(samples <- do.call(rbind, Map(NSsample, pp, det.limit=0.05, nNS=nNS, nS=nS, finalTime=10000, keepEvery=100)) )

saveRDS(as.tibble(samples), file=paste0(simName,".samples.gz"), compress = "gzip")

