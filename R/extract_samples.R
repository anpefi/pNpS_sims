#' Extract samples from a oncoSimulIndiv object
#'
#' @param x An oncoSimulIndiv object
#' @param det.limit A numeric value for the frquency threshold of a detected mutation
#' @param nNS Number of NS (Non-Synonymous) loci
#' @param nS Number of NS (Synonymous) loci
#' @return A tibble with different papameters for all the sampled times

extract_samples <- function(x, det.limit=0.05, nNS, nS, finalTime, keepEvery){
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
