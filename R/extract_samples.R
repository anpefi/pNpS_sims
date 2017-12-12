#' Extract samples from a oncoSimulIndiv object
#'
#' @description This function was adapted from OncoSimulR::get.mut.vec
#'
#' @param x An oncoSimulIndiv object
#' @param opt A list with the options of the simulations
#'   The needed elements are:
#'   det.limit
#'   nNS_pos
#'   nNS_neg
#'   nNS_neu
#'   nS
#'   sNS_pos
#'   sNS_neg
#' @return A tibble with different papameters for all the sampled times

extract_samples <- function(x, opt){

  if(class(opt)!="list"){
    stop("opt is not a list")
  }

  if(is.null(x$TotalPopSize)) {
    return(rep(NA, length(x$geneNames)))
  }

  tsamples <- x$pops.by.time[,1]
  n <- length(tsamples)
  m <- n * length(opt$det.limit)
  result <- tibble::tibble(t=rep(tsamples, each=length(opt$det.limit)),
                   TS=rep(0,m),
                   fitness=rep(NA,m),
                   d=rep(opt$det.limit, n),
                   N=rep(NA,m),S=rep(NA,m), drivers=rep(NA,m),
                   deletereous=rep(NA,m),
                   pN_pS=rep(NA,m))
  z <- 1
  for (i in 1:n){
    #if (!is.element(x$pops.by.time[i,1], tsamples)) break

    pop <- x$pops.by.time[i, -1]

    if(all(pop == 0)) {
      z <- z + 1
      break
    }
    nNS <- opt$nNS_pos+opt$nNS_neg+opt$nNS_neu
    popSize <- x$PerSampleStats[i, 1]
    freqs <- tcrossprod(pop, x$Genotypes)/popSize
    effects <- c(rep(opt$s_pos,opt$nNS_pos), rep(opt$s_neg,opt$nNS_neg))
    w <- prod(1 + freqs[1:opt$nNS_pos+opt$nNS_neg]*effects)


    for(d in opt$det.limit){
      variants <- as.numeric(freqs > d)
      drivers <- sum(variants[1:opt$nNS_pos])
      deletereous <- sum(variants[1:opt$nNS_neg])
      N <-sum(variants[1:nNS])
      S <-sum(variants[nNS+1:opt$nS])
      if (S>0) pNpS <- (N / nNS) / (S / opt$nS)
      else pNpS <- NA
      result[z,] <- c(tsamples[i],popSize,w,d,N,S,drivers,deletereous,pNpS)
      z <- z + 1
    }

  }

  return(result)

}
