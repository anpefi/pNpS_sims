#' Calculate the pN/pS across individuals (cohorts)
#'
#' @description Considering a gene (a set of variants dividend in non-synonymous and synonymous sites) for diferent samples/individuals this function returns the pN/pS vale within the gene across the samples
#'
#' @param variants numeric or logic vector of the variants of the considered
#'  gene. Each position represent a site, 0 wild, 1 mutated. Non-synonimous sites should be before the synonymous
#' @param nNS number of synonymous sites
#' @param nS number of synonymous sites
#' @return The value of the pN/pS for the analysed gene

calc_pNpS_cohort <- function(variants, nNS, nS){
  nsamples <- dim(variants)[1]
  N <-sum(variants[,1:nNS])
  S <-sum(variants[,nNS+1:nS])
  if (S>0) pNpS <- (N / (nNS)) / (S / (nS))
  else pNpS <- NA
  return(pNpS)
}
