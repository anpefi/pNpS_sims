#' Extract fitness of samples from a oncoSimulIndiv object
#'
#' @description For three genes cases
#'
#' @param x An oncoSimulIndiv object
#' @param opt A list with the options of the simulations
#'   The needed elements are:
#'   nNS_gene
#'   nS_gene
#'   s_pos
#'   s_neg
#' @return A vector of fitness for all the tumor and individual genes

extract_fitness_genes <- function(x, loci, opt, d = 1e-09){

  if(class(opt)!="list"){
    stop("opt is not a list")
  }

  if(is.null(x$TotalPopSize)) {
    return(rep(NA, length(x$geneNames)))
  }
  n <- nrow(x$pops.by.time)
  tsamples <- x$pops.by.time[n,1]
  pop <- x$pops.by.time[n, -1]

  if(all(pop == 0)) {
      z <- z + 1
      break
  }

  popSize <- x$PerSampleStats[n, 1]
  freqs <- tcrossprod(pop, x$Genotypes)/popSize
  gene_size <- opt$nNS_gene+opt$nS_gene

  #Individual fitness
  w_total <- prod(1 + freqs*loci)
  #Gene_fitness
  w_drv <- prod(1 + freqs[loci>0]*loci[loci>0])
  w_del <- prod(1 + freqs[loci<0]*loci[loci<0])
  w_neu <- prod(1 + freqs[loci==0]*loci[loci==0]) #This should be 1

  return(c(w_total,w_drv,w_del,w_neu))

}
