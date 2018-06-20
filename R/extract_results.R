#' Extract fitness of samples from a oncoSimulIndiv object
#'
#' @description For three genes cases
#'
#' @param x An oncoSimulIndiv object
#' @param loci A vector of the selection coefficents for each locus
#' @param opt A list with the options of the simulations
#'   The needed elements are:
#'   nNS_gene
#'   nS_gene
#'   s_pos
#'   s_neg
#' @return A vector of fitness for all the tumor and individual genes

extract_results <- function(x, loci, opt, d = 0.05){

  if(class(opt)!="list"){
    stop("opt is not a list")
  }

  if(is.null(x$TotalPopSize)) {
    return(rep(NA, length(x$geneNames)))
  }
  n <- nrow(x$pops.by.time)

  gene_size <- opt$nNS_gene+opt$nS_gene
  driver <- rep(c(TRUE,FALSE,FALSE),rep(gene_size,3))
  neutral <- rep(c(FALSE,FALSE,TRUE),rep(gene_size,3))
  deletereous <- rep(c(FALSE,TRUE,FALSE),rep(gene_size,3))
  NS_site <- rep(rep(c(TRUE,FALSE),c(opt$nNS_gene,opt$nS_gene)), 3)
  S_site <- !NS_site



  tumorGenotype <- apply(x$Genotypes[1:opt$nNS_gene,],2, sum) > 0
  result <- NULL
  for (i in 2:n){ #Avoid t=0
    pop <- x$pops.by.time[i, -1]
    popSize <- x$PerSampleStats[i, 1]
    tumorSize <- sum(pop[tumorGenotype])

    if(tumorSize < 1 && opt$s_pos > 0) next
    if (tumorSize > 0) nsize <- tumorSize
    else nsize <- popSize

    freqs <- tcrossprod(pop[tumorGenotype], x$Genotypes[,tumorGenotype])/nsize

    variants <- freqs > 0 #True variants, all of them
    sampled <- freqs > d  #sampled variants, those above the detection threshold
    clonal <- freqs > 0.95 #Fixed (with some error) variants

    #Number of mutations (variants)
    result <- rbind(result,
                 c( sum(variants[ NS_site & driver]),
                 sum(variants[ S_site & driver]),
                 sum(variants[ NS_site & neutral]),
                 sum(variants[ S_site & neutral]),
                 sum(variants[ NS_site & deletereous]),
                 sum(variants[ S_site & deletereous]),
                 sum(sampled[ NS_site & driver]),
                 sum(sampled[ S_site & driver]),
                 sum(sampled[ NS_site & neutral]),
                 sum(sampled[ S_site & neutral]),
                 sum(sampled[ NS_site & deletereous]),
                 sum(sampled[ S_site & deletereous]),
                 sum(clonal[ NS_site & driver]),
                 sum(clonal[ S_site & driver]),
                 sum(clonal[ NS_site & neutral]),
                 sum(clonal[ S_site & neutral]),
                 sum(clonal[ NS_site & deletereous]),
                 sum(clonal[ S_site & deletereous]),
                 prod(1 + freqs*loci),
                 prod(1 + freqs[driver]*loci[driver]),
                 prod(1 + freqs[deletereous]*loci[deletereous]),
                 prod(1 + freqs[neutral]*loci[neutral]), #This should be 1
                 tumorSize,
                 popSize - tumorSize,
                 x$pops.by.time[i,1] )
            )
  }
  return(result)

}
