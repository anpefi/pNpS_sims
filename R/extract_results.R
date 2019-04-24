#' Extract results from a oncoSimulIndiv object
#'
#' @description For the different  genes cases
#'
#' @param x An oncoSimulIndiv object
#' @param loci A vector of the selection coefficents for each locus
#' @param opt A list with the options of the simulations
#'
#' @return A data.frame for all the stats computed

extract_results <- function(x, loci, opt, d = 0.05){

  if(class(opt)!="list"){
    stop("opt is not a list")
  }

  if(is.null(x$TotalPopSize)) {
    return(rep(NA, length(x$geneNames)))
  }
  n <- nrow(x$pops.by.time)

  gene_size <- opt$nNS_gene+opt$nS_gene
  total_genes <- opt$n_pos + opt$n_neu + opt$n_neg

  driver <-
    rep(c(rep(TRUE,opt$n_pos),rep(FALSE,opt$n_neg),rep(FALSE,opt$n_neu)),
        rep(gene_size,total_genes))
  neutral <-
    rep(c(rep(FALSE,opt$n_pos),rep(FALSE,opt$n_neg),rep(TRUE,opt$n_neu)),
        rep(gene_size,total_genes))
  deletereous <-
    rep(c(rep(FALSE,opt$n_pos),rep(TRUE,opt$n_neg),rep(FALSE,opt$n_neu)),
        rep(gene_size,total_genes))
  NS_site <- rep(rep(c(TRUE,FALSE),c(opt$nNS_gene,opt$nS_gene)), total_genes)
  S_site <- !NS_site



  tumorGenotype <- apply(x$Genotypes[NS_site & driver,],2, sum) > 0
  result <- NULL
  for (i in 2:n){ #Avoid t=0
    pop <- x$pops.by.time[i, -1]
    popSize <- x$PerSampleStats[i, 1]
    tumorSize <- sum(pop[tumorGenotype])
    n_clones <- sum(pop[tumorGenotype]/tumorSize>0.05) #Avoid undeveloped mutants
    lg_clone <- max(pop[tumorGenotype])

    if(tumorSize < 1 && opt$s_pos > 0) next
    if (tumorSize > 0) {nsize <- tumorSize} else nsize <- popSize

    freqs <- tcrossprod(pop[tumorGenotype], x$Genotypes[,tumorGenotype])/nsize

    variants <- freqs > 0 #True variants, all of them
    sampled <- freqs > d  #sampled variants, those above the detection threshold
    clonal <- freqs == 1 #Fixed variants

    #Number of mutations (variants)
    result <- rbind(result,
                 c( sum(variants[ NS_site & driver])/opt$n_pos,
                    sum(variants[ S_site & driver])/opt$n_pos,
                    sum(variants[ NS_site & neutral])/opt$n_neu,
                    sum(variants[ S_site & neutral])/opt$n_neu,
                    sum(variants[ NS_site & deletereous])/opt$n_neg,
                    sum(variants[ S_site & deletereous])/opt$n_neg,
                    sum(sampled[ NS_site & driver])/opt$n_pos,
                    sum(sampled[ S_site & driver])/opt$n_pos,
                    sum(sampled[ NS_site & neutral])/opt$n_neu,
                    sum(sampled[ S_site & neutral])/opt$n_neu,
                    sum(sampled[ NS_site & deletereous])/opt$n_neg,
                    sum(sampled[ S_site & deletereous])/opt$n_neg,
                    sum(clonal[ NS_site & driver])/opt$n_pos,
                    sum(clonal[ S_site & driver])/opt$n_pos,
                    sum(clonal[ NS_site & neutral])/opt$n_neu,
                    sum(clonal[ S_site & neutral])/opt$n_neu,
                    sum(clonal[ NS_site & deletereous])/opt$n_neg,
                    sum(clonal[ S_site & deletereous])/opt$n_neg,
                    prod(1 + freqs*loci),
                    prod(1 + freqs[driver]*loci[driver]),
                    prod(1 + freqs[deletereous]*loci[deletereous]),
                    prod(1 + freqs[neutral]*loci[neutral]), #This should be 1
                    tumorSize, n_clones, lg_clone/tumorSize,
                    popSize - tumorSize,
                    x$pops.by.time[i,1] )
            )
  }
  out <- list(result=result,freqs=freqs)
  return(out)

}
