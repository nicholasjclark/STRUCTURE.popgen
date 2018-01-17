#'Compute observed heterozygosities and richness metrics for STRUCTURE clusters
#'
#'This function randomly assigns individuals to clusters using STRUCTURE
#'assignment probabilities provided in the qmatrix. It then uses resampling
#'to generate distributions of estimated heterozygosities and allelic richness metrics
#'per locus, per cluster.
#'
#'@importFrom magrittr %>%
#'
#'@param qmatrix Dataframe. The STRUCTURE qmatrix, with individuals in column 'ind'
#'@param allele.dat Dataframe. Raw allele data (sep = '/')
#'@param nclusters Integer. The number of clusters identified by STRUCTURE
#'@param nreps Integer. The number of times to repeat estimations. Default is 100
#'@param samp.snps Logical. If TRUE, a specified number of SNPs will be randomly sampled
#'in each iteration to account for uncertainty and speed computations. Default is FALSE
#'@param ncores (optional) Integer. Number of cores to split job across for parallel
#'computing. Default is 1 (no parallelisation)
#'@param nsamp (optional) Integer. Number of SNPs to sample during each iteration. If
#'samp.snps = TRUE and nsamp is not provided, default is to sample 50% of the loci
#'(i.e. ncol(allele.dat)/2) in each iteration
#'@param NA.symbol The code used to represent missing allele data
#'@return Matrix of mean and 95% quantiles of estimated allele frequencies per cluster
#'@export
obs.metrics.STRUCTURE.popgen = function(qmatrix,allele.dat,nclusters,nreps,
                                     NA.symbol,samp.snps,nsamp,ncores){

  #Specify default values for optional arguments
  if(missing(nreps)) {
    nreps = 100
  }

  if(missing(nsamp)){
    nsamp = ncol(allele.dat)/2
  }

  if(missing(ncores)){
    ncores = 1
  }

  if(missing(samp.snps)){
    samp.snps = FALSE
  }

  if(samp.snps==TRUE){
    #Initiate the parallel clusters
    cl <- makePSOCKcluster(ncores)
    setDefaultCluster(cl)

    #Export necessary data and variables to each cluster
    clusterExport(NULL, c('nreps','allele.dat','qmatrix','nclusters',
                          'NA.symbol','nsamp'),envir=environment())

    #Export necessary functions to each cluster
    clusterExport(NULL, c('STRUCTURE.popgen'))

    #Export necessary libraries
    clusterEvalQ(cl, library(adegenet))
    clusterEvalQ(cl, library(PopGenReport))
    clusterEvalQ(cl, library(hierfstat))

    test <- parLapply(NULL, seq(1:nreps), function(x) {
      STRUCTURE.popgen(qmatrix,allele.dat,nclusters,stat="allele.freqs",NA.symbol,
                       samp.snps=T,nsamp=nsamp)

    })
    raw.freqs = plyr::rbind.fill(test)
    stopCluster(cl)
  }

  if(samp.snps==FALSE){
    cl <- makePSOCKcluster(ncores)
    setDefaultCluster(cl)

    clusterExport(NULL, c('nreps','allele.dat','qmatrix','nclusters',
                          'NA.symbol'),envir=environment())

    clusterExport(NULL, c('STRUCTURE.popgen'))

    clusterEvalQ(cl, library(adegenet))
    clusterEvalQ(cl, library(PopGenReport))
    clusterEvalQ(cl, library(hierfstat))

    test <- parLapply(NULL, seq(1:nreps), function(x) {
      STRUCTURE.popgen(qmatrix,allele.dat,nclusters,stat="allele.freqs",NA.symbol,
                       samp.snps=F)
    })
    raw.freqs <- plyr::rbind.fill(rvest::pluck(test,'obs.richness'))
    raw.Ho <- plyr::rbind.fill(rvest::pluck(test,'Ho.stats'))
    raw.Hs <- plyr::rbind.fill(rvest::pluck(test,'Hs.stats'))

    stopCluster(cl)
  }

  summed.freqs = raw.freqs %>%
    dplyr::mutate_all(dplyr::funs(replace(., which(.==0), NA))) %>%
    dplyr::group_by(locus) %>%
    tidyr::gather(key=cluster,value=frequency,-locus)

  summed.Ho = raw.Ho %>%
    dplyr::mutate_all(dplyr::funs(replace(., which(.==0), NA))) %>%
    dplyr::group_by(locus) %>%
    tidyr::gather(key=cluster,value=Ho,-locus)

  summed.Hs = raw.Hs %>%
    dplyr::mutate_all(dplyr::funs(replace(., which(.==0), NA))) %>%
    dplyr::group_by(locus) %>%
    tidyr::gather(key=cluster,value=Hs,-locus)

  return(list(obs.richness = summed.freqs,obs.Ho = summed.Ho,obs.Hs=summed.Hs))
}
