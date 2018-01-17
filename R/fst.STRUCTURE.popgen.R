#'Estimate Fst values between population clusters
#'
#'This function assigns individuals to clusters using STRUCTURE
#'assignment probabilities provided in the qmatrix. It then calculates
#'pairwise Fst values between clusters.
#'

#'@importFrom magrittr %>%
#'
#'@param qmatrix Dataframe. The STRUCTURE qmatrix, with individuals in column 'ind'
#'@param allele.dat Dataframe/matrix containing '/' separated allele data,
#'with rownames representing individuals
#'@param nclusters (optional) Integer. The number of clusters identified by STRUCTURE
#'@param samp.snps Logical. If TRUE, a specified number of SNPs will be randomly sampled
#'in each iteration to account for uncertainty and speed computations. Default is FALSE
#'@param nsamp (optional) Integer. The number of SNPs to sample during each iteration
#'@param nreps Integer. The number of times to repeat estimations. Default is 100
#'@param ncores Integer. The number of cores to split job across. Default is 1 (no parallelisation)
#'@param NA.symbol The code used to represent missing allele data
#'@return Matrix of pairwise comparison confidence intervals
#'@export

fst.STRUCTURE.popgen = function(qmatrix,allele.dat,nclusters,
                                nreps,NA.symbol,
                                samp.snps,nsamp,ncores){

  #Specify default values for optional arguments
  if(missing(nreps)) {
    nreps = 100
  }

  if(missing(ncores)){
    ncores = 1
  }

  if(missing(samp.snps) & missing(nsamp)){
    samp.snps = FALSE
  }

  if(missing(samp.snps) & !missing(nsamp)){
    samp.snps = TRUE
  }

  if(missing(nsamp)){
    nsamp = ncol(allele.dat)/2
  }

  if(samp.snps==TRUE){
    #Initiate the parallel clusters
    cl <- makePSOCKcluster(ncores)
    setDefaultCluster(cl)

    #Export necessary data and variables to each cluster
    clusterExport(NULL, c('nreps','allele.dat','qmatrix','nclusters',
                          'NA.symbol','nsamp','samp.snps'),
                  envir=environment())

    #Export necessary functions to each cluster
    clusterExport(NULL, c('STRUCTURE.popgen'))

    #Export necessary libraries
    clusterEvalQ(cl, library(adegenet))
    clusterEvalQ(cl, library(PopGenReport))
    clusterEvalQ(cl, library(hierfstat))

    Fst.data <- parLapply(NULL, seq(1:nreps), function(x) {
      STRUCTURE.popgen(qmatrix=qmatrix,
                       allele.dat=allele.dat,
                       nclusters=nclusters,
                       stat="fst",
                       NA.symbol=NA.symbol,
                       samp.snps=samp.snps,
                       nsamp=nsamp)
    })

    raw.fst = plyr::rbind.fill(Fst.data)
    stopCluster(cl)
  }

  if(samp.snps==F){
  cl <- makePSOCKcluster(ncores)
  setDefaultCluster(cl)

  clusterExport(NULL, c('nreps','allele.dat','qmatrix','nclusters',
                        'NA.symbol'),envir=environment())

  clusterExport(NULL, c('STRUCTURE.popgen'))

  clusterEvalQ(cl, library(adegenet))
  clusterEvalQ(cl, library(PopGenReport))
  clusterEvalQ(cl, library(hierfstat))

  Fst.data <- parLapply(NULL, seq(1:nreps), function(x) {
    STRUCTURE.popgen(qmatrix=qmatrix,
                     allele.dat=allele.dat,
                     nclusters=nclusters,
                     stat="fst",
                     NA.symbol=NA.symbol,
                     samp.snps=FALSE)
  })

  raw.fst = plyr::rbind.fill(Fst.data)
  stopCluster(cl)
}

  summed.fst = raw.fst %>%
    dplyr::group_by(category.comparison)%>%
    dplyr::mutate_all(dplyr::funs(lower=quantile(.,probs=0.025),
                    med=median(.,na.rm=T),
                    upper=quantile(.,probs=0.975))) %>%
    dplyr::ungroup() %>%
    dplyr::select(dplyr::matches('category.comparison|lower|upper|med')) %>%
    .[, names(.)[order(names(.))]] %>% dplyr::distinct()

  return(summed.fst)
}

