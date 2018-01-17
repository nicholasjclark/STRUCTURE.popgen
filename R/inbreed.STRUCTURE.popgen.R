#'Subsample SNPs and estimate individual inbreeding coefficients
#'
#'This function assigns individuals to known (i.e. provided) populations
#'or to clusters using STRUCTURE assignment probabilities provided in the qmatrix.
#'It then calculates inbreeding coefficients for individuals and compares
#'cluster-level inbreeding with an ANOVA.
#'
#'@importFrom magrittr %>%
#'
#'
#'@param allele.dat Dataframe/matrix containing '/' separated allele data,
#'with rownames representing individuals
#'@param ind.populations (optional) Matrix. Population assignments for individuals,
#'with rownames representing individuals
#'@param qmatrix.prob (optional) Logical. If TRUE, a supplied STRUCTURE qmatrix
#'will be used to randomly assign individuals to clusters based on assignment probabilities.
#'Default is FALSE
#'@param samp.snps Logical. If TRUE, a specified number of SNPs will be randomly sampled
#'in each iteration to account for uncertainty and speed computations. Default is FALSE
#'@param qmatrix (optional) The STRUCTURE qmatrix, with individuals in column 'ind'. To be used
#'if qmatrix.prob = TRUE
#'@param nclusters (optional) Integer. The number of clusters identified by STRUCTURE
#'@param ncores (optional) Integer. The number of cores to split job across.
#'Default is 1 (no parallelisation)
#'@param nsamp (optional) Integer. The number of SNPs to sample during each iteration
#'@param nreps Integer. The number of times to repeat the subsampling and
#'estimate inbreeding. Default is 100
#'@param NA.symbol The code used to represent missing allele data
#'@return A list of two elements:\cr\cr
#'lin.coefs:   Tukey's population comparisons from repeated ANOVAs comparing individual
#'inbreeding estimates across populations.\cr\cr
#'raw.vals:   summaries of the raw inbreeding estimates for each individual
#'@export

inbreed.STRUCTURE.popgen = function(allele.dat,qmatrix.prob,
                                    samp.snps,nreps,ncores,nsamp,
                                    NA.symbol,qmatrix=NULL,
                                    ind.populations=NULL,
                                    nclusters=NULL){

  #Must supply either a qmatrix OR matrix of population assignment data
  if(is.null(ind.populations) & is.null(qmatrix)){
    stop('assignment data missing: please supply either a STRUCTURE qmatrix OR ind.populations assignment data')
  }

  if(!is.null(ind.populations) & !is.null(qmatrix)){
    stop('doubled assignment data: please supply either a STRUCTURE qmatrix OR ind.populations assignment data')
  }

  #Specify default values for arguments and state warnings
  if(missing(nreps)) {
    nreps = 100
  }

  if(missing(qmatrix.prob) & !is.null(qmatrix)) {
    warning('qmatrix supplied but qmatrix.prob not specified: setting to TRUE')
    qmatrix.prob = TRUE
  }

  if(missing(qmatrix.prob) & !is.null(ind.populations)) {
    qmatrix.prob = FALSE
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

  if(!is.null(qmatrix)){
    qmatrix = qmatrix
    ind.populations = qmatrix
  }

  if(!is.null(ind.populations)){
    ind.populations = ind.populations
  }

  if(!is.null(nclusters)){
    nclusters = nclusters
  }

n.allele.total = ncol(allele.dat)

#Initiate the parallel clusters
cl <- makePSOCKcluster(ncores)
setDefaultCluster(cl)

#Export necessary data and variables to each cluster
clusterExport(NULL, c('nreps','allele.dat','qmatrix','qmatrix.prob',
                      'nclusters','NA.symbol','nsamp','samp.snps','n.allele.total',
                      'ind.populations'), envir=environment())


#Export necessary libraries
clusterEvalQ(cl, library(adegenet))
clusterEvalQ(cl, library(PopGenReport))
clusterEvalQ(cl, library(hierfstat))

inbreeding.ML <- parLapply(NULL, seq(1:nreps), function(x) {
  if(samp.snps==TRUE){
    SNPs.samp = sample(c(1:n.allele.total),nsamp,replace=F)
    allele.dat <- allele.dat[,SNPs.samp]
    #Convert the SNP dataset to a genind object,
    #needed for calculations of population genetics statistics
    geno.inds <- df2genind(allele.dat, ploidy=2, sep="/",NA.char=NA.symbol)
    }

  if(samp.snps==FALSE){
    geno.inds <- df2genind(allele.dat, ploidy=2, sep="/",NA.char=NA.symbol)
    }

  if(qmatrix.prob==TRUE){
    nclusters=nclusters+1
    #Get labels of populations for inds based on STRUCTURE assignment probabilities
    ind.populations<-vector()
    for(i in 1:length(qmatrix[,1])){
      ind.populations[i]<-colnames(qmatrix[,c(2:nclusters)])[which(rmultinom(1,1,prob=qmatrix[i,c(2:nclusters)])>0)]
      }
    ind.populations <- data.frame(ind.populations)
    rownames(ind.populations) <- qmatrix$ind
  }

  if(qmatrix.prob==FALSE){
    ind.populations = ind.populations
    }

  #Include population information in the 'strata' slot
  ind.populations <- ind.populations[order(rownames(allele.dat)),]
  strata(geno.inds) <- data.frame(ind.populations)

  #Calculate the ML inbreeding estimate
  temp = adegenet::inbreeding(geno.inds, res.type='estimate')

  list(temp = temp,ind.populations=ind.populations)
  })

stopCluster(cl)

#Run linear models (ANOVAs) on each replicate and calculate Tukey's HSD groupings
lin.mods = lapply(seq_along(inbreeding.ML),function(x){
  dat = data.frame(temp=inbreeding.ML[[x]]$temp,
                   Population=inbreeding.ML[[x]]$ind.populations)
  mod = aov(temp~as.factor(Population),data=dat)
  tukes = TukeyHSD(mod)
  compares = tukes[[1]]
  list(coefs = compares[,1])
})

#Gather summary statistics from the Tukey's results
lin.coefs = lapply(lin.mods,'[',1) %>%
  data.frame %>%
  dplyr::bind_cols(as.data.frame(t(apply(.,1,quantile,c(0.025, 0.5,0.975))))) %>%
  dplyr::bind_cols(comparison = names(lin.mods[[1]]$coefs)) %>%
  dplyr::select('comparison','2.5%','50%','97.5%')

#Gather summary statistics from raw inbreeding estimates per individual
raw.vals = lapply(inbreeding.ML, '[', 1) %>%
  data.frame %>%
  dplyr::bind_cols(as.data.frame(t(apply(.,1,quantile,c(0.025, 0.5,0.975))))) %>%
  dplyr::bind_cols(individual = rownames(allele.dat)) %>%
  dplyr::select('individual','2.5%','50%','97.5%')

#Return summary matrices in a list
list(lin.coefs=lin.coefs,raw.vals=raw.vals)
}
