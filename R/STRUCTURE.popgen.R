#' Calculate population genetics parameters
#'
#' This function randomly assigns individuals to clusters using STRUCTURE
#' assignment probabilities provided in the qmatrix. The qmatrix must have
#' column 'ind', which stores individual ids.
#'
#' @importFrom magrittr %>%

#' @param qmatrix Dataframe. The STRUCTURE qmatrix, with individuals in column 'ind'
#' @param allele.dat Dataframe. Raw allele data (sep = '/')
#' @param nclusters Integer. The number of clusters identified by STRUCTURE
#' @param samp.snps Logical. If TRUE, a specified number SNPs will be randomly sampled
#'in each iteration to account for uncertainty and speed computations. Default is FALSE
#' @param nsamp (optional) Integer of SNPs to sample during each iteration
#' @param stat The population genetics statistic to return. Must be one of the following: 'ind.dist',
#' 'summary', 'allele.freqs',
#' 'heterozygosity', 'richness','fis','fst'. \cr\cr
#' Note that statistics 'ind.dist','summary' and 'heterozygosity' do not require a qmatrix, while
#' statistics 'allele.freqs', 'richness','fis' and 'fst' must have a supplied qmatrix.
#' @param NA.symbol The code used to represent missing allele data
#' @return Summaries based on the chosen population genetic statistic
#' @export

STRUCTURE.popgen=function(qmatrix=NULL,allele.dat,nclusters,
                          stat,NA.symbol,
                          samp.snps,nsamp){

  #Collect labels of the individuals
  ind <- rownames(allele.dat)

  #Replace any full stops "." in allele names. Names can't have "." in genind objects
  colnames(allele.dat) <- gsub("\\.", "_", colnames(allele.dat))

  #Replace any missing values with the NA.symbol
  allele.dat[is.na(allele.dat)] <- NA.symbol

  #Basic summary stats that don't require population assignments
  if(stat=="ind.dist"|stat=='summary'|stat=='heterozygosity'){
    Mydata1 <- df2genind(allele.dat, ploidy = as.integer(2),NA.char=NA.symbol,
                         ind.names = ind, sep = "/",
                         type="codom")

    if(stat=="summary"){
      Mydata2 <- genind2hierfstat(Mydata1,pop = rep(1, nrow(allele.dat)))
      summ.stats <- basic.stats(Mydata2, diploid = TRUE, digits = 2)
      output=list(locus.stats = summ.stats$perloc[,c(1:3,9)],
                  average.stats = summ.stats$overall[c(1:3,9)])
    }

    if(stat=="heterozygosity"){
      #Calculate obs. and exp. heterozygosity per locus
      div <- adegenet::summary(Mydata1)
      output=list(N.alleles = nAll(Mydata1),
                  HardyWeinberg.test = pegas::hw.test(Mydata1, B = 1000),
                  H.obs = div$Hobs, H.exp = div$Hexp,
                  Ttest = t.test(div$Hexp,div$Hobs,pair=T,
                                 var.equal=TRUE,alter="greater"))
    }

    if(stat=='ind.dist'){
      #Calculate individual genetic distance following Kosman and Leonard (2005)
      pca.n <- gd.kosman(Mydata1)

      #Calculate genetic distance matrix using normalised euclidian coordinates
      output <- dist(pca.n$geneticdist)
    }
    }

#Stats that do require population assignments must have a qmatrix supplied
  if(stat=='richness'|stat=='fis'|stat=='fst'|stat=="allele.freqs"){

    #Specify default values for optional arguments
    if(missing(samp.snps)){
      samp.snps = FALSE
    }

    if(missing(nsamp)){
      nsamp = ncol(allele.dat)/2
    }

    if(is.null(qmatrix)){
      stop('STRUCTURE qmatrix missing')
      }

    if(missing(nclusters)){
      nclusters = ncol(qmatrix)-1
      }

    #Use same order of individuals for STRUCTURE and allele data
    allele.dat = allele.dat[order(qmatrix$ind),]

    #Get population labels for individuals based on STRUCTURE assignment probabilities
    nclusters=nclusters+1
    population<-vector()
    for(i in 1:length(qmatrix[,1])){
      population[i]<-colnames(qmatrix[,c(2:nclusters)])[which(rmultinom(1,1,prob=qmatrix[i,c(2:nclusters)])>0)]
    }

  if(samp.snps==TRUE){
    n.allele.total = ncol(allele.dat)
    SNPs.samp = sample(c(1:n.allele.total),nsamp,replace=F)
    allele.dat <- allele.dat[,SNPs.samp]
    #Convert the SNP dataset to a genind object,
    #needed for calculations of population genetics statistics
    Mydata1 <- df2genind(allele.dat, ploidy = as.integer(2),NA.char=NA.symbol,
                         ind.names = ind, pop = population, sep = "/",
                         type="codom")
  }

  if(samp.snps==FALSE){
    Mydata1 <- df2genind(allele.dat, ploidy = as.integer(2),NA.char=NA.symbol,
                         ind.names = ind, pop = population, sep = "/",
                         type="codom")
  }

  if(stat=="richness"){
    allele.richness = allel.rich(Mydata1)
    output = data.frame(allele.richness$all.richness)
    output$locus = rownames(output)
  }

  if(stat=="fis"){
    Mydata2 <- genind2hierfstat(Mydata1)
    summ.stats <- basic.stats(Mydata2, diploid = TRUE, digits = 2)

    #Run an ANOVA on cluster-level Fis statistiscs
    Fis.stats = tidyr::gather(data.frame(summ.stats$Fis),cluster,Fis) %>%
      na.omit()
    mod = aov(Fis~as.factor(cluster),data=Fis.stats)
    tukes = TukeyHSD(mod)
    compares = tukes[[1]]
    output <- list(summ.stats=summ.stats,Fis.coefs = compares[,1],Raw.mod = compares)
  }

  if(stat=="allele.freqs"){
    Mydata2 <- genind2hierfstat(Mydata1)
    summ.stats <- nb.alleles(Mydata2, diploid = TRUE)
    all.freqs <- data.frame(summ.stats)
    colnames(all.freqs) <- levels(Mydata2$pop)
    Ho.stats <- data.frame(basic.stats(Mydata2, diploid = TRUE, digits = 2)$Ho)
    Hs.stats <- data.frame(basic.stats(Mydata2, diploid = TRUE, digits = 2)$Hs)
    all.freqs$locus <- rownames(all.freqs)
    Ho.stats$locus <- rownames(all.freqs)
    Hs.stats$locus <- rownames(all.freqs)
    output <- list(obs.richness = all.freqs, Ho.stats = Ho.stats,
                   Hs.stats = Hs.stats)
  }

  if(stat=="fst"){
    #Compare Fst between populations
    output=data.frame(pairwise.fst(Mydata1,res.type = "matrix"))
    output$category.comparison = rownames(output)
  }
  }

  return(output)
}
