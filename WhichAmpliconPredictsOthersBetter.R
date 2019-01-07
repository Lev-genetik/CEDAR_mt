# program to check the predicting value of enterotyping by a specific amplicon compared to others for number of 
# enterotypes (e.g. 2 to 5)
# INPUT:
# table sample-enterotype for v12, v34 and v56 amplicons
# OUTPUT:
# p-value of each amplicon to determine enterotypes for the other amplicons

getPredValue <- function(file.v12,file.v34,file.v56){
  source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/clustering_random_simulation_v3.R')
  library(psych)
  library(stats)
  library(MASS)
  ent.12 <- read.table(file.v12,header = T)
  ent.34 <- read.table(file.v34,header = T)
  ent.56 <- read.table(file.v56,header = T)

  compEnt_v12_v34 <- CompEntPred.pair(ent.12,ent.34)
  compEnt_v12_v56 <- CompEntPred.pair(ent.12,ent.56)
  compEnt_v34_v56 <- CompEntPred.pair(ent.34,ent.56)
  ident.1234 <- as.numeric(compEnt_v12_v34$Maximum_number_samples_clustering_together[1])
  ident.1256 <- as.numeric(compEnt_v12_v56$Maximum_number_samples_clustering_together[1])
  ident.3456 <- as.numeric(compEnt_v34_v56$Maximum_number_samples_clustering_together[1])
  
  # numbers of samples in each enterotype for each amplicon
  clus.numbers.12 <- as.numeric(table(ent.12[,2]))
  clus.numbers.34 <- as.numeric(table(ent.34[,2]))
  clus.numbers.56 <- as.numeric(table(ent.56[,2]))
  
simul.rand.1234 <- simulatn(clus.numbers.12,clus.numbers.34,10000,density = F)
simul.rand.1256 <- simulatn(clus.numbers.12,clus.numbers.56,10000,density = F)
simul.rand.3456 <- simulatn(clus.numbers.34,clus.numbers.56,10000,density = F)
# check that we can approximate the number of matching samples by normal distribution
# if(plots){
#   plot(density(simul.rand.1234$distr))
#   plot(density(simul.rand.1256$distr))
#   plot(density(simul.rand.3456$distr))
# }
# check that we can approximate by normal distr - skew and kurtisis <=1. Else, warning message
if((abs(describe(as.data.frame(simul.rand.1234$distr),trim=0)$skew)>1)|(abs(describe(as.data.frame(simul.rand.1234$distr),trim=0)$kurt)>1)){
  warning(paste('skewness and/or kurtosis are >1 for v12-v34 random simulation: here they are correspodingly',describe(as.data.frame(simul.rand.1234$distr),trim=0)$skew,' and ', describe(as.data.frame(simul.rand.1234$distr),trim=0)$kurt))
}
if((abs(describe(as.data.frame(simul.rand.1256$distr),trim=0)$skew)>1)|(abs(describe(as.data.frame(simul.rand.1256$distr),trim=0)$kurt)>1)){
  warning(paste('skewness and/or kurtosis are >1 for v12-v56 random simulation: here they are correspodingly',describe(as.data.frame(simul.rand.1256$distr),trim=0)$skew,' and ', describe(as.data.frame(simul.rand.1256$distr),trim=0)$kurt))
}
if((abs(describe(as.data.frame(simul.rand.3456$distr),trim=0)$skew)>1)|(abs(describe(as.data.frame(simul.rand.3456$distr),trim=0)$kurt)>1)){
  warning(paste('skewness and/or kurtosis are >1 for v12-v34 random simulation: here they are correspodingly',describe(as.data.frame(simul.rand.3456$distr),trim=0)$skew,' and ', describe(as.data.frame(simul.rand.3456$distr),trim=0)$kurt))
}

# so, we can use normal approximation
# and now, make it! Get mean and SD
mean.1234 <- fitdistr(simul.rand.1234$distr,densfun='normal')$estimate[1]
sd.1234 <- fitdistr(simul.rand.1234$distr,densfun='normal')$estimate[2]
mean.1256 <- fitdistr(simul.rand.1256$distr,densfun='normal')$estimate[1]
sd.1256 <- fitdistr(simul.rand.1256$distr,densfun='normal')$estimate[2]
mean.3456 <- fitdistr(simul.rand.3456$distr,densfun='normal')$estimate[1]
sd.3456 <- fitdistr(simul.rand.3456$distr,densfun='normal')$estimate[2]
# ident.1234 is the number of identical samples between enterotyping classification of v12 and v34
# etc.
pval.1234 <- pnorm(ident.1234,mean.1234,sd.1234,lower.tail = F)
pval.1256<-pnorm(ident.1256,mean.1256,sd.1256,lower.tail = F)
pval.3456 <-pnorm(ident.3456,mean.3456,sd.3456,lower.tail = F)
pred.value.12 <-exp(mean(log(c(pval.1234,pval.1256))))
pred.value.34 <-exp(mean(log(c(pval.1234,pval.3456))))  
pred.value.56 <-exp(mean(log(c(pval.1256,pval.3456))))  
result <- c(pred.value.12,pred.value.34,pred.value.56)
names(result) <- c('v1v2.predict.pval','v3v4.predict.pval','v5v6.predict.pval')
return(result)
}


# for files 1,2 and 3 containing tables sample-enterotype, makes comaprison PAIRWISE
# INPUT:
# files 1,2 and 3 containing tables sample-enterotype
# OUTPUT:
# a list: $ident.number = vector: v12, v34, v56, average
# $ident.percent = vector: v12, v34, v56, average
# $total.sample.number = total number of samples in each of the the files (must be equal!)
getPredNumber <- function (file.v12,file.v34,file.v56){
  source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Enterotyping_compar_functions_v5.R')
  ent.12 <- read.table(file.v12,header = T)
  ent.34 <- read.table(file.v34,header = T)
  ent.56 <- read.table(file.v56,header = T)
  # check that number of samples is the same in the 3 files
  if ((dim(ent.12)[1]!=dim(ent.34)[1])|(dim(ent.12)[1]!=dim(ent.56)[1])){
    stop()
    geterrmessage('number of samples is unequal in the 3 input files! It should be thoigh.')
  }

  # get number of elements identical between amplicons
  compEnt_v12_v34 <- CompEntPred.pair(ent.12,ent.34)
  compEnt_v12_v56 <- CompEntPred.pair(ent.12,ent.56)
  compEnt_v34_v56 <- CompEntPred.pair(ent.34,ent.56)
  # an internal quality check for CompEntPred.pair dunction
  if(compEnt_v12_v34$Maximum_number_samples_clustering_together[3]!=dim(ent.56)[1]){
    stop()
    geterrmessage('internal error in CompEntPred.pair: total number of samples in number of CompEntPred.pair is unequal to number of samples in input')
  }
  ident.number.1234 <- as.numeric(compEnt_v12_v34$Maximum_number_samples_clustering_together[1])
  ident.number.1256 <- as.numeric(compEnt_v12_v56$Maximum_number_samples_clustering_together[1])
  ident.number.3456 <- as.numeric(compEnt_v34_v56$Maximum_number_samples_clustering_together[1])
  # now, calculate average identity number between a particular amplicon enterotyping and enterotyping 
  # by the other 2 amplicons
  ident.number.v12 <- mean(c(ident.number.1234,ident.number.1256))
  ident.number.v34 <- mean(c(ident.number.1234,ident.number.3456))
  ident.number.v56 <- mean(c(ident.number.1256,ident.number.3456))
  # and finaly, the average of average identity number for the 3 amplicons 
  ident.number.average <- round(mean(c(ident.number.v12,ident.number.v34,ident.number.v56)),0)
  # construct the result vector 1 (identity number)
  ident.number <- c(ident.number.v12,ident.number.v34,ident.number.v56,ident.number.average)
  names(ident.number) <- c('v1v2','v3v4','v5v6','average.locat')
  
  # Now, calculate the average percentage of samples in v12 classification that are clasified similarly in v34 and v56
  ident.percent.v12 <- paste0(round(100*ident.number.v12/dim(ent.56)[1],2),'%')
  # the same for v34 and v56
  ident.percent.v34 <- paste0(round(100*ident.number.v34/dim(ent.56)[1],2),'%')
  ident.percent.v56 <- paste0(round(100*ident.number.v56/dim(ent.56)[1],2),'%')
  ident.percent.average <- paste0(round(100*ident.number.average/dim(ent.56)[1],2),'%')

  # construct the result vector 1 (identity percentage)
  ident.percent <- c(ident.percent.v12,ident.percent.v34,ident.percent.v56,ident.percent.average)
  names(ident.percent) <- c('v1v2','v3v4','v5v6','average.locat')
  total.sample.number <- dim(ent.56)[1]
  result <- list('ident.number'=ident.number,'ident.percent'=ident.percent,'total.sample.number'=total.sample.number)
  return(result)
}