# This is a program to get percentage of enterotype stability of samples based on the taxa abundance table
# input = QIIME output table
# N=1000 - number of iterations, can be changed
# clust.number = 3 do not change - can lead to errors
# accuracy is the threshold for a 'stable' enterotype. Can change it! 90% by default: the enterotype must be stable in 90% of subsamplings

getStablEntSubsampl <- function(input_file,N=1000,clust.number=3,accuracy=0.9){
  library(cluster)
  library(clusterSim)
  library(scatterplot3d)
  library(knb16slib)
  require(ade4)
  library(dplyr)
  library(data.table)
  if(clust.number!=3){warning('Enterotype number must be 3! otherwise program can work incorrectly.')}
  # First, make the steps similar to enterotyping Raes. Read table, remove rare bacteria
  #Download the example data and set the working directory
  data=read.table(input_file, header=T, row.names=1, dec=".", sep="\t")
  # select specific location - that can be .IL., .TR. or .RE.
  if (locat!='all'){
    data <- select(data,contains(locat))
  }
  # data=data[-1,]
  
  # # data=noise.removal(data, percent=0.01)
  noise.removal <- function(dataframe, percent=0.01, top=NULL){
    dataframe->Matrix
    bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent
    Matrix_1 <- Matrix[bigones,]
    # print(percent)
    return(Matrix_1)
  }
  
  if(noiserem){data=noise.removal(data, percent=0.01)}
  
  dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
    KLD <- function(x,y) sum(x *log(x/y))
    JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
    matrixColSize <- length(colnames(inMatrix))
    matrixRowSize <- length(rownames(inMatrix))
    colnames <- colnames(inMatrix)
    resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
    
    inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
    
    for(i in 1:matrixColSize) {
      for(j in 1:matrixColSize) { 
        resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                               as.vector(inMatrix[,j]))
      }
    }
    colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
    as.dist(resultsMatrix)->resultsMatrix
    attr(resultsMatrix, "method") <- "dist"
    return(resultsMatrix) 
  }
  
  pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
    require(cluster)
    cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
    return(cluster)
  }
  
  # Now, perform permutations
  print("Permatyping permutations started")

  
  # only for permatyping
  # make a matrix matrix_cocluster. It will (after permutations) contain number of times each pair of samples has same emterotype
  # matrix_cocluster <- matrix(0,dim(data)[2],dim(data)[2])
  
  # Make ent_prob_matrix = matrix with columns of sample name and rows corresp to how many times the sample
  # is Bacteroides, Prevotella and enterotype_3
  ent_prob_matrix <- matrix(0,nrow = 4,ncol = dim(data)[2])
  rownames(ent_prob_matrix) <- c('g__Bacteroides','g__Prevotella','enterotype_3','total')
  colnames(ent_prob_matrix) <- colnames(data)
  
  for (i in 1:N){
    print(paste0('permutation ',i))
    safety_counter = 0
    # select half of the samples
    sample_i <- sample(colnames(data), round(dim(data)[2]/2,0))
    data_subsample <- data[,sample_i]
    data_subsample.dist=dist.JSD(data_subsample)
    data_subsample.cluster=pam.clustering(data_subsample.dist, k=clust.number)
    # now we get a vector sample-enterotype (samples are in names(vector))
    names(data_subsample.cluster) <- sample_i
    
    obs.pca=dudi.pca(data.frame(t(data_subsample)), scannf=F, nf=10)
    # I do not know why 19! But in original program nf=19...
    obs.bet=bca(obs.pca, fac=as.factor(data_subsample.cluster), scannf=F, nf=19)
    
    first_enterotype <- names(data_subsample.cluster[data_subsample.cluster==1])
    second_enterotype <- names(data_subsample.cluster[data_subsample.cluster==2])
    third_enterotype <- names(data_subsample.cluster[data_subsample.cluster==3])
    # make sure that we have 1 cluster driven by Bacteroides and 1 by Prevotella - else skip the iteration!
    tmp_chimera <- paste(GetDrivers(obs.bet[1], n=as.numeric(j))[[2]][[1]][[1]][[as.numeric(j)]],GetDrivers(obs.bet[2], n=as.numeric(j))[[2]][[1]][[1]][[as.numeric(j)]],GetDrivers(obs.bet[3], n=as.numeric(j))[[2]][[1]][[1]][[as.numeric(j)]],sep=';')
    if ((!tmp_chimera %like% 'g__Bacteroides')|(!tmp_chimera %like% 'g__Prevotella')){
      print(paste0('Iteration ',i,': no g__Bacteroides and/or g__Prevotella enterotype'))
      next
    }
    # Fill in the ent_prob_matrix matrix!
    # first enterotype can be either Bacteroides or Prevotella or other. Depending on that, fill in the ent_prob_matrix
    if(GetDrivers(obs.bet[1], n=as.numeric(j))[[2]][[1]][[1]][[as.numeric(1)]] %like% 'g__Bacteroides'){
      ent_prob_matrix["g__Bacteroides",first_enterotype] <- ent_prob_matrix["g__Bacteroides",first_enterotype]+1
    } else if (GetDrivers(obs.bet[1], n=as.numeric(j))[[2]][[1]][[1]][[as.numeric(j)]] %like% 'g__Prevotella'){
      ent_prob_matrix["g__Prevotella",first_enterotype] <- ent_prob_matrix["g__Prevotella",first_enterotype]+1   
    } else{
      ent_prob_matrix["enterotype_3",first_enterotype] <- ent_prob_matrix["enterotype_3",first_enterotype]+1            
      safety_counter = safety_counter+1
    }
    # now, for the second enterotype
    if(GetDrivers(obs.bet[2], n=as.numeric(j))[[2]][[1]][[1]][[as.numeric(1)]] %like% 'g__Bacteroides'){
      ent_prob_matrix["g__Bacteroides",second_enterotype] <- ent_prob_matrix["g__Bacteroides",second_enterotype]+1
    } else if (GetDrivers(obs.bet[2], n=as.numeric(j))[[2]][[1]][[1]][[as.numeric(j)]] %like% 'g__Prevotella'){
      ent_prob_matrix["g__Prevotella",second_enterotype] <- ent_prob_matrix["g__Prevotella",second_enterotype]+1   
    } else{
      ent_prob_matrix["enterotype_3",second_enterotype] <- ent_prob_matrix["enterotype_3",second_enterotype]+1      
      safety_counter = safety_counter+1
      
    }
    # now, for the third enterotype
    if(GetDrivers(obs.bet[3], n=as.numeric(j))[[2]][[1]][[1]][[as.numeric(1)]] %like% 'g__Bacteroides'){
      ent_prob_matrix["g__Bacteroides",third_enterotype] <- ent_prob_matrix["g__Bacteroides",third_enterotype]+1
    } else if (GetDrivers(obs.bet[3], n=as.numeric(j))[[2]][[1]][[1]][[as.numeric(j)]] %like% 'g__Prevotella'){
      ent_prob_matrix["g__Prevotella",third_enterotype] <- ent_prob_matrix["g__Prevotella",third_enterotype]+1   
    } else{
      ent_prob_matrix["enterotype_3",third_enterotype] <- ent_prob_matrix["enterotype_3",third_enterotype]+1
      safety_counter = safety_counter+1
    }
    if (safety_counter>1){
      geterrmessage('Internal error! non-Bacteroides/Prevotella enterotype seen > 1 time at this iteration!')
    }
    ent_prob_matrix["total",sample_i] <- ent_prob_matrix["total",sample_i]+1
    
  }
  
  # Now make a table sample-enterotype: columns = sample, enterotype, stable_enterotype, quality
  # Now we have the info for each sample - in how many cases it got to a certain enterotype. 
  # Need to select samples that have been assigned the same enterotype >accuracy% of cases (90% po umolch)
  sample.enterotype <- matrix(0,dim(ent_prob_matrix)[2],4)
  colnames(sample.enterotype) <- c('Sample','Enterotype','Stable_enterotype','Stability')
  if(!is.numeric(ent_prob_matrix)){ent_prob_matrix <- as.numeric(ent_prob_matrix)}
    for (i in 1:dim(ent_prob_matrix)[2]){
      enterotype <- which.max(ent_prob_matrix[1:3,i])
      if (ent_prob_matrix['total',i]==0){next()}
      quality <- signif(ent_prob_matrix[enterotype,i]/ent_prob_matrix['total',i],digits = 3)
      if(quality >= accuracy){stable_enterotype = names(enterotype)}else{stable_enterotype = NA}
      sample.enterotype[i,] <- c(colnames(ent_prob_matrix)[i],names(enterotype),stable_enterotype,quality)
    }
  return(sample.enterotype)
}




# AIM
# The program to check how the probability of getting v12-v34-v56/v12-v56 stable enterotype (file2) depending on 
# enterotype stability by subsampling in a specific amplicon (say, v1v2, file1)
# enterotype probability obtained by getStablEntSubsampl function
# INPUT
# file1 = a table of enterotype probabilities done by subsampling ('Sample','Enterotype','Stable_enterotype','Stability')
# file2 = sample-stable enterotype across locations
# folder = folder for drawing graph, suggested folder = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype/Enterotype_Probabilities/Graphs_enterotype_stability'
# graph.file = 'v1v2_ent_stability_by_sampl_vs_v12v34v56.jpg' - need to change depending on amplicon. If empty, no graph
# window = size of rolling window (po umolch 200)
# OUTPUT
# 1. Enterotype probabilities sorted in descending order, and for a rolling window (po umolch width 200) we see in output
# the frequency of a stable enterotype. E.g. it can be >0.6 at the beginning and <0.4 at the end
# 2. If graph.file name indicated, the graph is drawn in the file and placed to folder "folder"

subsVsLoc <- function(file1 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype/Enterotype_Probabilities/V1V2_3e_sample-enterotype_probabilities.csv',file2='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype/Stable_enterotypes.csv', folder = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype/Enterotype_Probabilities/Graphs_enterotype_stability',graph.file = '',window=200){
  library(data.table)
  library(zoo)
  setwd(folder)
  ent.with.prob <- read.table(file1)
  ent.stable <- read.table(file2)
  
  # need to delete .V12 etc. from file1
  ent.with.prob[,1] <- gsub('\\.V[1-6][1-6]','',ent.with.prob[,1])
  
  # sort the ent.with.prob according to probability of enterope
  ent.with.prob <- ent.with.prob[order(ent.with.prob$Stability,decreasing = T),]
  # make a vector of 787 values 0 = absent in stable enterotypes, 1 = present
  identity.vector <- ent.with.prob$Sample %in% ent.stable$stable.samples
  stability.window <- rollmean(identity.vector,k=window)

  if(graph.file!=''){
  ylab1 = paste0('Probability to be among to stable enterotypes for the sample, window size ', window)
  jpeg(filename = graph.file,width=1000,height = 600)
  plot(stability.window, xlab = 'Sample enterotype stability index in descend order', ylab = ylab1)
  dev.off()
  }
  print(paste0("Here the enterotype probabilities obtained for a certain amplicon are sorted in descending order and than for these samples frequnency of across-location-stable enterotype calculated - average per rolling window, window size ",window))
  # construct a vector of average enterotype stability for first n elements of identity.vector = first n samples ranked according to enterotype probabilities
  ent.stability.besttoworst <- c(rep(0,times=length(identity.vector)))
  for (i in 1:length(identity.vector)){
    ent.stability.besttoworst[i] <- mean(identity.vector[1:i])
  }
  return(list('stability.window'=stability.window,'identity.vector'=identity.vector,'ent.stability.besttoworst.sum.first.n.samples'=ent.stability.besttoworst))
}
  