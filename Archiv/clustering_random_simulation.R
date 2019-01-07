#check if the enterotyping with v1v2 v3v4 and v5v6 is not noise: do samples tend to cluster similarly?
#ent.Raes12_ etc are results of /home/lev-genetik/Desktop/Projects/liege/src/Lev/Enterotyping_Raes_noise_removed_v5.R
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Enterotyping_compar_functions_v2.R')
freq_12 <- getEntFreq(ent.Raes12_[[2]],as.vector(ent.Raes12_[[1]][,1]))
freq_34 <- getEntFreq(ent.Raes34_[[2]],as.vector(ent.Raes34_[[1]][,1]))
freq_56 <- getEntFreq(ent.Raes56_[[2]],as.vector(ent.Raes56_[[1]][,1]))

# Simulates n=sum(pr1)=sum(pr2) samples randomly clustered to k=length(pr1)=length(pr2) clusters (numbers of elements in each cluster are in pr1 and pr2). 
# result[[1]]: the distribution of number of elements getting to similar cluster (taking into account the 'best' of all cluster number correspondences)
# result[[2]], result[[3]] - the number of identical elements such that probability to get higher number of identical elements is <=5% and <=1% respectively
# Get maximal number of elements in the same cluster in 2 random series (taking into account all possible cluster correspondences)
simulatn <- function(pr1,pr2,permut.n){
  if(sum(pr1)!=sum(pr2)){
   stop()
    geterrmessage('Error: arrays for simulation of number of identical elements have different sizes!')
  }else {n<-sum(pr1)}
  if(length(pr1)!=length(pr2)){
    stop()
    geterrmessage('Error: arrays for simulation of number of identical elements have different number of clusters!')
  }else{k = length(pr1)}

      #construct a vector of all possible combinations of cluster correspondences
  correspondence <- t(array(unlist(permn(k)), dim = c(k, factorial(k))))
  
  #construct 2 vectors of clustering (according to each custer number) and shuffle them randomly
  v1 <- integer(n)
  v2 <- integer(n)
  st <- 1
  fin <- 0
  for (i in 1:k){
    fin <-fin + pr1[i]
    v1[st:fin] <- i
    st <- fin + 1
  }
  st <- 1
  fin <- 0
  for (i in 1:k){
    fin <-fin + pr2[i]
    v2[st:fin] <- i
    st <- fin + 1
  }  
  #get the result
  result <- integer(n)
  for (i in 1:permut.n){
    v1s <- sample(v1)
    v2s <- sample(v2)
    # here need to compare the vectors using specific comparison rule: i corresponds to correspondence[cycle_number,i]
    # and get maximum number of identical vectors over all possible rules in correspondce
    result.current.perm <- 0
    for (j in 1:nrow(correspondence)){
      count <- 0
      for (l in 1:length(v1s)){
        if(v1s[l]==correspondence[j,v2s[l]]){
          count = count + 1
        }
      }
      result.current.perm = max(result.current.perm,count)
    }
    result[i] <- result.current.perm
  }
  distr <- table(result)
  summ <- sum(distr)
  count <- 0
  # numbers of identical elements corresponding to 5 percentile and 1 percentile
  five_pers <- summ/20
  one_pers <- summ/100
  # 5 percentile and 1 percentile
  five_percentile <- 0
  one_percentile <- 0
  # choose the percentiles = the number of clusters such that probability to get higher number of clusters is <=5% and <=1% respectively
  for(i in length(distr):1){
    count <- count + distr[i]
    if((count > five_pers)&(five_percentile==0)){
      five_percentile <- labels(distr[i])
    }
    if((count > one_pers)&(one_percentile==0)){
      one_percentile <- labels(distr[i])
    }
  }
  result1 <- list(distr,five_percentile,one_percentile)
}

simulation_V12_V34 <- simulatn(freq_12,freq_34,permut.n = 10000)
simulation_V12_V56 <- simulatn(freq_12,freq_56,permut.n = 10000)
simulation_V34_V56 <- simulatn(freq_34,freq_56,permut.n = 10000)
