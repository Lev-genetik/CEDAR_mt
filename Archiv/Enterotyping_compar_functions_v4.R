
# get number of samples for each enterotype (enterotype calculated according to its first driver)
# input: a data frame sample-enterotype e.g. ent.Raes12_[[2]] and a vector with frist drivers, e.g. as.vector(ent.Raes12_[[1]][,1])
#the function works fine! (checked 14.08.18)
getEntFreq <- function(df1,drivers){
  result <- integer(length(drivers))
  for (i in 1:nrow(df1)){
    for (j in 1:length(drivers)){
      if (df1[i,2]==drivers[j]){
        result[j] = result[j]+1
        # print('sovpalo!')
        break
      } else if (j==length(drivers)){
        print('Error: no match for the enterotype in one of the rows')
      }
    }
  }
  # summa <- sum(result)
  # for(i in 1:length(result)){
  #   result[i] <- round(result[i]/summa,4)
  # }
  names(result) <- drivers
  return(result)
}

concordance <- function(df1,df2){
  library(pracma)
  
  n_row_df1 <- as.integer(nrow(df1))
  n_row_df2 <- as.integer(nrow(df2))
  result <- integer(n_row_df1)
  
  for (n in 1:n_row_df1){
    current_sample <- df1[n,1]
    ent1_tmp <- df1[n,2]
    # print(current_sample)
    
    if(current_sample %in% df2$sample){
      ent2_tmp <- df2[df2$sample==current_sample,2]
      # print(ent2_tmp)
      if (strcmp(as.character(ent1_tmp),as.character(ent2_tmp))){
        result[n] <- 'concordant'
        print(paste0('concordant, n= ',n, result[n]))
      }else{result[n] <- 'discordant'
      print(paste0('discordant, n= ',n, result[n]))
      }
    } else{print('some data not available for this sample')}
  }
  summary <- table(result)
  con <- suppressWarnings(as.numeric(summary('concordant')))
  dis <- suppressWarnings(as.numeric(summary('discordant')))
  concord <- con/(con + dis)
  return(result)
}

#works (16.08.18)
#compare enterotypes for a pair of loci taking into account all possible cluster-to-cluster correspondences
# input: df1, df2 - data.frame sample-enterotype. Samples should be identical, but may have '.V12' '.V34' or '.V56' tags
# driv1, driv2 - vectors containing first driver of each enterotype in df1 and df2 correspondingly
# file.name - optional: name of the file for plot of number of samples clustering similar vs cluster correspondence
# returns a list of $Maximum_number_samples_clustering_together number of samples that got to the same clusters
# $Optimal_cluster_correspondences correspondence of clusters (first drivers shown)
# Number_of_samples_clustering_together - each row of the data frame has 2 columns: cluster correspondence 
# (e.g. 1234 -> 3214) and number of samples in similar clusters for it
# Plot of $Number_of_samples_clustering_together is made in /media/lev-genetik/980E73270E72FD96/Liege/Results/Clustering/Compare/
# if 5th argument not empty
# Example: tmp <- CompEntPred.pair(ent.Raes12_[[2]],ent.Raes34_[[2]],as.vector(ent.Raes12_[[1]][,1]),as.vector(ent.Raes34_[[1]][,1]),'CompEntPred_V12_V34.jpg')
CompEntPred.pair <- function(df1,df2,driv1.df,driv2.df,file.name){
  library('stringr')
  library('combinat')
  library('plyr')
  if (length(driv1)!=length(driv2)){
    stop()
    geterrmessage('Error: number of enterotypes is not the same in the 2 data frames')
  }
  driv1 <- driv1.df[,1]
 

  # need to use anyDuplicated to check if for any enterotypes we have the same drivers. In that case, add second driver to driv1!
  if(anyDuplicated(driv1)){
    i=2
    while((anyDuplicated(driv1))&(i<=ncol(driv1.df))){
     for (t in 1:nrow(driv1.df)){
      driv1[t] <- paste(driv1[t],driv1.df[t,i],sep='+++')
     }    
    }
    
  }
  # replace enterotypes by numbers
  for (i in 1:length(driv1)){
      # dr <- driv1[i]
      # print(paste0('current driver: ',driv1[i]))
      df1$enterotype <- mapvalues(df1$enterotype,driv1[i],as.character(i))
  }
  
  driv2 <- driv2.df[,1] 
  # need to use anyDuplicated to check if for any enterotypes we have the same drivers. In that case, add second driver to driv2!
  if(anyDuplicated(driv2)){
    i=2
    while((anyDuplicated(driv2))&(i<=ncol(driv2.df))){
      for (t in 1:nrow(driv2.df)){
       driv2[t] <- paste(driv2[t],driv2.df[t,i],sep='+++')
      }    
    }
  }
  for (i in 1:length(driv2)){
    # dr <- driv2[i]
    df2$enterotype <- mapvalues(df2$enterotype,driv2[i],as.character(i))
  }

  
  # delete '.V12','.V34' and '.V56' in the sample names
  df1$sample <- str_replace(df1$sample,'\\.V[0-9]{2}','')
  df2$sample <- str_replace(df2$sample,'\\.V[0-9]{2}','')
  
  #construct a vector of all possible combinations of cluster correspondences
  correspondence <- t(array(unlist(permn(length(driv1))), dim = c(length(driv1), factorial(length(driv1)))))
  
  # calculate maximal number of samples in the same cluster whatever cluster correspondence
  identity.number <- 0
  numb.ident.clustered.samples <- character(length = nrow(correspondence)) 
  cluster.order <- character(length = nrow(correspondence)) 
  for (j in 1:nrow(correspondence)){
    count.conc <- 0
    count.disc <- 0
    count.na <- 0
    for (l in 1:nrow(df1)){
      sample.curr <- df1[l,1]
# this 'if' checks for df1 samples absent in df2. PROGRAM IS MUCH FASTER IF WE DISABLE IT
      # if((is.na(table(unlist(df2))[sample.curr]))&(j==1)){
      #   print(paste0('Warning!',table(unlist(df2))[sample.curr]))
      #   count.na = count.na + 1
      #   next()
      # }else if (is.na(table(unlist(df2))[sample.curr])){
      #   count.na = count.na + 1
      #   next()
      # }
      # get enterotype of sample.curr in df2
      ent.df2.curr <- as.numeric(levels(df2[df2[,1]==sample.curr,2])[df2[df2[,1]==sample.curr,2]])
      if (!is.numeric(ent.df2.curr)){
        stop('enterotype not numeric')
        geterrmessage()
      }
      if(as.numeric(levels(df1[l,2])[df1[l,2]])==correspondence[j,ent.df2.curr]){
        count.conc = count.conc + 1
      } else{
        count.disc = count.disc + 1
      }
    }
    count = count.conc+count.disc+count.na
    if(count!=nrow(df1)){
      print(paste0('correspondence combination = ',correspondence[j,]))
      stop('Error while counting matches in the 2 data frames')
      geterrmessage()
    }
    numb.ident.clustered.samples[j] <- as.character(count.conc)
    cluster.order[j] <- as.character(paste(correspondence[j,],collapse=""))
    if (count.conc > identity.number){
      identity.number <-count.conc
      enterotype.correspondence <- rbind(driv2,driv1[correspondence[j,]])
      print(paste0('Optimal correspondence of clusters: ',j))
    }else if (count.conc==identity.number){
      # print(paste0('Similar number of identical elements for 2 different cluster correspondence variants. 1,2,3,4 to ',correspondence[j,],', number of similarly clustered elements',identity.number))
    }
  }
  total.samples <- nrow(df1)
  identity.percent <- paste(round(100*identity.number/total.samples, 2), "%", sep="")
  print(identity.percent)
  identity.result <- cbind(identity.number,identity.percent,total.samples)
  print(identity.result)
  identity.all.clust.corresp <- cbind(cluster.order,numb.ident.clustered.samples)
  result <- list(identity.result,enterotype.correspondence,identity.all.clust.corresp)
  names(result) <- c('Maximum_number_samples_clustering_together', 'Optimal_cluster_correspondences', 'Number_of_samples_clustering_together')

  if(!is.null(file.name)){
    setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Clustering/Compare')
    jpeg(filename = file.name)
    plot(c(1:nrow(correspondence)),as.numeric(identity.all.clust.corresp[,2]),type='l',xlab='Cluster correspondence',ylab='Number of samples in corresponding clusters')
    dev.off()
  }
  return(result)
  # else if (!is.null(file.name)){
  #   stop('File name (last argument) must contain jpeg extension!')
  #   geterrmessage()
  # }
}




#NOT WORKING YET
#compare enterotype predictions.
#n number of data frames with 2 columns: sample and enterotype
#make a distance matrix between the samples (distance=number of identical enterotypes)
CompEntPred.matrix <- function(pred.list){
  library(pracma)
    #get rid of .V12, .V34 and .V56
  if (is.list(pred.list)==F){print ('Error: argument should be a list')
    stop()
    geterrmessage('Error: argument should be a list')
  }
  for (i in 1:length(pred.list)){
    pred.list[[i]]$sample <- gsub('.V12','',pred.list[[i]]$sample)
    pred.list[[i]]$sample <- gsub('.V34','',pred.list[[i]]$sample)
    pred.list[[i]]$sample <- gsub('.V56','',pred.list[[i]]$sample)
  }
  
  # fill in the resul matrixmatrix
  result <- matrix(0,length(pred.list),length(pred.list))
  for (i in 1:length(pred.list)){
    for (j in 1:length(pred.list)){
      #now get concordance percent for 2 data frames
      result(i,j) <- round(concordance(pred.list[[i]],pred.list[[j]]),3)
    }
  }  
  return(result)  
  
}