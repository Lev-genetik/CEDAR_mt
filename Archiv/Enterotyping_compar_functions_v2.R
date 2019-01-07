
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
CompEntPred.pair <- function(df1,df2,driv1,driv2,file.name){
  library('stringr')
  library('combinat')
  if (length(driv1)!=length(driv2)){
    stop()
    geterrmessage('Error: number of enterotypes is not the same in the 2 data frames')
  }
  
  # replace enterotypes by numbers
  for (i in 1:length(driv1)){
    dr <- driv1[i]
    levels(df1$enterotype) <- c(1:length(driv1),levels(df1$enterotype))
    df1$enterotype[df1$enterotype==dr] <- i
  }
  for (i in 1:length(driv2)){
    dr <- driv2[i]
    levels(df2$enterotype) <- c(1:length(driv2),levels(df2$enterotype))
    df2$enterotype[df2$enterotype==dr] <- i
  }  
  
  # delete '.V12','.V34' and '.V56' in the sample names
  df1$sample <- str_replace(df1$sample,'\\.V[0-9]{2}','')
  df2$sample <- str_replace(df2$sample,'\\.V[0-9]{2}','')
  
  #construct a vector of all possible combinations of cluster correspondences
  correspondence <- t(array(unlist(permn(length(driv1))), dim = c(length(driv1), factorial(length(driv1)))))
  
  # calculate maximal number of samples in the same cluster whatever cluster correspondence
  identity.df <- 0
  numb.ident.clustered.samples <- character(length = nrow(correspondence)) 
  cluster.order <- character(length = nrow(correspondence)) 
  for (j in 1:nrow(correspondence)){
    count.conc <- 0
    count.disc <- 0
    for (l in 1:nrow(df1)){
      sample.curr <- df1[l,1]
      ent.df2.curr <- as.numeric(df2[df2[,1]==sample.curr,2])
      if (!is.numeric(ent.df2.curr)){
        stop('enterotype not numeric')
        geterrmessage()
      }
      if(df1[l,2]==correspondence[j,ent.df2.curr]){
        count.conc = count.conc + 1
      } else{
        count.disc = count.disc + 1
      }
    }
    count = count.conc+count.disc
    if(count!=nrow(df1)){
      print(paste0('correspondence combination = ',correspondence[j,]))
      stop()
      geterrmessage('Error: counting matches in the 2 data frames')
    }
    numb.ident.clustered.samples[j] <- as.character(count.conc)
    cluster.order[j] <- as.character(paste(correspondence[j,],collapse=""))
    if (count.conc > identity.df){
      identity.df <-count.conc
      enterotype.correspondence <- rbind(driv1[correspondence[1,]],driv2[correspondence[j,]])
    }else if (count.conc==identity.df){
      # print(paste0('Similar number of identical elements for 2 different cluster correspondence variants. 1,2,3,4 to ',correspondence[j,],', number of similarly clustered elements',identity.df))
    }
  }
  identity.all.clust.corresp <- cbind(cluster.order,numb.ident.clustered.samples)
  result <- list(identity.df,enterotype.correspondence,identity.all.clust.corresp)
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