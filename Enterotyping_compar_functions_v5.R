# v5.    here, changed CompEntPred.pair function: input - only data frame sample-enterotype and 
# output file name

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

# old program... 
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

#works (09.2018)
#compare enterotypes for a pair of loci taking into account all possible cluster-to-cluster correspondences
# input:
# df1, df2 - data.frame sample-enterotype. Samples should be identical, but may have '.V12' '.V34' or '.V56' tags
# file.name - optional: name of the file for plot of number of samples clustering similar vs cluster correspondence
# Output:
# returns a list of $Maximum_number_samples_clustering_together number of samples that got to the same clusters
# $Optimal_cluster_correspondences correspondence of clusters (first drivers shown)
# Number_of_samples_clustering_together - each row of the data frame has 2 columns: cluster correspondence 
# (e.g. 1234 -> 3214) and number of samples in similar clusters for it
# Plot of $Number_of_samples_clustering_together is made in /media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Compare_enterotypings/
# if 5th argument not empty
# Example: tmp <- CompEntPred.pair(ent.Raes12_[[2]],ent.Raes34_[[2]],as.vector(ent.Raes12_[[1]][,1]),as.vector(ent.Raes34_[[1]][,1]),'CompEntPred_V12_V34.jpg')
CompEntPred.pair <- function(df1,df2,file.name=NULL){
  library('stringr')
  library('combinat')
  library('plyr')
  if (ncol(df1)!=2){
    stop()
    geterrmessage('Error: number of columns in the input data frame 1 should be 2!')
  }
  if (ncol(df2)!=2){
    stop()
    geterrmessage('Error: number of columns in the input data frame 2 should be 2!')
  }
  if(!all.equal(colnames(df1),c('sample','enterotype'))){
    warning('The colnames of the input data frame must be sample and enterotype, now changing')
    colnames(df1)=c('sample','enterotype')
  }
  if(!all.equal(colnames(df2),c('sample','enterotype'))){
    warning('The colnames of the input data frame must be sample and enterotype, now changing')
    colnames(df2)=c('sample','enterotype')
  }
  
  driv1 <- levels(df1[,2])
  driv2 <- levels(df2[,2])
  if (length(driv1)!=length(driv2)){
    stop()
    geterrmessage('Error: number of enterotypes is not the same in the 2 data frames')
  }
  
  # replace enterotypes by numbers
  # print(paste0('enterotypes before adding levels: ',df1$enterotype))
  for (i in 1:length(driv1)){
    dr <- driv1[i]
    # print(paste0('current driver: ',driv1[i]))
    df1$enterotype <- mapvalues(df1$enterotype,driv1[i],as.character(i))
  }
  # options(max.print = 13)
  # print(paste0('df1 enterotypes',df1$enterotype))
  
  for (i in 1:length(driv2)){
    dr <- driv2[i]
    df2$enterotype <- mapvalues(df2$enterotype,driv2[i],as.character(i))
  }
  # print(paste0('df2 enterotypes',df2$enterotype))
  # print(df2$enterotype)
  
  # delete '.V12','.V34' and '.V56' in the sample names
  df1$sample <- str_replace(df1$sample,'\\.V[0-9]{2}','')
  df2$sample <- str_replace(df2$sample,'\\.V[0-9]{2}','')
  
  # delete samples that are different between df1 and df2
  df1 <- df1[df1$sample %in% intersect(df1$sample,df2$sample),]
  df2 <- df2[df2$sample %in% intersect(df1$sample,df2$sample),]
  
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
    setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Compare_enterotypings')
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



# Aim: get correlations BETWEEN DRIVERS for different enterotype pairs (all in all 9 enterotypes as Ent1-V12 to Ent3-V56)
# for 3 enterotype clusters only
# result: matrix of 9x9
# check=FALSE. If check=TRUE, we check correlations between enterotypes in the same location using all 
# 92-99 bacteria. To compare correlation with the one obtained using the 68 geni common between
# all amplicons. I think <0.1 difference is okay.
# all in all, we have 92 drivers - in case of 3 enterotypes
# remove.driver=NULL by deafult, we can make it e.g. remove.driver=c('g__Bacteroides','g__Prevotella') to remove these 2 
# drivers from the analysis to check Michel's idea.
getEntCorr <- function(file1='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt', file2='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt', file3='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt', correl.method = 'spearman', locat='all', rnd = 3, check=F,remove.driver=NULL) {
 library(stringr)
  source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Enterotyping_Raes_noise_removed_v9.R')
  
  if (!is.numeric(rnd)){
    stop()
    geterrmessage('rnd argument should be integer! It is the number of digits for values kept in the result matrix')
  }
  
  # make enterotyping
  ent.Raes12.3e <- Enterotyping_Raes(file1, '/media/lev-genetik/980E73270E72FD96/Liege/Results/Clustering/V1V2', clust.number=3, noiserem=T,  locat = locat, driv.number = 1000, plots=F)
  ent.Raes34.3e <- Enterotyping_Raes(file2, '/media/lev-genetik/980E73270E72FD96/Liege/Results/Clustering/V3V4', clust.number=3, noiserem=T,  locat = locat, driv.number = 1000, plots=F)
  ent.Raes56.3e <- Enterotyping_Raes(file3, '/media/lev-genetik/980E73270E72FD96/Liege/Results/Clustering/V5V6', clust.number=3, noiserem=T,  locat = locat, driv.number = 1000, plots=F)
  
  # get bacteria that are present in all enterotypings (V1V2, V3V4 and V5V6)
  names12 <- colnames(ent.Raes12.3e$obs_bet_result$tab)
  names34 <- colnames(ent.Raes34.3e$obs_bet_result$tab)
  names56 <- colnames(ent.Raes56.3e$obs_bet_result$tab)
  names_intersect <- intersect(names56,intersect(names12,names34))
  # if remove.driver!= NULL, we need to remove the drivers
  if (!is.null(remove.driver)){
    for (i in 1:length(remove.driver)){
      names_intersect <- names_intersect[!names_intersect %in% grep(remove.driver[i],names_intersect,value=T)]
    }
  }
  
  # get vectors of average bacterial abundance per eneterotype. Bacteria in the order done by bca
  e12_1 <- as.numeric(as.vector(ent.Raes12.3e$obs_bet_result$tab[1,]))
  e12_2 <- as.numeric(as.vector(ent.Raes12.3e$obs_bet_result$tab[2,]))
  e12_3 <- as.numeric(as.vector(ent.Raes12.3e$obs_bet_result$tab[3,]))
  e34_1 <- as.numeric(as.vector(ent.Raes34.3e$obs_bet_result$tab[1,]))
  e34_2 <- as.numeric(as.vector(ent.Raes34.3e$obs_bet_result$tab[2,]))
  e34_3 <- as.numeric(as.vector(ent.Raes34.3e$obs_bet_result$tab[3,]))
  e56_1 <- as.numeric(as.vector(ent.Raes56.3e$obs_bet_result$tab[1,]))
  e56_2 <- as.numeric(as.vector(ent.Raes56.3e$obs_bet_result$tab[2,]))
  e56_3 <- as.numeric(as.vector(ent.Raes56.3e$obs_bet_result$tab[3,]))
  names(e12_1) <- names12
  names(e12_2) <- names12
  names(e12_3) <- names12
  names(e34_1) <- names34
  names(e34_2) <- names34
  names(e34_3) <- names34
  names(e56_1) <- names56
  names(e56_2) <- names56
  names(e56_3) <- names56
  
  # out of the abundance vectors, only data on the 68 [in case all samples enterotyped together and no geni removed] 
  # bacteria geni common within amplicons should remain
  e12_1i <- e12_1[names_intersect]
  e12_2i <- e12_2[names_intersect]
  e12_3i <- e12_3[names_intersect]
  e34_1i <- e34_1[names_intersect]
  e34_2i <- e34_2[names_intersect]
  e34_3i <- e34_3[names_intersect]
  e56_1i <- e56_1[names_intersect]
  e56_2i <- e56_2[names_intersect]
  e56_3i <- e56_3[names_intersect]
  
  
  # get the correlation matrix
  enterotypes.df <- cbind(e12_1i,e12_2i,e12_3i,e34_1i,e34_2i,e34_3i,e56_1i,e56_2i,e56_3i)
  ent.correlation <- matrix(0,ncol(enterotypes.df),ncol(enterotypes.df))
  for (i in 1:ncol(enterotypes.df)){
    for (j in 1:ncol(enterotypes.df)){
      ent.correlation[i,j] <- round(cor(enterotypes.df[,i],enterotypes.df[,j],method=correl.method),rnd)
    }
  }
  
  # Need to name all the 9 rows. 1st driver: take only the part of the string after '_'
  rownames(ent.correlation) <- c(rep('ND',dim(ent.correlation)[1]))
  rownames(ent.correlation)[1] <- paste0('V1V2_ent1: ', word(ent.Raes12.3e$drivers_n[1,1],start=-1,sep='_'))
  rownames(ent.correlation)[2] <- paste0('V1V2_ent2: ', word(ent.Raes12.3e$drivers_n[2,1],start=-1,sep='_'))
  # Here a manual correction due to no specific genus for V1V2_ent3 1st driver
  rownames(ent.correlation)[3] <- paste0('V1V2_ent3: ', word(ent.Raes12.3e$drivers_n[3,1],start=-15,sep='_')[10])
  
  rownames(ent.correlation)[4] <- paste0('V3V4_ent1: ', word(ent.Raes34.3e$drivers_n[1,1],start=-1,sep='_'))
  rownames(ent.correlation)[5] <- paste0('V3V4_ent2: ', word(ent.Raes34.3e$drivers_n[2,1],start=-1,sep='_'))
  rownames(ent.correlation)[6] <- paste0('V3V4_ent3: ', word(ent.Raes34.3e$drivers_n[3,1],start=-1,sep='_'))
  
  rownames(ent.correlation)[7] <- paste0('V5V6_ent1: ', word(ent.Raes56.3e$drivers_n[1,1],start=-1,sep='_'))
  rownames(ent.correlation)[8] <- paste0('V5V6_ent2: ', word(ent.Raes56.3e$drivers_n[2,1],start=-1,sep='_'))
  rownames(ent.correlation)[9] <- paste0('V5V6_ent3: ', word(ent.Raes56.3e$drivers_n[3,1],start=-1,sep='_'))
  
  colnames(ent.correlation) <- abbreviate(rownames(ent.correlation),minlength = 8)
  
  
  # If check=TRUE only - get precise correlations of enterotypes within each amplicon
  if (check){
    # first define a matrix of correlations: we already have row and column names!
    ent.correlation.precise <- ent.correlation
    ent.correlation.precise[,] <- NA
    
    
    enterotypes12.df <- cbind(e12_1,e12_2,e12_3)
    ent.correlation12 <- round(cor(enterotypes12.df,method=correl.method),rnd)
    enterotypes34.df <- cbind(e34_1,e34_2,e34_3)
    ent.correlation34 <- round(cor(enterotypes34.df,method=correl.method),rnd)
    enterotypes56.df <- cbind(e56_1,e56_2,e56_3)
    ent.correlation56 <- round(cor(enterotypes56.df,method=correl.method),rnd)
    ent.correlation.precise[1:3,1:3] <- ent.correlation12
    ent.correlation.precise[4:6,4:6] <- ent.correlation34
    ent.correlation.precise[7:9,7:9] <- ent.correlation56
    
    
    # just to double-check
    ent.correlation12a <- matrix(0,ncol(enterotypes12.df),ncol(enterotypes12.df))
    for (i in 1:ncol(enterotypes12.df)){
      for (j in 1:ncol(enterotypes12.df)){
        ent.correlation12a[i,j] <- round(cor(enterotypes12.df[,i],enterotypes12.df[,j],method=correl.method),rnd)
      }
    }
   
    result <- list(ent.correlation,ent.correlation.precise,ent.correlation12a)
    names(result) <- c('correlation matrix for all enterotypes','precise correlation matrices enterotypes inside amplicon','double-check V12 amplicon')
  } else{
    result <- ent.correlation
  }
  return(result)
}
