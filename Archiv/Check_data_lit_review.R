# check our data by lit review and Michel's tasks

# Get species from QIIME file with specific location and % total abundance > 0.01% of different threshold
getSpec <- function(input_file, locat = 'all', percent_threshold = 0.01){
  library(dplyr)
  data=read.table(input_file, header=T, row.names=1, dec=".", sep="\t")
  if (locat!='all'){
    data <- dplyr::select(data,contains(locat))
  }
  noise.removal <- function(dataframe, top=NULL){
    
    dataframe->Matrix
    bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent_threshold
    Matrix_1 <- Matrix[bigones,]
    
    # dataframe->Matrix
    # rowsums_dolya <- rowSums(Matrix)*100/(sum(rowSums(Matrix)))
    # bigones <- rowsums_dolya[rowsums_dolya > percent_threshold]
    # Matrix_1 <- Matrix[bigones,]
    return(Matrix_1)
  }
  
  data=noise.removal(data)
  result <- rownames(data)
  return(result)
}

# Get Shannon entropy of the samples
# INPUT:
# input_file = a file with QIIME outut, locat = 'all','.IL','.TR' or '.RE', percent_threshold = noise filtration 
# percentage, def. 0.01%
# OUTPUT:
# a vector of Shannon entropy per sample
entrShannon <- function(input_file, locat = 'all', percent_threshold = 0.01){
  library(dplyr)
  data=read.table(input_file, header=T, row.names=1, dec=".", sep="\t")
  if (locat!='all'){
    data <- dplyr::select(data,contains(locat))
  }
  noise.removal <- function(dataframe, top=NULL){
    dataframe->Matrix
    bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent_threshold
    Matrix_1 <- Matrix[bigones,]
    return(Matrix_1)
  }
  data=noise.removal(data)
  data_col_sum <- colSums(data)
  fraction.table <- data
  for (i in 1:dim(data)[2]){
    fraction.table[,i] <- data[,i]/data_col_sum[i]
  }
  entr.shannon <- c(rep(0,dim(data)[2]))
  for (i in 1:dim(fraction.table)[2]){
    for (j in 1:dim(fraction.table)[1]){
      if(fraction.table[j,i]!=0){
      entr.shannon[i] <- entr.shannon[i]-fraction.table[j,i]*logb(fraction.table[j,i])
      }
    }
  }
  names(entr.shannon) <- colnames(data)
  return(entr.shannon)
}


# Compare taxa in three locations (or amplicons)
library(stringr)
library(ggplot2)
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Functions_general_non_enterotyping.R')
compTaxa <- function(file1,loc1,file2,loc2,file3,loc3,percent_threshold = 0.01,remove.outliers = 0.01){
  df1 <- readQiimeSmart(input_file = file1, locat = loc1, percent_threshold = percent_threshold, remove.outliers = remove.outliers)
  df2 <- readQiimeSmart(input_file = file2, locat = loc2, percent_threshold = percent_threshold, remove.outliers = remove.outliers)
  df3 <- readQiimeSmart(input_file = file3, locat = loc3, percent_threshold = percent_threshold, remove.outliers = remove.outliers)
  # appendix
  locat_trim <- function(locat){
    if(locat=='.IL'){result = 'IL'}
    else if(locat=='.TR'){result = 'TR'}
    else if(locat=='.RE'){result = 'RE'}
    else if(locat=='all'){result = 'all'}
    return(result)
  }
  
  # in each column, need to obtain fraction of each bacteria
  # df1
  df1_col_sum <- colSums(df1)
  df1_norm <- df1
  for (i in 1:dim(df1)[2]){
    df1_norm[,i] <- df1[,i]/df1_col_sum[i]
  }
  # df2
  df2_col_sum <- colSums(df2)
  df2_norm <- df2
  for (i in 1:dim(df2)[2]){
    df2_norm[,i] <- df2[,i]/df2_col_sum[i]
  }
  # df3
  df3_col_sum <- colSums(df3)
  df3_norm <- df3
  for (i in 1:dim(df3)[2]){
    df3_norm[,i] <- df3[,i]/df3_col_sum[i]
  }
  
  # calculate average abundance and SD per taxon
  # get matrix taxon-average-SD
  # v1v2
  df1_taxon_aver_abund <- matrix(0,dim(df1)[1],4)
  tag_df1 <- paste0(word(file1,start=-2,sep='/'),'_',locat_trim(loc1))
  for (i in 1:dim(df1)[1]){
  df1_taxon_aver_abund[i,1] <- tag_df1
  df1_taxon_aver_abund[i,2] <- rownames(df1_norm)[i] 
  df1_taxon_aver_abund[i,3] <- mean(as.numeric(df1_norm[i,]))
  df1_taxon_aver_abund[i,4] <- sd(as.numeric(df1_norm[i,]))
  }
  # v3v4
  df2_taxon_aver_abund <- matrix(0,dim(df2)[1],4)
  tag_df2 <- paste0(word(file2,start=-2,sep='/'),'_',locat_trim(loc2))
  for (i in 1:dim(df2)[1]){
    df2_taxon_aver_abund[i,1] <- tag_df2
    df2_taxon_aver_abund[i,2] <- rownames(df2_norm)[i] 
    df2_taxon_aver_abund[i,3] <- mean(as.numeric(df2_norm[i,]))
    df2_taxon_aver_abund[i,4] <- sd(as.numeric(df2_norm[i,]))
  }
  # v5v6
  df3_taxon_aver_abund <- matrix(0,dim(df3)[1],4)
  tag_df3 <- paste0(word(file3,start=-2,sep='/'),'_',locat_trim(loc3))
  for (i in 1:dim(df3)[1]){
    df3_taxon_aver_abund[i,1] <- tag_df3
    df3_taxon_aver_abund[i,2] <- rownames(df3_norm)[i] 
    df3_taxon_aver_abund[i,3] <- mean(as.numeric(df3_norm[i,]))
    df3_taxon_aver_abund[i,4] <- sd(as.numeric(df3_norm[i,]))
  }
  # sort the matrices
  df1_taxon_aver_abund <- df1_taxon_aver_abund[order(as.numeric(df1_taxon_aver_abund[,3]),decreasing = TRUE),]
  df2_taxon_aver_abund <- df2_taxon_aver_abund[order(as.numeric(df2_taxon_aver_abund[,3]),decreasing = TRUE),]
  df3_taxon_aver_abund <- df3_taxon_aver_abund[order(as.numeric(df3_taxon_aver_abund[,3]),decreasing = TRUE),]
  # save first N elements
  N=10
  df1_taxon_aver_abund <- df1_taxon_aver_abund[1:N,]
  df2_taxon_aver_abund <- df2_taxon_aver_abund[1:N,]
  df3_taxon_aver_abund <- df3_taxon_aver_abund[1:N,]
  
  
  # construct the final table for plotting
  taxa_abund_table <- rbind(df1_taxon_aver_abund,df2_taxon_aver_abund,df3_taxon_aver_abund)
  colnames(taxa_abund_table) <- c('Samples','Taxon','Mean_abundance','SD')
  taxa_abund_table[,3] <- as.numeric(taxa_abund_table[,3])
  taxa_abund_table.df <- as.data.frame(taxa_abund_table)
  taxa_abund_table.df[,3] <- taxa_abund_table[,3]
  taxa_abund_table.df[,4] <- taxa_abund_table[,4]
  # setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Check_liter_review/Location_compar')
  # jpeg('tmp.jpg',width = 3800, height = 3800)
  ggplot(taxa_abund_table.df[,c(1,2,3)],aes(x = Samples, y = Mean_abundance, fill = Taxon))+
  geom_col(width = 0.5)+
    theme(legend.text=element_text(size=15),axis.text.x = element_text(size=100))
  # dev.off()

  # setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Check_liter_review/Location_compar')
  # jpeg('tmp.jpg',width = 3800, height = 3800)
  # ggplot(taxa_abund_table.df[,c(1,2,3)],aes(x = Samples, y = Mean_abundance, fill = Taxon))+
  #   geom_bar(stat = 'identity', na.rm = T, width = 0.1)
  # dev.off()
  # 
  
  
    
  # ggplot() +
  # geom_bar(data = Example,
  # aes(x = X_Axis, y = Percent, fill = Stack_Group), stat = 'identity', width = 0.5) 
  #  
}
