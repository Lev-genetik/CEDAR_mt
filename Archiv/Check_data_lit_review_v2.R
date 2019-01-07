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
entrShannon <- function(input_file, locat = 'all', percent_threshold = 0.01,remove.outliers = 0.01){
  source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Functions_general_non_enterotyping.R')
  library(dplyr)
  data <- readQiimeSmart(input_file = input_file, locat = locat, percent_threshold = percent_threshold, remove.outliers = remove.outliers)
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
# INPUT:
# (need to open connection before running the program and close after!)
# 3 files, 3 locations (e.g. file=file2=file3 (v1v2), loc1 = '.IL', loc2 = '.TR',loc3 = '.RE')
# percent_threshold - to remove noise
# remove.outliers - p-value threshold to remove outliers farest from the other samples
# OUTPUT:
# stacked barplot with lenend 
compTaxa <- function(file1,loc1,file2,loc2,file3,loc3,percent_threshold = 0.01,remove.outliers = 0.01){
  library(stringr)
  library(ggplot2)
  source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Functions_general_non_enterotyping.R')
  df1 <- readQiimeSmart(input_file = file1, locat = loc1, percent_threshold = percent_threshold, remove.outliers = remove.outliers)
  df2 <- readQiimeSmart(input_file = file2, locat = loc2, percent_threshold = percent_threshold, remove.outliers = remove.outliers)
  df3 <- readQiimeSmart(input_file = file3, locat = loc3, percent_threshold = percent_threshold, remove.outliers = remove.outliers)
  # appendix - this is for labels only
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
  ggplot(taxa_abund_table.df[,c(1,2,3)],aes(x = Samples, y = Mean_abundance, fill = Taxon))+
  geom_col(width = 0.5)+
    theme(legend.text=element_text(size=15),axis.text.x = element_text(size=100))
 
  # pomoika nizhe:
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


# Plot [v1v2] bacteria abundance vs [v3v4] bacteria abundance for a given individual
# AIM:
# make a plot of taxa abundance difference across amplicons and label the taxa that have >10-fold difference
# log scale
# Note: do not forget to open (before running the program) and close (after) the file of output
# Note: this program has a small bug. The names of outlying taxa are correct (I checked), but the position does not
# correspond to the outlying points
# INPUT: 
# file1, file2 = files to compare.
# percent_threshold = noise filtration percentage for file1 - to get rid of too rare taxa. 
# For file2, we have no noise filtration because need all taxa
# locat = location to see - e.g. '.IL'
# sample_name = name of the sample to use for comparison
# remove.outliers = samples standing too far (p-value < remove.outliers) are removed, po umolch 0.01
# OUTPUT:
# a plot of abundance_amplicon1 vs abundance_amplicon2 in log10 coordinates
plotBactFract <- function(file1, file2, sample_name = 'IPC293', 
                          locat='all', percent_threshold=0.01, remove.outliers=0.01, min.coverage = 10){
  source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Functions_general_non_enterotyping.R')
  library(stringr)
  df1 <- readQiimeSmart(input_file = file1, locat = locat, percent_threshold = percent_threshold, remove.outliers = remove.outliers,to.scale = TRUE)
  # a df1 but with read counts 
  df1_helper <- readQiimeSmart(input_file = file1, locat = locat, percent_threshold = percent_threshold, remove.outliers = remove.outliers,to.scale = FALSE)
  df2 <- readQiimeSmart(input_file = file2, locat = locat, percent_threshold = 0, remove.outliers = 0,to.scale = TRUE)
  # X and Y axis names
  x_label <- paste0('Taxon fraction in ', word(file1,start = -2,sep='/'))
  y_label <- paste0('Taxon fraction in ', word(file2,start = -2,sep='/'))
  # get bacteria abundance values for sample_name in v1v2
  if(dim(df1[colnames(df1)%like%sample_name])[2]==1){
  x_coord <- df1[colnames(df1)%like%sample_name]
  x_coord_new <- x_coord[,1]
  names(x_coord_new) <- rownames(x_coord)
  # here, the x_coord_helper is a technical vector, needed to delete taxa with low abundance < min.coverage
  x_coord_helper <- df1_helper[colnames(df1)%like%sample_name]
  x_coord_new_helper <- x_coord_helper[,1]
  # now remove from analysis taxa with <n read counts for our sample
  x_coord <- x_coord_new[which(x_coord_new_helper>min.coverage)]
  }else if (dim(df1[colnames(df1)%like%sample_name])[2]>1){
    stop(as.character(paste0(sample_name,' is found > 1 times in ',file1,' file')))
    # warning('multiple sample matches in file1! The first one taken')
    # x_coord <- df1[colnames(df1)%like%sample_name]
    # x_coord <- x_coord[,1]
  } else {
    stop(as.character(paste0(sample_name,' not found in the ',file1,' file')))
  }
  # check that there are no multiple sample_name's in df2
  if (dim(df2[colnames(df2)%like%sample_name])[2]>1){
    stop(as.character(paste0(sample_name,' is found > 1 times in ',file2,' file')))
    # warning('multiple sample matches in file1! The first one taken')
    # x_coord <- df1[colnames(df1)%like%sample_name]
    # x_coord <- x_coord[,1]
  } else if (dim(df2[colnames(df2)%like%sample_name])[2]<1){
    stop(as.character(paste0(sample_name,' not found in the ',file2,' file')))
  }
  # compose a vector for df2 with abundances of the same taxa from df1[,sample_name] in df2
  df2_sample_name <- df2[colnames(df2)%like%sample_name]
  y_coord <- x_coord 
  black_label = integer(0) # this will contain numbers of genera(taxa) that are not present in df2 - to remove!
  for (i in 1:length(x_coord)){
    if(names(x_coord)[i] %in% rownames(df2)){
    y_coord[i] <- df2[rownames(df2) %in% names(x_coord)[i],colnames(df2)%like%sample_name]
    } else{
      black_label <- c(black_label,i)
      y_coord[i] = 0
    }
  }
  if(length(black_label)!=0){
  taxon_missing_df2 <- names(x_coord[black_label])
  warning(paste0('The following taxa went missing in file2: ',str_c(taxon_missing_df2,sep = ';')))
  x_coord <- x_coord[-black_label]
  y_coord <- y_coord[-black_label]
  }

  x_coord_log <- log(x_coord,base = 10)
  x_coord_log[x_coord_log==-Inf] <- -8
  y_coord_log <- log(y_coord,base = 10)
  y_coord_log[y_coord_log==-Inf] <- -8
  
  # now, draw the plot
  data_ = as.data.frame(cbind(x_coord_log,y_coord_log))
  data_inconsistent <- data_[data_$x_coord_log>-8,] 
  data_inconsistent <- data_inconsistent[data_inconsistent$y_coord_log>-8,]
  data_inconsistent <- data_inconsistent[abs(data_inconsistent$x_coord_log-data_inconsistent$y_coord_log)>1,]
  #logical vector if the taxon in data_ has inconsistency across amplicons
  in_data_inconsistent = rownames(data_) %in% rownames(data_inconsistent) 
  # a set of specific labels for the plot: only family and genus
  labels1 = row.names(data_[which(in_data_inconsistent),])
  labels1 <- word(labels1,start=-2,end=-1,sep=';')
  # draw the plot
  par(mar=c(7,7,7,7))
  plot(data_,main = paste0('Bacteria fraction agreement for ', sample_name), xlim = c(-8,0), ylim = c(-8,0), 
       xlab = x_label, ylab = y_label, cex.main = 2.5, 
       cex.lab=3,cex.axis=1.5)
  title(sub = expression('log'[10]*'scale'),cex.sub=2,line=6)
  
  lines(rbind(c(1,1),c(-13,-13)),col = 'red')
  lines(rbind(c(0,1),c(-14,-13)),col = 'blue')
  lines(rbind(c(1,0),c(-13,-14)),col = 'blue')
  with(data_,text(x_coord_log[which(in_data_inconsistent)]~y_coord_log[which(in_data_inconsistent)],labels = labels1, pos = 1,cex=1,col='red'))
  # with(data_,text(x_coord_log~y_coord_log,labels = row.names(data_), pos = 4)) #all labels
  # with(data_inconsistent,text(x_coord_log[in_data_inconsistent]~y_coord_log[in_data_inconsistent],labels = row.names(data_[in_data_inconsistent,]), pos = 3))
  # with(data_inconsistent,text(x_coord_log[1:dim(data_inconsistent)[1]]~y_coord_log[1:dim(data_inconsistent)[1]],labels = row.names(data_[1:dim(data_inconsistent)[1],]), pos = 3,cex=0.4,col='red'))
  
  

}