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

# There is a problem with this program regarding deleting low(1-9) read cunts. In v1 and v2 can be better. No time to solve. For Sunday!
plotBactFract <- function(file1, file2, files = T, sample_name = 'IPC293', 
                          locat='all', percent_threshold=0.01, remove.outliers=0.01, min.coverage = 10){
  source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Functions_general_non_enterotyping.R')
  library(stringr)
  # if we need to read the files. Otherwise, file1 and file2 are data frames
  if (files){
  df1 <- readQiimeSmart(input_file = file1, locat = locat, percent_threshold = percent_threshold, remove.outliers = remove.outliers,to.scale = FALSE)
  # a df1 but with read counts 
  # df1_helper <- readQiimeSmart(input_file = file1, locat = locat, percent_threshold = percent_threshold, remove.outliers = remove.outliers,to.scale = FALSE)
  df2 <- readQiimeSmart(input_file = file2, locat = locat, percent_threshold = 0, remove.outliers = 0,to.scale = FALSE)
  } else if (!is.data.frame(file1)){stop('file1 must be a data frame of taxa counts! If it is a file, use files = T argument.')
  } else if (!is.data.frame(file2)){stop('file1 must be a data frame of taxa counts! If it is a file, use files = T argument.')
  } else{
    df1 = file1
    df2 = file2
  }
  # X and Y axis names
  x_label <- paste0('Taxon fraction in ', word(file1,start = -2,sep='/'))
  y_label <- paste0('Taxon fraction in ', word(file2,start = -2,sep='/'))
  # get bacteria abundance values for sample_name in v1v2
  if(dim(df1[colnames(df1)%like%sample_name])[2]==1){
  # Define the threshold of min.coverage reads
  threshold.df1.log10 <- log(min.coverage/sum(df1[colnames(df1)%like%sample_name]),base = 10)
  threshold.df2.log10 <- log(min.coverage/sum(df2[colnames(df2)%like%sample_name]),base = 10)
  # and now can scale the numbers in df1 and df2
  df1 <- to.scale(df1)
  df2 <- to.scale(df2)
  x_coord <- df1[colnames(df1)%like%sample_name]
  x_coord_new <- x_coord[,1]
  names(x_coord_new) <- rownames(x_coord)
  }else if (dim(df1[colnames(df1)%like%sample_name])[2]>1){
    stop(as.character(paste0(sample_name,' is found > 1 times in ',file1,' file')))
   } else {
    stop(as.character(paste0(sample_name,' not found in the ',file1,' file')))
  }
  # check that there are no multiple sample_name's in df2
  if (dim(df2[colnames(df2)%like%sample_name])[2]>1){
    stop(as.character(paste0(sample_name,' is found > 1 times in ',file2,' file')))
  } else if (dim(df2[colnames(df2)%like%sample_name])[2]<1){
    stop(as.character(paste0(sample_name,' not found in the ',file2,' file')))
  }
  # compose a vector for df2 with abundances of the same taxa from df1[,sample_name] in df2
  df2_sample_name <- df2[colnames(df2)%like%sample_name]
  y_coord <- x_coord 
  black_label = integer(0) # this will contain numbers of genera(taxa) that are not present in df2 - to remove!
  for (i in 1:length(x_coord)){
    if(rownames(x_coord)[i] %in% rownames(df2)){
    y_coord[i] <- df2[rownames(df2) %in% names(x_coord)[i],colnames(df2)%like%sample_name]
    } else{
      black_label <- c(black_label,i)
      y_coord[i] = 0
    }
  }
  # x_coord <- as.vector(x_coord)
  # y_coord <- as.vector(y_coord)
  if(length(black_label)!=0){
  taxon_missing_df2 <- rownames(x_coord)[black_label]
  warning(paste0('The following taxa went missing in file2: ',str_c(taxon_missing_df2,sep = ';')))
  x_coord_new <- as.data.frame(x_coord[-black_label,])
  rownames(x_coord_new) <- rownames(x_coord)[-black_label]
  y_coord_new <- as.data.frame(y_coord[-black_label,])
  rownames(y_coord_new) <- rownames(x_coord)[-black_label]
  x_coord <- x_coord_new
  y_coord <- y_coord_new
  }

  x_coord_log <- log(x_coord,base = 10)
  x_coord_log[x_coord_log==-Inf] <- -8
  y_coord_log <- log(y_coord,base = 10)
  y_coord_log[y_coord_log==-Inf] <- -8
  
  # now, draw the plot
  data_ = as.data.frame(cbind(x_coord_log,y_coord_log))
  # data_inconsistent <- data_[data_$x_coord_log>-8,] 
  # data_inconsistent <- data_inconsistent[data_inconsistent$y_coord_log>-8,]
  data_inconsistent <- data_inconsistent[abs(data_inconsistent$x_coord_log-data_inconsistent$y_coord_log)>1,]
  #remove from the inconsistent samples the ones that are either -8,<threshold.df1.log10 or threshold.df1.log10<,-8
  too_low_reads1 <- data_[data_$x_coord==-8,]
  too_low_reads1 <- too_low_reads1[too_low_reads1$y_coord<y_coord_log]
  too_low_reads2 <- data_[data_$y_coord==-8,]
  too_low_reads2 <- too_low_reads1[too_low_reads1$x_coord<x_coord_log]
  data_inconsistent <- data_inconsistent[!(rownames(data_inconsistent) %in% rownames(too_low_reads1)),]
  data_inconsistent <- data_inconsistent[!(rownames(data_inconsistent) %in% rownames(too_low_reads2)),]
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


# How many reads per sample - make graphs
# INPUT:
# input_file = qiime output file (input_type='QIIME', po umolch)      
# OR a tab-delimited file with sample names in column 1 and total read counts in col. 2 (input_type = 'total.coverage')
# locat = location, e.g. '.IL', percent_threshold = 0, remove.outliers = 0 to take all reads and samples into account
# cumulative = FALSE for graphs showing number of samples with cov = +-(usredneniye/2) k reads
# usredneniye = size of the rolling coverage window, see also above, po umolch = 16.
# to.append = if TRUE, line is added to the existing plot (color = color, po umolch 'red')
# y_axis_length_factor = limits of y axis are from 0 to max(y)*y_axis_length_factor, po umolch = 1.1
# x_axis_length_factor = limits of y axis are from 0 to max(x)*x_axis_length_factor, po umolch = 1.03
# my.subset = a subset of samples to plot (a vector)
# adjust.factor = the value to multiply all numbers of reads per sample. Can be used e.g. when appending 2 graphs
# that have different total number of samples. Do not change it otherwise!
# OUTPUT:
# a plot of number of reads for each location as well as all locations together
exploreCoverage <- function(input_type='QIIME', input_file, locat = 'all', percent_threshold = 0, remove.outliers = 0, cumulative = TRUE,
                                  usredneniye=16,to.append=FALSE,color='black',y_axis_length_factor = 1.1, x_axis_length_factor = 1.03,
                            my.subset=NULL,adjust.factor=1){
  library(stringr)
  library(zoo)
  library(data.table)
  source(file.path(progr.dir,'/Functions_general_non_enterotyping.R'))
  if(input_type=='QIIME'){
  df1 <- readQiimeSmart(input_file = input_file, locat = locat, percent_threshold = percent_threshold, remove.outliers = remove.outliers,to.scale = FALSE)
  colnames(df1) <- gsub('.V[0-9][0-9]','',colnames(df1))
  # if we need to display a specific subset of our samples
  if(!is.null(my.subset)){
    if(!is.vector(my.subset)){geterrmessage('my.subset must be NULL or a vector')
    }else {
        df1 <- df1[colnames(df1) %in% my.subset]
      }
  }
  colSummi <- colSums(df1)
  } else {
    if(!is.null(my.subset)){
      warning('filtration has not been implemented yet for input_type = total.coverage!')
    }
    df1 <- read.table(input_file)
    # Here we need tp get rid of .V12 in df1.. but is difficult: factors
    # df1[,1] <- gsub('.V[0-9][0-9]','',df1[,1])
    # chose samples of the target location only
    if(locat!='all'){df1 <- df1[df1[,1] %like% locat,]
    }else{df1 <-df1[df1[,1] %like% '[0-9][0-9]',]} #We assume that samples have names with at least 2 digits.. to get rid of non-sample columns
    # colSummi: need to make this vector
    # again factors...
    colSummi <- rep(0,times=dim(df1)[1])
    for (i in 1:dim(df1)[1]){
    sample.curr <- as.character(df1[i,1])
    colSummi[i] <- as.numeric(levels(df1[df1[,1]==sample.curr,2])[df1[df1[,1]==sample.curr,2]])
    }
  }
  max.coverage <- ceiling(max(colSummi)/1000)*1000
  coverage <- seq(0,max.coverage,by=1000)
  numb.of.samples <- rep(0,times=length(coverage)) #values will be changed, do not worry! only length matters
  data <- cbind(coverage,numb.of.samples)
  if (adjust.factor!=1){print(paste0('Beware: numbers of samples multiplied by ',adjust.factor))}
  # here we calculate number of samples with coverage >=n
  for(i in 1:dim(data)[1]){
    data[i,2] <- adjust.factor*length(colSummi[colSummi>=data[i,1]]) # adjust the 'number of samples' to 1000 in silico
  }
  
  
  # for plots
  if(locat == 'all') {locat <- '.all'} #techical for plot subtitle
  N <- length(colSummi)
  if(!is.null(my.subset)){N <- length(colSummi[names(colSummi) %in% my.subset])}
  subtitle = paste0(word(input_file,start=-2,sep='/'),' ampl. coverage; location=',str_sub(locat,start = 2),', N=',N)
  y_lim <- c(0,max(data[,2])*y_axis_length_factor)
  x_lim <- c(min(data[,1]),max(data[,1])*x_axis_length_factor)
  if(cumulative){
    if(!to.append){
      plot(data, ylim =  y_lim, xlim =  x_lim, xlab = 'coverage >=, reads', ylab='Number of samples', 
         main = subtitle,cex.main = 3, cex.lab = 5, cex.axis = 3,col.main='grey',col=color)
      grid(nx=20,lty=6)
    }else{
      lines(data,col=color)
    }
    # title(sub = subtitle, cex.sub = 4, line = -5, col.sub = 'grey')
    
  } else{
    # get samples covered n 000 to (n+1) 000
    for(i in 0:(dim(data)[1]-2)){
      data[dim(data)[1]-i,2] <- data[(dim(data)[1]-i)-1,2] - data[dim(data)[1]-i,2]
    }
    data <- data[-1,]
    x_values <- seq(usredneniye/2,dim(data)[1]-usredneniye/2,by=1)
    y_values <- rollsum(data[,2],k=usredneniye)
    y_lim <- c(0,max(y_values)*y_axis_length_factor)
    x_lim <- c(min(x_values),max(x_values)*x_axis_length_factor)
    if(!to.append){
    plot(x=x_values, y=y_values,  ylim =  y_lim, xlim =  x_lim, ylab = 'Number of samples',
         xlab=paste0('Coverage, th. reads +- ',usredneniye/2), 
         main = subtitle, cex.main = 3, cex.lab = 5, cex.axis = 2,type = 'l',col.main='grey',col=color,xaxt='n')
    axis(1, at = seq(0,dim(data)[1]-usredneniye/2, by = 2),cex.axis=3)
    grid(ny=NULL,nx=50,lty=6)
    }else{
      lines(x=x_values, y=y_values,type='l',col=color)
    }
  }
 }


# Correlate number of reads per sample with enterotype stability for a single amplicon - location pair
corrCoverEnterotype <- function(qiime_file, locat = 'all', stable_enterotype_file){
  source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Functions_general_non_enterotyping.R')
  data <- readQiimeSmart(input_file = qiime_file, locat = locat, percent_threshold = 0, remove.outliers = 0,to.scale = FALSE)
  stable.enterotypes <- read.table(stable_enterotype_file)
  colSummi <- colSums(data)
  stable <- rep(F,times=dim(data)[2])
  colnames(data) <- gsub('.V[0-9][0-9]','',colnames(data))
  stable <- colnames(data) %in% stable.enterotypes[,1]
  # a data frame with col 2 = enterot stability, 3 = read depth for the individual
  result <- cbind(colnames(data),stable,colSummi)
  # check difference of the means by 2 methods (normality not checked here)
  wilc_res = signif(wilcox.test(as.numeric(result[result[,2]==T,3]),as.numeric(result[result[,2]==F,3]))$p.value,digits = 3)
  t_res = signif(t.test(as.numeric(result[result[,2]==T,3]),as.numeric(result[result[,2]==F,3]))$p.value,digits = 3)
  # get mean coverage by group
  mean.stable = round(mean(as.numeric(result[result[,2]==T,3])),digits = 0)
  mean.unstable = round(mean(as.numeric(result[result[,2]==F,3])),digits = 0)
  output <- c('Mann-Whitney_res' = wilc_res, 't.test_res' = t_res, 'mean_cov_stable' = mean.stable, 
             'mean_cov_unstable' = mean.unstable)
  return(output)
}

# Plot sample coverage vs enterotype
# INPUT:
# qiime_file = QIIME output with reads per taxon
# locat = location to compare
# stable_enterotype_file = file with enterotypes stable across locations
# OUTPUT:
# a plot of coverage density function (across samples) for both stable and unstable samples
plotCoverEnterotype <- function(qiime_file, locat = 'all', stable_enterotype_file){
  source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Functions_general_non_enterotyping.R')
  data <- readQiimeSmart(input_file = qiime_file, locat = locat, percent_threshold = 0, remove.outliers = 0,to.scale = FALSE)
  stable.enterotypes <- read.table(stable_enterotype_file)
  colSummi <- colSums(data)
  stable <- rep(F,times=dim(data)[2])
  colnames(data) <- gsub('.V[0-9][0-9]','',colnames(data))
  stable <- colnames(data) %in% stable.enterotypes[,1]
  # a data frame with col 2 = enterot stability, 3 = read depth for the individual
  result <- cbind(colnames(data),stable,colSummi)
  if(locat == 'all') {locat <- '.all'} #techical for plot subtitle
  main_ <- paste0(word(qiime_file,start=-2,sep='/'),', location: ',str_sub(locat,start = 2))
  ylim_max = 1.1*max(max(density(as.numeric(result[result[,2]==F,3]))$y),max(density(as.numeric(result[result[,2]==T,3]))$y))
  plot(density(as.numeric(result[result[,2]==F,3])),main=main_,ylab = 'Density',xlab = 'Microbiome coverage',
       xlim = c(0,max(as.numeric(result[,3]))), ylim = c(0,ylim_max), cex.main = 3, cex.lab = 3, cex.axis = 2,col='red')
  lines(density(as.numeric(result[result[,2]==T,3])))
  grid(nx=20,lty=6)
}


# How different are abilities of amplicons to catch taxa
# INPUT
# QIIME files file1 and file2
# sample.cov.filter = for file1 and file2, samples with at least sample.cov.filter[1] and sample.cov.filter[2]
# reads considered. Po umolch = 10k reads for each file
# p.value = 0.01 - p-value for returning only taxa that are described differently across amplicons. Bonferoni
# correction applied in the program (p.value = 0.01/N e.g.)
# OUTPUT:
# a data frame of taxa that have statistically signif (t-test) different abundance across the 2 amplicons.
# columns: average abindance df1, df2, difference, sd1, sd2, non-corrected p-value
# rows: taxa
checkTaxaAcrossAmpl <- function(file1, file2, locat='all', log.scale = FALSE, percent_threshold=0.01, 
                                remove.outliers=0.01, sample.cov.filter=c(10000,10000),p.value = 0.01,
                                select.taxa = 't-test'){
  source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Functions_general_non_enterotyping.R')
  library(stringr)
  if(!is.vector(sample.cov.filter)){geterrmessage('sample.cov.filter must be a vector of 2 elements')}
  df1 <- readQiimeSmart(input_file = file1, locat = locat, percent_threshold = percent_threshold, 
                        remove.outliers = remove.outliers, sample.cov.filter = sample.cov.filter[1], to.scale = TRUE)
  df2 <- readQiimeSmart(input_file = file2, locat = locat, percent_threshold = percent_threshold, 
                        remove.outliers = remove.outliers, sample.cov.filter = sample.cov.filter[2], to.scale = TRUE)
  #in df1 and df2, remove .V12 etc.
  colnames(df1) <- gsub('.V[0-9][0-9]','',colnames(df1))
  colnames(df2) <- gsub('.V[0-9][0-9]','',colnames(df2))
  # get samples and taxa present in both df and df2
  common.samples <- intersect(colnames(df1),colnames(df2))
  common.taxa <- intersect(rownames(df1),rownames(df2))
  test <- numeric(length(common.taxa))
  # define mean and SD taxon fraction vectors to be filled in
  taxa.fraction.df1.mean <- numeric(length(common.taxa))
  taxa.fraction.df1.sd <- numeric(length(common.taxa))
  taxa.fraction.df2.mean <- numeric(length(common.taxa))
  taxa.fraction.df2.sd <- numeric(length(common.taxa))
  # log transformation if required
  if(log.scale){
    # the abindances <10^-8 set to 10^-8 to avoid -Inf
    df1[df1<0.00000001] <- 0.00000001
    df2[df2<0.00000001] <- 0.00000001
    df1 <- log10(df1)
    df2 <- log10(df2)
  }
  for (i in 1:length(common.taxa)){
    test[i] <- t.test(as.numeric(subset(df1,rownames(df1) %in% common.taxa[i])),as.numeric(subset(df2,rownames(df2) %in% common.taxa[i])))$p.value
    names(test)[i] <- common.taxa[i]
    taxa.fraction.df1.mean[i] <- mean(as.numeric(subset(df1,rownames(df1) %in% common.taxa[i])))
    taxa.fraction.df1.sd[i] <- sd(as.numeric(subset(df1,rownames(df1) %in% common.taxa[i])))
    taxa.fraction.df2.mean[i] <- mean(as.numeric(subset(df2,rownames(df2) %in% common.taxa[i])))
    taxa.fraction.df2.sd[i] <- sd(as.numeric(subset(df2,rownames(df2) %in% common.taxa[i])))
    names(taxa.fraction.df1.mean)[i] <- common.taxa[i]
    names(taxa.fraction.df2.mean)[i] <- common.taxa[i]
    names(taxa.fraction.df1.sd)[i] <- common.taxa[i]
    names(taxa.fraction.df2.sd)[i] <- common.taxa[i]
  }
  difference <- taxa.fraction.df2.mean-taxa.fraction.df1.mean
  taxa.fraction.info <- cbind(taxa.fraction.df1.mean,taxa.fraction.df2.mean,difference,taxa.fraction.df1.sd,taxa.fraction.df2.sd,test)
  if (select.taxa == 't-test'){
  result <- taxa.fraction.info[taxa.fraction.info[,6]<p.value/length(common.taxa),]
  } else if (select.taxa == 'all'){
    result <- taxa.fraction.info
  } else{
    stop('select.taxa must be t-test or all')
  }
  return(result)
}