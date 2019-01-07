# Get stable enterotypes across V12 V34 and V56, V12-V34, V34-V56 and V12-V56 and the unstable rest
# INPUT
# 3 files = sample-enterotype
# correspondence = NA is usually okay. In case we get inconsistency across locations, we get warning here and
# correspondence should be made manually
# OUTPUT
# 5 data frames of stable enterotypes across V12 V34 and V56, V12-V34, V34-V56 and V12-V56 and the unstable rest

getStabAcrAmpl <- function(file1, file2, file3, correspondence=NA){
  library(data.table)
  source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Enterotyping_compar_functions_v5.R')
  source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/GetStablyEnterotypedSamples.R')

  # read the files 
  df_1 <- read.table(file1)
  df_2 <- read.table(file2)
  df_3 <- read.table(file3)
  
  if(is.na(correspondence)){
  # design the correspondence matrix
  corr_13 <- CompEntPred.pair(df_1,df_3)$Optimal_cluster_correspondences
  corr_12 <-CompEntPred.pair(df_1,df_2)$Optimal_cluster_correspondences
  corr_23 <-CompEntPred.pair(df_2,df_3)$Optimal_cluster_correspondences
  
# need to correct the row order of corr. As it is inversed in the output of CompEntPred.pair
  corr_13 <- corr_13[c(2,1),]
  corr_12 <- corr_12[c(2,1),]
  corr_23 <- corr_23[c(2,1),]
  
  # we assume that v5v6 (df3) is the best amplicon so it most likely has drivers of Bacteroides and Prevotella
  # this is df1 driver of ent 1.   
  c1_1 <- corr_13[1,corr_13[2,]%like%'g__Bacteroides']
  # this is df3 driver of ent 1.   
  c1_3 <- corr_13[2,corr_13[2,]%like%'g__Bacteroides']
  c2_1 <- corr_13[1,corr_13[2,]%like%'g__Prevotella']
  c2_3 <- corr_13[2,corr_13[2,]%like%'g__Prevotella']
  c3_1 <- levels(df_1[,2])[!levels(df_1[,2])%in%c(c1_1,c2_1)]
  c3_3 <- levels(df_3[,2])[!levels(df_3[,2])%in%c(c1_3,c2_3)]
  # and now finalize the column 2 of the correspondence matrix (v3v4) from corr_23
  c1_2 <- corr_23[1,corr_23[2,]==c1_3]
  c2_2 <- corr_23[1,corr_23[2,]==c2_3]
  c3_2 <- corr_23[1,corr_23[2,]==c3_3]
  # now check that the correspondence 1-3 and 2-3 are in aggreement with 1-2
  if (corr_12[2,corr_12[1,]==c1_1]!=c1_2){
    warning('enterotype 1 drivers are not conserved across v12-v56, v34-v56 and v12-v34 amplicon pairs! Need to 
            set driver correspondence matrix manually! see kokon')
  }
  if (corr_12[2,corr_12[1,]==c2_1]!=c2_2){
    warning('enterotype 2 drivers are not conserved across v12-v56, v34-v56 and v12-v34 amplicon pairs! Need to 
            set driver correspondence matrix manually! see kokon')
  }
  if (corr_12[2,corr_12[1,]==c3_1]!=c3_2){
    warning('enterotype 3 drivers are not conserved across v12-v56, v34-v56 and v12-v34 amplicon pairs! Need to 
            set driver correspondence matrix manually! see kokon')
  }
  correspondence <- matrix(c(c1_1,c2_1,c3_1,c1_2,c2_2,c3_2,c1_3,c2_3,c3_3),3,3)
  colnames(correspondence) <- c('v12','v34','v56')
  rownames(correspondence) <- c('ent1','ent2','ent3')
  }
  
  
  result.v12v34v56 <- getStableSamples(df1=df_1,df2=df_2,df3=df_3,correspondence)

    # make correspondences for v12-v34-v12, v12-v56-v12, v34-v56-v34 and v12-v12-v12
  correspondence.1234 <- as.matrix(cbind(correspondence[,1],correspondence[,2],correspondence[,1]))
  correspondence.1256 <- as.matrix(cbind(correspondence[,1],correspondence[,3],correspondence[,1]))
  correspondence.3456 <- as.matrix(cbind(correspondence[,2],correspondence[,3],correspondence[,2]))
  correspondence.1212 <- as.matrix(cbind(correspondence[,1],correspondence[,1],correspondence[,1]))
  
  # make lists of samples with enterotypes conserved per 2 amplicons
  result.v12v34 <- getStableSamples(df1=df_1,df2=df_2,df3=df_1,correspondence.1234)
  result.v12v56 <- getStableSamples(df1=df_1,df2=df_3,df3=df_1,correspondence.1256)
  result.v34v56 <- getStableSamples(df1=df_2,df2=df_3,df3=df_2,correspondence.3456)
  # and now substarct the samples with enterotypes conserved per all 3 amplicons
  result.v12v34 <- result.v12v34[!result.v12v34[,1]%in%result.v12v34v56[,1],]
  result.v12v56 <- result.v12v56[!result.v12v56[,1]%in%result.v12v34v56[,1],]
  result.v34v56 <- result.v34v56[!result.v34v56[,1]%in%result.v12v34v56[,1],]
  
  # make list of samples with enterotypes not conserved across amplicons at all
  result.nonconserved <- getStableSamples(df1=df_1,df2=df_1,df3=df_1,correspondence.1212)
  result.chimera.1234.1256.3456.123456 <- rbind(result.v12v34,result.v12v56,result.v34v56,result.v12v34v56)
  result.nonconserved <- result.nonconserved[!result.nonconserved[,1]%in%result.chimera.1234.1256.3456.123456[,1],]
  #transform matrices into data frames 
  result.v12v34v56 <- as.data.frame(result.v12v34v56)
  result.v12v34 <- as.data.frame(result.v12v34)
  result.v12v56 <- as.data.frame(result.v12v56)
  result.v34v56 <- as.data.frame(result.v34v56)
  result.nonconserved <- as.data.frame(result.nonconserved)
  # get the final result
  result.total <- list(result.v12v34v56,result.v12v34,result.v12v56,result.v34v56,result.nonconserved,correspondence)
  names(result.total) <- c('result.v12v34v56','result.v12v34','result.v12v56','result.v34v56','result.nonconserved','driver.correspondence')
  return(result.total)
}

