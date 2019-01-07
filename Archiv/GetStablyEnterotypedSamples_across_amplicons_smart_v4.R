# v4. Not only 3 enterotypes, but also 2 enterotypes possible
# Get stable enterotypes across V12 V34 and V56, V12-V34, V34-V56 and V12-V56 and the unstable rest
# INPUT
# 3 files = sample-enterotype
# correspondence = NA is usually okay. In case we get inconsistency across locations, we get warning here and
# correspondence should be made manually
# OUTPUT
# 5 data frames of stable enterotypes across V12 V34 and V56, V12-V34, V34-V56 and V12-V56 and the unstable rest

getStabAcrAmpl <- function(file1=file.path(home.dir,'Results/Enterotyping/2018_11_30_enterotyping/Sample-enterotype/V1V2_2e_RE_sample-enterotype.csv'),
                           file2=file.path(home.dir,'Results/Enterotyping/2018_11_30_enterotyping/Sample-enterotype/V3V4_2e_RE_sample-enterotype.csv'), 
                           file3=file.path(home.dir,'Results/Enterotyping/2018_11_30_enterotyping/Sample-enterotype/V5V6_2e_RE_sample-enterotype.csv'), 
                           correspondence=NA, files=TRUE){
  library(data.table)
  source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Enterotyping_compar_functions_v5.R')
  source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/GetStablyEnterotypedSamples_v2.R')


  # read the files 
  if(files){
    df_1 <- read.table(file1)
    df_2 <- read.table(file2)
    df_3 <- read.table(file3)
  }else{
    df_1 <- file1
    df_2 <- file2
    df_3 <- file3
  }
  # check that we have 2 or 3 enterotypes
  if ((length(levels(df_1[,2]))!=length(levels(df_2[,2])))|(length(levels(df_1[,2]))!=length(levels(df_3[,2])))){
    stop('Number of enterotypes not the same across df1, df2 and df3')
    geterrmessage()
  }
  number.enterotypes <- length(levels(df_1[,2]))
  if(number.enterotypes>=4){
    stop('The program is designed for analysis of 2 or 3 enterotypes only!')
    geterrmessage()
  }
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
  if(number.enterotypes==3){
  c3_1 <- levels(df_1[,2])[!levels(df_1[,2])%in%c(c1_1,c2_1)]
  c3_3 <- levels(df_3[,2])[!levels(df_3[,2])%in%c(c1_3,c2_3)]
  }
  # and now finalize the column 2 of the correspondence matrix (v3v4) from corr_23
  c1_2 <- corr_23[1,corr_23[2,]==c1_3]
  c2_2 <- corr_23[1,corr_23[2,]==c2_3]
  if(number.enterotypes==3){c3_2 <- corr_23[1,corr_23[2,]==c3_3]}
  # now check that the correspondence 1-3 and 2-3 are in aggreement with 1-2
  if (corr_12[2,corr_12[1,]==c1_1]!=c1_2){
    warning('enterotype 1 drivers are not conserved across v12-v56, v34-v56 and v12-v34 amplicon pairs! Need to 
            set driver correspondence matrix manually! see kokon')
  }
  if (corr_12[2,corr_12[1,]==c2_1]!=c2_2){
    warning('enterotype 2 drivers are not conserved across v12-v56, v34-v56 and v12-v34 amplicon pairs! Need to 
            set driver correspondence matrix manually! see kokon')
  }
  if(number.enterotypes==3){
  if (corr_12[2,corr_12[1,]==c3_1]!=c3_2){
    warning('enterotype 3 drivers are not conserved across v12-v56, v34-v56 and v12-v34 amplicon pairs! Need to 
            set driver correspondence matrix manually! see kokon')
  }
  }
  if(number.enterotypes==3){
  correspondence <- matrix(c(c1_1,c2_1,c3_1,c1_2,c2_2,c3_2,c1_3,c2_3,c3_3),3,3)
  colnames(correspondence) <- c('v12','v34','v56')
  rownames(correspondence) <- c('ent1','ent2','ent3')
  } else if(number.enterotypes==2){
    correspondence <- matrix(c(c1_1,c2_1,c1_2,c2_2,c1_3,c2_3),2,3)
    colnames(correspondence) <- c('v12','v34','v56')
    rownames(correspondence) <- c('ent1','ent2')
  }
  }
  
  result.v12v34v56 <- getStableSamples(df1=df_1,df2=df_2,df3=df_3,correspondence=correspondence,
                                       number.enterotypes=number.enterotypes)[[1]]
  unstable.v12v34v56 <- getStableSamples(df1=df_1,df2=df_2,df3=df_3,correspondence=correspondence,
                                         number.enterotypes=number.enterotypes)[[2]]
    # make correspondences for v12-v34-v12, v12-v56-v12, v34-v56-v34 and v12-v12-v12
  correspondence.1234 <- as.matrix(cbind(correspondence[,1],correspondence[,2],correspondence[,1]))
  correspondence.1256 <- as.matrix(cbind(correspondence[,1],correspondence[,3],correspondence[,1]))
  correspondence.3456 <- as.matrix(cbind(correspondence[,2],correspondence[,3],correspondence[,2]))
  correspondence.1212 <- as.matrix(cbind(correspondence[,1],correspondence[,1],correspondence[,1]))
  
  # make lists of samples with enterotypes conserved per 2 amplicons
  result.v12v34 <- getStableSamples(df1=df_1,df2=df_2,df3=df_1,correspondence.1234,number.enterotypes)[[1]]
  unstable.v12v34 <- getStableSamples(df1=df_1,df2=df_2,df3=df_1,correspondence.1234,number.enterotypes)[[2]]
  result.v12v56 <- getStableSamples(df1=df_1,df2=df_3,df3=df_1,correspondence.1256,number.enterotypes)[[1]]
  unstable.v12v56 <- getStableSamples(df1=df_1,df2=df_3,df3=df_1,correspondence.1256,number.enterotypes)[[2]]
  result.v34v56 <- getStableSamples(df1=df_2,df2=df_3,df3=df_2,correspondence.3456,number.enterotypes)[[1]]
  unstable.v34v56 <- getStableSamples(df1=df_2,df2=df_3,df3=df_2,correspondence.3456,number.enterotypes)[[2]]
  
  # and now substarct the samples with enterotypes conserved per all 3 amplicons
  result.v12v34 <- result.v12v34[!result.v12v34[,1]%in%result.v12v34v56[,1],]
  result.v12v56 <- result.v12v56[!result.v12v56[,1]%in%result.v12v34v56[,1],]
  result.v34v56 <- result.v34v56[!result.v34v56[,1]%in%result.v12v34v56[,1],]
  
  # # make list of samples with enterotypes not conserved across amplicons at all
  # result.nonconserved <- getStableSamples(df1=df_1,df2=df_1,df3=df_1,correspondence.1212,number.enterotypes)[[1]]
  # result.chimera.1234.1256.3456.123456 <- rbind(result.v12v34,result.v12v56,result.v34v56,result.v12v34v56)
  # result.nonconserved <- result.nonconserved[!result.nonconserved[,1]%in%result.chimera.1234.1256.3456.123456[,1],]
  
  #transform matrices into data frames 
  result.v12v34v56 <- as.data.frame(result.v12v34v56)
  result.v12v34 <- as.data.frame(result.v12v34)
  result.v12v56 <- as.data.frame(result.v12v56)
  result.v34v56 <- as.data.frame(result.v34v56)
  # result.nonconserved <- as.data.frame(result.nonconserved)
  # get the final result
  result.total <- list(result.v12v34v56,result.v12v34,result.v12v56,result.v34v56,unstable.v12v34v56,unstable.v12v34,
                       unstable.v12v56,unstable.v34v56,correspondence)
  names(result.total) <- c('stable.v12v34v56','stable.v12v34','stable.v12v56','stable.v34v56','unstable.v12v34v56','unstable.v12v34',
                           'unstable.v12v56','unstable.v34v56','driver.correspondence')
  return(result.total)
}

# this program gives final list of samples and enterotypes conserved for 2 enterotypes
# INPUT:
# 3 files = sample-enterotype
# correspondence = NA is usually okay. In case we get inconsistency across locations, we get warning here and
# correspondence should be made manually
# logical arguments:
# stable.three = if we take samples stable across all 3 amplicons (=TRUE)
# unstable.three= if we take  samples that are unstable across 3 amplicons. If fasle, these samples are deleted 
# from output even if stable across 2 amplicons (=FALSE)
# stable.two = if we take samples stable across 2 amplicons (=TRUE)
# stable.one - if we take samples 'stable' across only 1 amplicon (which is all of them with known enterotype)
getStabAcrAmpl.final.2ent <- function(file1, file2, file3, stable.three=T, stable.two=T, stable.one=T, unstable.three = F,
                                      unstable.two = F, correspondence=NA, files=TRUE)
  stable.total <- data.frame("stable.samples"=character(),"enterotypes"=character(),"enterotypes.short"=character())
  stable.three <- data.frame("stable.samples"=character(),"enterotypes"=character(),"enterotypes.short"=character())
  unstable.three <- data.frame("stable.samples"=character(),"enterotypes"=character(),"enterotypes.short"=character())
  stable.two <- data.frame("stable.samples"=character(),"enterotypes"=character(),"enterotypes.short"=character())
  unstable.two <- data.frame("stable.samples"=character(),"enterotypes"=character(),"enterotypes.short"=character())
  only.one.ampl.known.ent <- data.frame("stable.samples"=character(),"enterotypes"=character(),"enterotypes.short"=character())
  stability.across.ampl <- getStabAcrAmpl(file1, file2, file3)
  if(stable.three){result.three <- rbind(result,stability.across.ampl$stable.v12v34v56)}
  if(!unstable.three){result.unstable.three <- stability.across.ampl$unstable.v12v34v56}
  if(stable.two){
    result <- rbind(result,stability.across.ampl$stable.v12v34)
    result <- rbind(result,stability.across.ampl$stable.v12v56)
    result <- rbind(result,stability.across.ampl$stable.v34v56)
  }

  if(!unstable.three){
    result
  }
  
  