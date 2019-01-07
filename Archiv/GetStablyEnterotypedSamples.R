# get table sample-enterotype for only samples that cluster to the same clusters (3 enterotypes)
# df1, df2, df3 = data frames of sample-enterotype
# correspondence = data frame of correspondence of Drivers to each other
# row1: ent1 Unique Driver for v12, v34 and v56
# row2: ent2 Unique Driver for v12, v34 and v56
# row3: ent3 Unique Driver for v12, v34 and v56

getStableSamples <- function(df1,df2,df3,correspondence){
   if((dim(df1)[2]!=2)|(dim(df2)[2]!=2)|(dim(df3)[2]!=2)){
   stop()
   geterrmessage('Error: number of columns in the input data frames should be 2!')
 } 
  if(dim(df1)[1]!=dim(df2)[1]){
    print('Warning: number of samples in data frames differs')
  }
  # delete '.V34' etc
  df1[,1] <- gsub('.V[1-6][1-6]','',df1[,1])
  df2[,1] <- gsub('.V[1-6][1-6]','',df2[,1])
  df3[,1] <- gsub('.V[1-6][1-6]','',df3[,1])
  
  # check that df2 and df3 has all samples from df1
  if((length(intersect(df1[,1],df2[,1]))!=dim(df1)[1])|(length(intersect(df1[,1],df3[,1]!=dim(df1)[1])))){
    print('Warning: not all samples from df1 present in df2 and/or df3')
  }
  # add 3rd row to df1,df2 and df3: universal cluster number (1,2 or 3)
    df1[df1[,2]==correspondence[1,1],3] <- 'a'
    df2[df2[,2]==correspondence[1,2],3] <- 'a'
    df3[df3[,2]==correspondence[1,3],3] <- 'a'
    df1[df1[,2]==correspondence[2,1],3] <- 'b'
    df2[df2[,2]==correspondence[2,2],3] <- 'b'
    df3[df3[,2]==correspondence[2,3],3] <- 'b'
    df1[df1[,2]==correspondence[3,1],3] <- 'c'
    df2[df2[,2]==correspondence[3,2],3] <- 'c'
    df3[df3[,2]==correspondence[3,3],3] <- 'c'
    # print(df1[df1[,1]=='IPC147.RE',])
    # print(df2[df2[,1]=='IPC147.RE',])
    # print(df3[df3[,1]=='IPC147.RE',])

    # get the samples where enterotypes are the same for all 3 16S parts 
    # step1.Oobtain vectors of ent1, ent2 and ent3 for all 3 data frames
    df1.ent1 <- df1[df1[,3]=='a',1]
    df1.ent2 <- df1[df1[,3]=='b',1]
    df1.ent3 <- df1[df1[,3]=='c',1]
    df2.ent1 <- df2[df2[,3]=='a',1]
    df2.ent2 <- df2[df2[,3]=='b',1]
    df2.ent3 <- df2[df2[,3]=='c',1]    
    df3.ent1 <- df3[df3[,3]=='a',1]
    df3.ent2 <- df3[df3[,3]=='b',1]
    df3.ent3 <- df3[df3[,3]=='c',1]   
    
       
    # step2. get samples that have ent1 for all 3 16S regions
    ent1 <- intersect(df1.ent1,df2.ent1)
    ent1 <- intersect(ent1,df3.ent1)
    ent2 <- intersect(df1.ent2,df2.ent2)
    ent2 <- intersect(ent2,df3.ent2)
    ent3 <- intersect(df1.ent3,df2.ent3)
    ent3 <- intersect(ent3,df3.ent3)
 
    # get the resultant data frame
    stable.samples <- c(ent1,ent2,ent3)
    ent1.name <- paste(correspondence[1,],collapse='|')
    ent2.name <- paste(correspondence[2,],collapse='|')
    ent3.name <- paste(correspondence[3,],collapse='|')
    

    enterotypes.tmp <- rep(list(character(3)),7)
    enterotypes.tmp[[1]] <- c(rep(ent1.name,length(ent1)))
    enterotypes.tmp[[2]] <- c(rep(ent2.name,length(ent2)))
    enterotypes.tmp[[3]] <- c(rep(ent3.name,length(ent3)))
 
    # Now make enterotype names short - Bacteroides, Prevotella and 'Other'.
    # first, check that they are such. If yes, we make enterotypes - Bacteroides, Prevotella and enterotype_3
    # if no, we make enterotypes - enterotype_1, enterotype_2 and enterotype_3
    ent.names.bacteroides <- grep('g__Bacteroides',c(ent1.name,ent2.name,ent3.name))
    ent.names.prevotella <- grep('g__Prevotella',c(ent1.name,ent2.name,ent3.name))
    if(length(ent.names.bacteroides)+length(ent.names.prevotella)==2){
    ent.names.other <- intersect(grep('g__Prevotella',c(ent1.name,ent2.name,ent3.name),invert=T),grep('g__Bacteroides',c(ent1.name,ent2.name,ent3.name),invert=T))
    
    enterotypes.tmp.short <- rep(list(character(3)),3)
    enterotypes.tmp.short[[ent.names.bacteroides]] <- c(rep('Bacteroides',length(enterotypes.tmp[[ent.names.bacteroides]])))
    enterotypes.tmp.short[[ent.names.prevotella]] <- c(rep('Prevotella',length(enterotypes.tmp[[ent.names.prevotella]])))
    enterotypes.tmp.short[[ent.names.other]] <- c(rep('enterotype_3',length(enterotypes.tmp[[ent.names.other]])))
    } else {
     warning('the first two enterotypes are not Bacteroides and Prevotella. The short names assigned are enterotype_1, enterotype_2 and enterotype_3')      
      enterotypes.tmp.short <- rep(list(character(3)),3)
      enterotypes.tmp.short[[1]] <- c(rep('enterotype_1',length(enterotypes.tmp[[1]])))
      enterotypes.tmp.short[[2]] <- c(rep('enterotype_2',length(enterotypes.tmp[[2]])))
      enterotypes.tmp.short[[3]] <- c(rep('enterotype_3',length(enterotypes.tmp[[3]])))
    }
    enterotypes <- c(enterotypes.tmp[[1]],enterotypes.tmp[[2]],enterotypes.tmp[[3]])
    enterotypes.short <- c(enterotypes.tmp.short[[1]],enterotypes.tmp.short[[2]],enterotypes.tmp.short[[3]])

    if(length(enterotypes)!=length(stable.samples)){
      stop()
      geterrmessage('Internal error #1 - length of enterotypes and samples not the same!')
    }
   
     df.result <-cbind(stable.samples,enterotypes,enterotypes.short)
     return(df.result)
}