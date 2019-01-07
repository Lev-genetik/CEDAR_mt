# Input: a data frame sample-stable enterotype
# Output: how stable the enterotype is within locations
# IPC001: TRUE (if stable within 3 locations)
# IPC002: FALSE (if non-stable within 3 locations) (only for samples where we have 3 locations)
# than the same for samples that have only 2 locations

# works for data frames that contain at least 1 sample info about 3 locations
getStabStats <- function(df1){
  library('stringr')
  library('data.table')
  if((!is.data.frame(df1))|(dim(df1)[2]!=2)){
    stop()
    geterrmessage('The argument must be a data frame of 2 columns: sample-eterotype')
  }
  # print('2 stage')
  sample.numbers <- str_sub(df1[,1],start = 4,end=-4)
  first.occurence <- !duplicated(str_sub(df1[,1],start = 4,end=-4))
  samples.unique <- sample.numbers[first.occurence]
  # result.three.loc <- integer(27)
  # # generate names for the enterotype combinations
  # letters <- c('B','B','B','P','P','P','O','O','O')
  # combs.tmp <- combn(letters,3)
  # combs.vector <- character(dim(combs.tmp)[2])
  # for (i in 1:dim(combs.tmp)[2]){
  #   combs.vector[i] <- str_c(combs.tmp[,i])
  # }
  
# get enterotypes for all samples that have info for 3 locations: e.g. row 1 = sample 1, colons are locations IL|TR|RE
# get enterotypes for all samples that have info for 2 locations: e.g. row 1 = sample 1, colons are locations IL|TR|RE
  result.3 <- matrix('',ncol = 3,nrow = 0)
  result.2 <- matrix('',ncol = 2,nrow = 0)
  count.3=1
  count.2=1
  for (i in 1:length(samples.unique)){
    if (!is.null(df1[df1[,1] %like% samples.unique[i],])&(!is.vector(df1[df1[,1] %like% samples.unique[i],]))){
    if(dim(df1[df1[,1] %like% samples.unique[i],])[1]==3){
     m <- as.vector(df1[df1[,1] %like% samples.unique[i],2])
     result.3 <- rbind(result.3,m)
     rownames(result.3)[count.3]=samples.unique[i]
     count.3 <- count.3 + 1
    }
      else if(dim(df1[df1[,1] %like% samples.unique[i],])[1]==2){
     m <- as.vector(df1[df1[,1] %like% samples.unique[i],2])
     result.2 <- rbind(result.2,m)
     rownames(result.2)[count.2]=samples.unique[i]
     count.2 <- count.2 + 1

    }
    }
  }
  # generate vectors containing enterotypes for all 3/2 locations
  conc.3 <- character(dim(result.3)[1])
  conc.2 <- character(dim(result.2)[1])
  for (i in 1:dim(result.3)[1]){
    conc.3[i] <- paste(result.3[i,],collapse = '-')
  }
  for (i in 1:dim(result.2)[1]){
    conc.2[i] <- paste(result.2[i,],collapse = '-')
  }
  # Result: how many Bact-Bact-Bact etc.
  three.loc.table <- table(conc.3)
  two.loc.table <- table(conc.2)
  result.list <- list(three.loc.table,two.loc.table,result.3)
  names(result.list) <- c('three.loc.table','two.loc.table','three.loc.enterotypes')
  return(result.list)
  }

# simulation of random number of stable enterotypes across all locations
# N is the total number of samples that have info about enterotypes at all locations
# there.loc.enterotypes is output of getStabStats, list[[3]]. 


getNoiseStabStats <- function(df1,there.loc.enterotypes,N=30,shuffle=100){
  
  library('stringr')
  library('data.table')
  if((!is.data.frame(df1))|(dim(df1)[2]!=2)){
    stop()
    geterrmessage('The argument must be a data frame of 2 columns: sample-eterotype')
  }
  
  # get all samples from IL, TR and RE separately
  samples.IL <- df1[df1[,1] %like% '.IL',]
  samples.TR <- df1[df1[,1] %like% '.TR',]
  samples.RE <- df1[df1[,1] %like% '.RE',]
  
  # now calculate proportion of all enterotypes in each location
  Bact.IL <- (dim(samples.IL[samples.IL[,2] %like% 'Bacteroides',])[1])/(dim(samples.IL)[1])
  Prev.IL <- (dim(samples.IL[samples.IL[,2] %like% 'Prevotella',])[1])/(dim(samples.IL)[1])
  ent3.IL <- (dim(samples.IL[samples.IL[,2] %like% 'enterotype_3',])[1])/(dim(samples.IL)[1])
  Bact.TR <- (dim(samples.TR[samples.TR[,2] %like% 'Bacteroides',])[1])/(dim(samples.TR)[1])
  Prev.TR <- (dim(samples.TR[samples.TR[,2] %like% 'Prevotella',])[1])/(dim(samples.TR)[1])
  ent3.TR <- (dim(samples.TR[samples.TR[,2] %like% 'enterotype_3',])[1])/(dim(samples.TR)[1])
  Bact.RE <- (dim(samples.RE[samples.RE[,2] %like% 'Bacteroides',])[1])/(dim(samples.RE)[1])
  Prev.RE <- (dim(samples.RE[samples.RE[,2] %like% 'Prevotella',])[1])/(dim(samples.RE)[1])
  ent3.RE <- (dim(samples.RE[samples.RE[,2] %like% 'enterotype_3',])[1])/(dim(samples.RE)[1])
  
  # Now calculate expectance of each each enterotype=location combination
  BBB <- Bact.IL*Bact.RE*Bact.TR*N
  BBP <- Bact.IL*Bact.RE*Prev.TR*N
  BB3 <- Bact.IL*Bact.RE*ent3.TR*N
  B33 <- Bact.IL*ent3.RE*ent3.TR*N
  PPP <- Prev.IL*Prev.TR*Prev.RE*N
  e333 <- ent3.IL*ent3.TR*ent3.RE*N
  
  # Now get permutations of enterotypes and calculate stats - write to an array
  BBB.simul <- integer(shuffle)
  PPP.simul <- integer(shuffle)
  e333.simul <- integer(shuffle)
  BB3.simul <- integer(shuffle)
  B33.simul <- integer(shuffle)
  
  conc <- character(dim(there.loc.enterotypes)[1])
  there.loc.enterotypes.new <- there.loc.enterotypes
  # making the permutations
  for (i in 1:shuffle){
    # making the permutations per se
    there.loc.enterotypes.new[,1] <- sample(there.loc.enterotypes[,1])
    there.loc.enterotypes.new[,2] <- sample(there.loc.enterotypes[,2])
    there.loc.enterotypes.new[,3] <- sample(there.loc.enterotypes[,3])
    
    for (j in 1:dim(there.loc.enterotypes.new)[1]){
      conc[j] <- paste(there.loc.enterotypes.new[j,],collapse = '-')
    }
    table.current <- table(conc)
    
    # recording the values of each combination in this permutation
    BBB.simul[i] <- max(as.integer((as.list(table.current))$'Bacteroides-Bacteroides-Bacteroides'),0)
    PPP.simul[i] <- max(as.integer((as.list(table.current))$'Prevotella-Prevotella-Prevotella'),0)
    e333.simul[i] <- max(as.integer((as.list(table.current))$'enterotype_3-enterotype_3-enterotype_3'),0)
    BB3.simul[i] <- max(as.integer((as.list(table.current))$'Bacteroides-Bacteroides-enterotype_3'),0)
    B33.simul[i] <- max(as.integer((as.list(table.current))$'Bacteroides-enterotype_3-enterotype_3'),0)

    
  }
  BBB.stats <- table(BBB.simul)
  PPP.stats <- table(PPP.simul)
  e333.stats <- table(e333.simul)
  BB3.stats <- table(BB3.simul)
  B33.stats <- table(B33.simul)
  result <- list(BBB.stats,PPP.stats,e333.stats,BB3.stats,B33.stats)
  names(result) <- c('BBB.stats','PPP.stats','e333.stats','BB3.stats','B33.stats')
  return(result)
} 

getCondProb.of.enterotypes <- function(df1){
  library('stringr')
  library('data.table')
  if((!is.data.frame(df1))|(dim(df1)[2]!=2)){
    stop()
    geterrmessage('The argument must be a data frame of 2 columns: sample-eterotype')
  }
  
  # get all samples from IL, TR and RE separately
  samples.IL <- df1[df1[,1] %like% '.IL',]
  samples.TR <- df1[df1[,1] %like% '.TR',]
  samples.RE <- df1[df1[,1] %like% '.RE',]
  
  # get rid of everything except sample numbers
  samples.IL[,1] <- str_sub(samples.IL[,1],start = 4,end=-4)
  samples.TR[,1] <- str_sub(samples.TR[,1],start = 4,end=-4)
  samples.RE[,1] <- str_sub(samples.RE[,1],start = 4,end=-4)

 get_cond_prob <- function(ent_cons,loc_cons,ent_cause,loc_cause){
    # comments are written as if the task was to evaluate conditional probability of Bacteroides in TR in case of Bacteroides in IL
    #get all Bacteroides in IL 
    if(loc_cause=='IL'){
      samples.loc_cause.ent_cause <- samples.IL[samples.IL[,2]==ent_cause,]
    }else if (loc_cause=='TR'){
      samples.loc_cause.ent_cause <- samples.TR[samples.TR[,2]==ent_cause,]
    }else if (loc_cause=='RE'){
      samples.loc_cause.ent_cause <- samples.RE[samples.RE[,2]==ent_cause,]  
    }

    if(loc_cons=='IL'){
      samples.loc_cons <- samples.IL
    } else if (loc_cons=='TR'){
      samples.loc_cons <- samples.TR
    } else if (loc_cons=='RE'){
      samples.loc_cons <- samples.RE
    }
    # intersect samples.IL.Bacteroides with samples.TR
    common_samples.loc_cause.ent_cause.samples.loc_cons <- intersect(samples.loc_cause.ent_cause[,1],samples.loc_cons[,1])
    # now get number of Bacteroides in TR out of samples with Bacteroides in IL
    bact=0
    prev=0
    ent3=0
    for (i in 1:length(common_samples.loc_cause.ent_cause.samples.loc_cons)){
      if (samples.loc_cons[samples.loc_cons[,1]==common_samples.loc_cause.ent_cause.samples.loc_cons[i],2]=='Bacteroides'){
        bact = bact + 1
      } else if (samples.loc_cons[samples.loc_cons[,1]==common_samples.loc_cause.ent_cause.samples.loc_cons[i],2]=='Prevotella'){
        prev = prev + 1
      } else if (samples.loc_cons[samples.loc_cons[,1]==common_samples.loc_cause.ent_cause.samples.loc_cons[i],2]=='enterotype_3'){
        ent3 = ent3 + 1
      }
    }
    sum1 = bact+prev+ent3
    if (length(common_samples.loc_cause.ent_cause.samples.loc_cons)!= sum1){
      stop()
      geterrmessage('Enterotypes should be Bacteroides, Prevotella and enterotype_3 only!')
    }
    # Now get percentage of ent_consequence in loc_consequence among the elements with ent_cause in loc_cause
    if(ent_cons=='Bacteroides'){
      p_ent_cons_ent_cause <- bact/sum1
    } else if(ent_cons=='Prevotella'){
      p_ent_cons_ent_cause <- prev/sum1
    } else if(ent_cons=='enterotype_3'){
      p_ent_cons_ent_cause <- ent3/sum1
    }
    rezultat <- list(p_ent_cons_ent_cause,length(common_samples.loc_cause.ent_cause.samples.loc_cons))
    names(rezultat) <- c('enterotype_conditional_probability','number_samples_for_evaluation_of_it')
    # Now get probability of ent_consequence in loc_consequence
    # p_ent_cons.loc_cons <- length(samples.loc_cons[samples.loc_cons[,2]==ent_cause,2])/(dim(samples.loc_cons)[1])
    return(rezultat)      
  }
 result.B1 <- matrix(0,3,3)
 rownames(result.B1) <- c('Bact_IL','Bact_TR','Bact_Re')
 colnames(result.B1) <- c('Bact_IL','Bact_TR','Bact_Re')
 result.B1[1,1] <- round(get_cond_prob(ent_cons = 'Bacteroides',loc_cons = 'IL',ent_cause='Bacteroides',loc_cause='IL')[[1]],3)
 result.B1[2,1] <- round(get_cond_prob(ent_cons = 'Bacteroides',loc_cons = 'TR',ent_cause='Bacteroides',loc_cause='IL')[[1]],3)
 result.B1[3,1] <- round(get_cond_prob(ent_cons = 'Bacteroides',loc_cons = 'RE',ent_cause='Bacteroides',loc_cause='IL')[[1]],3)
 
 result.B1[1,2] <- round(get_cond_prob(ent_cons = 'Bacteroides',loc_cons = 'IL',ent_cause='Bacteroides',loc_cause='TR')[[1]],3)
 result.B1[2,2] <- round(get_cond_prob(ent_cons = 'Bacteroides',loc_cons = 'TR',ent_cause='Bacteroides',loc_cause='TR')[[1]],3)
 result.B1[3,2] <- round(get_cond_prob(ent_cons = 'Bacteroides',loc_cons = 'RE',ent_cause='Bacteroides',loc_cause='TR')[[1]],3)
 
 result.B1[1,3] <- round(get_cond_prob(ent_cons = 'Bacteroides',loc_cons = 'IL',ent_cause='Bacteroides',loc_cause='RE')[[1]],3)
 result.B1[2,3] <- round(get_cond_prob(ent_cons = 'Bacteroides',loc_cons = 'TR',ent_cause='Bacteroides',loc_cause='RE')[[1]],3)
 result.B1[3,3] <- round(get_cond_prob(ent_cons = 'Bacteroides',loc_cons = 'RE',ent_cause='Bacteroides',loc_cause='RE')[[1]],3)
 
 result.B2 <- matrix(0,3,3)
 rownames(result.B2) <- c('Bact_IL','Bact_TR','Bact_Re')
 colnames(result.B2) <- c('Bact_IL','Bact_TR','Bact_Re')
 result.B2[1,1] <- round(get_cond_prob(ent_cons = 'Bacteroides',loc_cons = 'IL',ent_cause='Bacteroides',loc_cause='IL')[[2]],3)
 result.B2[2,1] <- round(get_cond_prob(ent_cons = 'Bacteroides',loc_cons = 'TR',ent_cause='Bacteroides',loc_cause='IL')[[2]],3)
 result.B2[3,1] <- round(get_cond_prob(ent_cons = 'Bacteroides',loc_cons = 'RE',ent_cause='Bacteroides',loc_cause='IL')[[2]],3)
 
 result.B2[1,2] <- round(get_cond_prob(ent_cons = 'Bacteroides',loc_cons = 'IL',ent_cause='Bacteroides',loc_cause='TR')[[2]],3)
 result.B2[2,2] <- round(get_cond_prob(ent_cons = 'Bacteroides',loc_cons = 'TR',ent_cause='Bacteroides',loc_cause='TR')[[2]],3)
 result.B2[3,2] <- round(get_cond_prob(ent_cons = 'Bacteroides',loc_cons = 'RE',ent_cause='Bacteroides',loc_cause='TR')[[2]],3)
 
 result.B2[1,3] <- round(get_cond_prob(ent_cons = 'Bacteroides',loc_cons = 'IL',ent_cause='Bacteroides',loc_cause='RE')[[2]],3)
 result.B2[2,3] <- round(get_cond_prob(ent_cons = 'Bacteroides',loc_cons = 'TR',ent_cause='Bacteroides',loc_cause='RE')[[2]],3)
 result.B2[3,3] <- round(get_cond_prob(ent_cons = 'Bacteroides',loc_cons = 'RE',ent_cause='Bacteroides',loc_cause='RE')[[2]],3)

 # Now get result for Prevotella - the same
 
 result.P1 <- matrix(0,3,3)
 rownames(result.P1) <- c('Prev_IL','Prev_TR','Prev_Re')
 colnames(result.P1) <- c('Prev_IL','Prev_TR','Prev_Re')
 result.P1[1,1] <- round(get_cond_prob(ent_cons = 'Prevotella',loc_cons = 'IL',ent_cause='Prevotella',loc_cause='IL')[[1]],3)
 result.P1[2,1] <- round(get_cond_prob(ent_cons = 'Prevotella',loc_cons = 'TR',ent_cause='Prevotella',loc_cause='IL')[[1]],3)
 result.P1[3,1] <- round(get_cond_prob(ent_cons = 'Prevotella',loc_cons = 'RE',ent_cause='Prevotella',loc_cause='IL')[[1]],3)
 
 result.P1[1,2] <- round(get_cond_prob(ent_cons = 'Prevotella',loc_cons = 'IL',ent_cause='Prevotella',loc_cause='TR')[[1]],3)
 result.P1[2,2] <- round(get_cond_prob(ent_cons = 'Prevotella',loc_cons = 'TR',ent_cause='Prevotella',loc_cause='TR')[[1]],3)
 result.P1[3,2] <- round(get_cond_prob(ent_cons = 'Prevotella',loc_cons = 'RE',ent_cause='Prevotella',loc_cause='TR')[[1]],3)
 
 result.P1[1,3] <- round(get_cond_prob(ent_cons = 'Prevotella',loc_cons = 'IL',ent_cause='Prevotella',loc_cause='RE')[[1]],3)
 result.P1[2,3] <- round(get_cond_prob(ent_cons = 'Prevotella',loc_cons = 'TR',ent_cause='Prevotella',loc_cause='RE')[[1]],3)
 result.P1[3,3] <- round(get_cond_prob(ent_cons = 'Prevotella',loc_cons = 'RE',ent_cause='Prevotella',loc_cause='RE')[[1]],3)
 
 result.P2 <- matrix(0,3,3)
 rownames(result.P2) <- c('Prev_IL','Prev_TR','Prev_Re')
 colnames(result.P2) <- c('Prev_IL','Prev_TR','Prev_Re')
 result.P2[1,1] <- round(get_cond_prob(ent_cons = 'Prevotella',loc_cons = 'IL',ent_cause='Prevotella',loc_cause='IL')[[2]],3)
 result.P2[2,1] <- round(get_cond_prob(ent_cons = 'Prevotella',loc_cons = 'TR',ent_cause='Prevotella',loc_cause='IL')[[2]],3)
 result.P2[3,1] <- round(get_cond_prob(ent_cons = 'Prevotella',loc_cons = 'RE',ent_cause='Prevotella',loc_cause='IL')[[2]],3)
 
 result.P2[1,2] <- round(get_cond_prob(ent_cons = 'Prevotella',loc_cons = 'IL',ent_cause='Prevotella',loc_cause='TR')[[2]],3)
 result.P2[2,2] <- round(get_cond_prob(ent_cons = 'Prevotella',loc_cons = 'TR',ent_cause='Prevotella',loc_cause='TR')[[2]],3)
 result.P2[3,2] <- round(get_cond_prob(ent_cons = 'Prevotella',loc_cons = 'RE',ent_cause='Prevotella',loc_cause='TR')[[2]],3)
 
 result.P2[1,3] <- round(get_cond_prob(ent_cons = 'Prevotella',loc_cons = 'IL',ent_cause='Prevotella',loc_cause='RE')[[2]],3)
 result.P2[2,3] <- round(get_cond_prob(ent_cons = 'Prevotella',loc_cons = 'TR',ent_cause='Prevotella',loc_cause='RE')[[2]],3)
 result.P2[3,3] <- round(get_cond_prob(ent_cons = 'Prevotella',loc_cons = 'RE',ent_cause='Prevotella',loc_cause='RE')[[2]],3)
 
 result.E1 <- matrix(0,3,3)
 rownames(result.E1) <- c('Ent3_IL','Ent3_TR','Ent3_Re')
 colnames(result.E1) <- c('Ent3_IL','Ent3_TR','Ent3_Re')
 result.E1[1,1] <- round(get_cond_prob(ent_cons = 'enterotype_3',loc_cons = 'IL',ent_cause='enterotype_3',loc_cause='IL')[[1]],3)
 result.E1[2,1] <- round(get_cond_prob(ent_cons = 'enterotype_3',loc_cons = 'TR',ent_cause='enterotype_3',loc_cause='IL')[[1]],3)
 result.E1[3,1] <- round(get_cond_prob(ent_cons = 'enterotype_3',loc_cons = 'RE',ent_cause='enterotype_3',loc_cause='IL')[[1]],3)
 
 result.E1[1,2] <- round(get_cond_prob(ent_cons = 'enterotype_3',loc_cons = 'IL',ent_cause='enterotype_3',loc_cause='TR')[[1]],3)
 result.E1[2,2] <- round(get_cond_prob(ent_cons = 'enterotype_3',loc_cons = 'TR',ent_cause='enterotype_3',loc_cause='TR')[[1]],3)
 result.E1[3,2] <- round(get_cond_prob(ent_cons = 'enterotype_3',loc_cons = 'RE',ent_cause='enterotype_3',loc_cause='TR')[[1]],3)
 
 result.E1[1,3] <- round(get_cond_prob(ent_cons = 'enterotype_3',loc_cons = 'IL',ent_cause='enterotype_3',loc_cause='RE')[[1]],3)
 result.E1[2,3] <- round(get_cond_prob(ent_cons = 'enterotype_3',loc_cons = 'TR',ent_cause='enterotype_3',loc_cause='RE')[[1]],3)
 result.E1[3,3] <- round(get_cond_prob(ent_cons = 'enterotype_3',loc_cons = 'RE',ent_cause='enterotype_3',loc_cause='RE')[[1]],3)
 
 result.E2 <- matrix(0,3,3)
 rownames(result.E2) <- c('Bact_IL','Bact_TR','Bact_Re')
 colnames(result.E2) <- c('Bact_IL','Bact_TR','Bact_Re')
 result.E2[1,1] <- round(get_cond_prob(ent_cons = 'enterotype_3',loc_cons = 'IL',ent_cause='enterotype_3',loc_cause='IL')[[2]],3)
 result.E2[2,1] <- round(get_cond_prob(ent_cons = 'enterotype_3',loc_cons = 'TR',ent_cause='enterotype_3',loc_cause='IL')[[2]],3)
 result.E2[3,1] <- round(get_cond_prob(ent_cons = 'enterotype_3',loc_cons = 'RE',ent_cause='enterotype_3',loc_cause='IL')[[2]],3)
 
 result.E2[1,2] <- round(get_cond_prob(ent_cons = 'enterotype_3',loc_cons = 'IL',ent_cause='enterotype_3',loc_cause='TR')[[2]],3)
 result.E2[2,2] <- round(get_cond_prob(ent_cons = 'enterotype_3',loc_cons = 'TR',ent_cause='enterotype_3',loc_cause='TR')[[2]],3)
 result.E2[3,2] <- round(get_cond_prob(ent_cons = 'enterotype_3',loc_cons = 'RE',ent_cause='enterotype_3',loc_cause='TR')[[2]],3)
 
 result.E2[1,3] <- round(get_cond_prob(ent_cons = 'enterotype_3',loc_cons = 'IL',ent_cause='enterotype_3',loc_cause='RE')[[2]],3)
 result.E2[2,3] <- round(get_cond_prob(ent_cons = 'enterotype_3',loc_cons = 'TR',ent_cause='enterotype_3',loc_cause='RE')[[2]],3)
 result.E2[3,3] <- round(get_cond_prob(ent_cons = 'enterotype_3',loc_cons = 'RE',ent_cause='enterotype_3',loc_cause='RE')[[2]],3)
 
 consistency <- matrix(0,3,2)
 rownames(consistency) <- c('Bacteroides','Prevotella','Enterotype 3')
 colnames(consistency) <- c('Av.cons_by_loc','Av.cons_by_sample')
 consistency[1,1] <- (sum(result.B1)-3)/6
 consistency[2,1] <- (sum(result.P1)-3)/6
 consistency[3,1] <- (sum(result.E1)-3)/6
 consistency[1,2] <- (result.B1[2,1]*result.B2[2,1]+result.B1[3,1]*result.B2[3,1]+result.B1[1,2]*result.B2[1,2]+result.B1[2,3]*result.B2[2,3]+result.B1[1,3]*result.B2[1,3]+result.B1[3,2]*result.B2[3,2])/(result.B2[2,1]+result.B2[3,1]+result.B2[1,2]+result.B2[3,2]+result.B2[1,3]+result.B2[2,3])
 consistency[2,2] <- (result.P1[2,1]*result.P2[2,1]+result.P1[3,1]*result.P2[3,1]+result.P1[1,2]*result.P2[1,2]+result.P1[2,3]*result.P2[2,3]+result.P1[1,3]*result.P2[1,3]+result.P1[3,2]*result.P2[3,2])/(result.P2[2,1]+result.P2[3,1]+result.P2[1,2]+result.P2[3,2]+result.P2[1,3]+result.P2[2,3])
 consistency[3,2] <- (result.E1[2,1]*result.E2[2,1]+result.E1[3,1]*result.E2[3,1]+result.E1[1,2]*result.E2[1,2]+result.E1[2,3]*result.E2[2,3]+result.E1[1,3]*result.E2[1,3]+result.E1[3,2]*result.E2[3,2])/(result.E2[2,1]+result.E2[3,1]+result.E2[1,2]+result.E2[3,2]+result.E2[1,3]+result.E2[2,3])
 
 consistency <- round(consistency[,],4)
 
 result.list <- list(result.B1,result.B2,result.P1,result.P2,result.E1,result.E2,consistency)
 names(result.list) <- c('Bacteroides_cond_prob','Bacteroides_sample_#(e.g.[2,1]=total_samples_TR_that_match_to_Bacteroides_in_IL)','Prevotella_cond_prob','Prevotella_sample_#(e.g.[2,1]=total_samples_TR_that_match_to_Prevotella_in_IL)','Enterotype_3_cond_prob','Enterotype_3_sample_#(e.g.[2,1]=total_samples_TR_that_match_to_enterotype_3_in_IL)','Average enterotype stability diff. locations')
 return(result.list)
}
