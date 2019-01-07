# The programs tests if the enterotypes are correlated with Principal components
# ALL PORGRAMS HERE
# Input: 
# a file1 with col1 - number (automatically read as rowname), col2 - sample ID, col3 - enterotype_Driver_full, col4 -enterotype number - short
# a file2 with Principal component decomposition of the gene expression data. Column 1 and 2 = sample ID
# (matching to file1), columns 3+ = PCs. 1-200 respectively
# n = number of PCs to check
# locat = enterotypes of which location need to gather?
# Output: 
# see in each program description


# corrEnterExpr.initiation function
# Output: 
# [[1]] a matrix with info of normality and variance equality for each PC.
# matrix[,3]  - 1 if t-test can be used and 0 otherwise
# after running this program the result is that we can not use parametric tests for most of the PCs
# [[2]] the file with PCs and enterotypes added
corrEnterExpr.initiation <- function(file1='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype/Stable_enterotypes.csv',file2,locat,n=10){
  library('data.table')
  library('stringr')
  library(stats)
  enterotypes <- read.table(file1, header = T)
  colnames(enterotypes) <- c('stable.samples','enterotypes','enterotypes.short')
  PCs <- read.table(file2,header = T)

 # check that dimensions of files are nx3 and mxm+2
  if(dim(enterotypes)[2]!=3){
    stop()
    geterrmessage('number of columns in file1 should be 3: sampleID, enterotype_long, enterotype_short')
  }
  if((dim(PCs)[2]-dim(PCs)[1])!=2){
    stop()
    geterrmessage('number of columns in PCs table (file2) should be number of rows + 2. First 2 columns are sample name')
  }
  
  # check location and filter table: get only lines containing the location
  if (locat=='IL'){
    enterotypes <- enterotypes[enterotypes[,1] %like% '\\.IL',]
    enterotypes[,1] <- gsub('\\.IL','',enterotypes[,1])
  }else if (locat=='TR'){
    enterotypes <- enterotypes[enterotypes[,1] %like% '\\.TR',]
    enterotypes[,1] <- gsub('\\.TR','',enterotypes[,1])
  }else if (locat=='RE'){
    enterotypes <- enterotypes[enterotypes[,1] %like% '\\.RE',]
    enterotypes[,1] <- gsub('\\.RE','',enterotypes[,1])
  }
  
  # mod.PCs is a data table of PCs with enterotypes written to the first column instead 
  # of FID. NEED TO add eneterotype of 'sample.1' if 'sample' is absent!
  mod.PCs <- PCs
  match1 <- match(PCs$IID, enterotypes$stable.samples)
  mod.PCs[,1] <- enterotypes[match1,3]
  colnames(mod.PCs[1]) <- 'enterotype'
  # # Than, get a match for the 'vector of' rows with NAs in column 1 from mod.PCs both in 
  # # enterotypes and in the original mod.PCs table (as when we take only rows with NAs we 
  # # forget their numbers)   
  # matching.NA.enterotypes <- match(mod.PCs[is.na(mod.PCs[,1]),2],enterotypes.stable.samples.short)
  # matching.NA.mod.PCs <- match(mod.PCs[is.na(mod.PCs[,1]),2],mod.PCs[,2])
  # # Now, for each element of the 'vector of' rows with NAs that has correspondence to enterotypes' 
  # # table, get corresponding numbers of rows in mod.PCs table
  # match1 <- cbind(matching.NA.enterotypes,matching.NA.mod.PCs)
  # # leave only rows with some match to enterotype
  # match1 <- match1[!is.na(matching.NA.enterotypes),]
  # # finally, fill in the enterotypes
  # mod.PCs[match1[,2],1] <- enterotypes$enterotypes.short[match1[,1]]
  # 
  # # checked 5 random elemets of the mod.PCs table - enterotypes are correct!
  # 
  
  # Now check if residuals in the model are normally distributed
  normaliok <- matrix(0,n,3)
  for (i in 3:(n+2)){
  # PC #1 is in row 3. So, PC number = i-2
  # get p-values for shapiro test for normality of PC i-2 for each enterotype
    shapiro1 <- shapiro.test(mod.PCs[mod.PCs[,1]=='Bacteroides',i])$p.value
    shapiro2 <- shapiro.test(mod.PCs[mod.PCs[,1]=='Prevotella',i])$p.value
    # get p-values for F Test to Compare Variances for PC i-2 each enterotype pair
    var12 <- var.test(mod.PCs[mod.PCs[,1]=='Bacteroides',i],mod.PCs[mod.PCs[,1]=='Prevotella',i])$p.value
    normaliok[i-2,1] <- min(shapiro1,shapiro2)
    normaliok[i-2,2] <- var12
    normaliok[i-2,3] <- if (min(normaliok[i-2,1],normaliok[i-2,2])>=0.01){1}else{0}
  }
  normaliok <- round(normaliok,5)
  rownames(normaliok) <- c(3:(n+2))
  colnames(normaliok) <- c('p.value.normality','p.value.st_dev_equal','can.use.t.test')
  result <- list(mod.PCs,normaliok)
  names(result) <- c('PCs_Julia_with_enterotypes','can_we_use_t_test')
  return(result)
}

# Test for association of enterotype with PCs. 
# Originally Mann Whitney only, now added - 't.test'
# PC.number = the number of principal component for which we make the test
corrEnterExpr.test <- function(file_1, file_2,locat_,PC.number=1,test.type = 't.test'){
  i <- PC.number + 2
  test.p.val <- numeric(1)
  sample.size <- integer(2)
  average <- numeric(2)
  st.dev <- numeric(2)
  PCs <- corrEnterExpr.initiation(file1=file_1, file2 = file_2, locat = locat_)$PCs_Julia_with_enterotypes
  # depending on the test type, perform the statistic test
  if(test.type == 'wilcox') {
  test.p.val[1] <- wilcox.test(PCs[PCs[,1]=='Bacteroides',i],PCs[PCs[,1]=='Prevotella',i])$p.value
  # test.p.val[2] <- wilcox.test(PCs[PCs[,1]=='Bacteroides',i],PCs[PCs[,1]=='enterotype_3',i])$p.value
  # test.p.val[3] <- wilcox.test(PCs[PCs[,1]=='Prevotella',i],PCs[PCs[,1]=='enterotype_3',i])$p.value
  } else if (test.type == 't.test'){
    test.p.val[1] <- t.test(PCs[PCs[,1]=='Bacteroides',i],PCs[PCs[,1]=='Prevotella',i])$p.value
    # test.p.val[2] <- t.test(PCs[PCs[,1]=='Bacteroides',i],PCs[PCs[,1]=='enterotype_3',i])$p.value
    # test.p.val[3] <- t.test(PCs[PCs[,1]=='Prevotella',i],PCs[PCs[,1]=='enterotype_3',i])$p.value
  } else{
    stop('test.type must be wilcox or t.test')
    geterrmessage()
  }

  sample.size[1] <- sum(!is.na(PCs[PCs[,1]=='Bacteroides',i]))
  sample.size[2] <- sum(!is.na(PCs[PCs[,1]=='Prevotella',i]))
  # sample.size[3] <- sum(!is.na(PCs[PCs[,1]=='enterotype_3',i]))
  average[1] <- mean(PCs[PCs[,1]=='Bacteroides',i],na.rm = T)
  average[2] <- mean(PCs[PCs[,1]=='Prevotella',i],na.rm = T)
  # average[3] <- mean(PCs[PCs[,1]=='enterotype_3',i],na.rm = T)
  st.dev[1] <- sd(PCs[PCs[,1]=='Bacteroides',i],na.rm = T)
  st.dev[2] <- sd(PCs[PCs[,1]=='Prevotella',i],na.rm = T)
  # st.dev[3] <- sd(PCs[PCs[,1]=='enterotype_3',i],na.rm = T)
  # names(test.p.val) <- c('Bact-Prev','Bact-ent3','Prev-ent3')
  additional.info <- cbind(sample.size,average,st.dev)
  rownames(additional.info) <- c('Bacteroides','Prevotella')
  result <- list(test.p.val,additional.info)
  if(test.type == 'wilcox'){
    names(result) <- c('Mann.Whitney.pvalues','additional.info')
  } else if (test.type == 't.test'){
    names(result) <- c('T.test.pvalues','additional.info')
  }
  return(result)
}


# perform corrEnterExpr.test for first n PCs
corrEnterExpr.test.n <- function(file_1=file_1, file_2=file_2,locat_=locat_,n=3,test.type = 't.test'){
  associat <- list()
  for (i in 1:n){
    associat[[i]] <- corrEnterExpr.test(file_1=file_1,file_2=file_2,locat_ = locat_,PC.number=i,test.type = test.type)
  }
  if(test.type == 'wilcox'){
    associat.pvalues <- as.data.frame(associat[[1]]$Mann.Whitney.pvalues)
  }else if (test.type == 't.test'){
    associat.pvalues <- as.data.frame(associat[[1]]$T.test.pvalues)
  }else{
    stop('test.type must be wilcox or t.test')
    geterrmessage()
  }
  if(n>=2){
    for (i in 2:n){
    associat.pvalues <- cbind(associat.pvalues,as.data.frame(associat[[i]][[1]]))
    }
  }
  
  numbers <- seq(1,n)
  names(associat.pvalues)  <- paste0('PC',numbers,'.',test.type,'.pval')
  return(list(p.values=associat.pvalues, more.info=associat))
}