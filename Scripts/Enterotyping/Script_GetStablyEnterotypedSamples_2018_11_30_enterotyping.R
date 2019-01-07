# Get stable and consensus enterotypes across amplicons for IL, TR and RE separately
source(file.path(progr.dir,'GetStablyEnterotypedSamples_across_amplicons_smart_v4.R'))
ent.stability.stats <- getStableEntAcrossLoc()

# Get consensus enterotypes across locations given consensus ent for each location, including:
# 1. concordant consensus IL + consensus TR + consensus RE
# 2. no consensus enterotype for 1 of the locations but concordant for other 2 locations
# 3. consensus enterotype present for 1 location only
absolutely.stable.enterotypes.2e <- getStabAcrAmpl.final.2ent(file1=file.path(home.dir,'Results/Enterotyping/2018_11_30_enterotyping/Sample-enterotype/Consensus_enterotypes_IL.csv'),
                          file2=file.path(home.dir,'Results/Enterotyping/2018_11_30_enterotyping/Sample-enterotype/Consensus_enterotypes_TR.csv'),
                          file3=file.path(home.dir,'Results/Enterotyping/2018_11_30_enterotyping/Sample-enterotype/Consensus_enterotypes_RE.csv'),
                          stable.three=T, stable.two=T, stable.one=T, unstable.three = F, unstable.two = F, 
                          correspondence=NA,type='across.locations')
# in column 2, we need only k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Bacteroidaceae.g__Bacteroides
# and k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Prevotellaceae.g__Prevotella
absolutely.stable.enterotypes.2e[,2] <- ''
absolutely.stable.enterotypes.2e[absolutely.stable.enterotypes.2e[,3]=='Bacteroides',2] <- 
  'k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Bacteroidaceae.g__Bacteroides'
absolutely.stable.enterotypes.2e[absolutely.stable.enterotypes.2e[,3]=='Prevotella',2] <- 
  'k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Prevotellaceae.g__Prevotella'
# check that for the replicates we have same results
replicates <- absolutely.stable.enterotypes.2e[absolutely.stable.enterotypes.2e[,1]%like%'\\.1',]
replicates[,1] <- gsub('\\.[0-9]','',replicates[,1])
samples.unique.in.replicates <- character(0)
replica.original.concordance = 0
replica.original.discordance = 0
for (i in 1:dim(replicates)[1]){
  if(is.na(match(replicates[,1],absolutely.stable.enterotypes.2e[,1])[i])){
    samples.unique.in.replicates <- c(samples.unique.in.replicates,replicates[i,])
  } else if(absolutely.stable.enterotypes.2e[match(replicates[,1],absolutely.stable.enterotypes.2e[,1])[i],3]==replicates[i,3]){
    replica.original.concordance = replica.original.concordance+1
  } else if(absolutely.stable.enterotypes.2e[match(replicates[,1],absolutely.stable.enterotypes.2e[,1])[i],3]!=replicates[i,3]){
    replica.original.discordance = replica.original.discordance+1
  }
}
print(paste0('Fraction of cross-amplicon-and-location-consensus enterotypes consistent across biological replicates: ',
             signif(replica.original.concordance/(replica.original.concordance+replica.original.discordance),digits = 3),
             ', N = ',dim(replicates)[1]))
# delete all replicates from the absolutely.stable.enterotypes.2e table
absolutely.stable.enterotypes.2e <- absolutely.stable.enterotypes.2e[!absolutely.stable.enterotypes.2e[,1]%like%'\\.1',]
# add the replicates that do not correspond to any original sample to the absolutely.stable.enterotypes.2e table
absolutely.stable.enterotypes.2e <- rbind(absolutely.stable.enterotypes.2e,replicates[replicates[,1]%in%samples.unique.in.replicates,])
# and finally sort the table 
absolutely.stable.enterotypes.2e <- absolutely.stable.enterotypes.2e[order(as.numeric(str_sub(absolutely.stable.enterotypes.2e[,1],start = 4))),]
# absolutely.stable.enterotypes.2e <- absolutely.stable.enterotypes.2e[,c(1,3)]
rownames(absolutely.stable.enterotypes.2e) <- NULL
colnames(absolutely.stable.enterotypes.2e) <- c('sample','enterotype','enterotype.short')
setwd(file.path(home.dir,'Results/Enterotyping/2018_11_30_enterotyping/Sample-enterotype'))
write.table(absolutely.stable.enterotypes.2e,'CEDAR_consensus_enterotypes_across_ampl_and_loc.csv')
