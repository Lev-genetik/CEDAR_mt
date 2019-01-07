# This script gets absolutely stable enterotypes: across both locations (v12+v56) and amplicons

# this is part of Script_CheckEnterotypeStability_v3.R 

source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/CheckEnterotypeStability_v3.R')
# get enterotypes stable over different 16S locations
stable.enterotypes <- read.table('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype/Stable_enterotypes_v12_v56.csv')
# sample-enterotype table
stable.enterotypes <- subset(stable.enterotypes, select=c(1,3))
# get a table of samples that have amplicon-switch-stable enterotype across each of locations 
StabStats <- getStabStats(stable.enterotypes)

# now write stably enterotyped samples to a file
# get samples that have stable enterotype for IL=TR=RE
three.loc.ent <- StabStats$three.loc.enterotypes
three.loc.ent <- three.loc.ent[three.loc.ent[,1]==three.loc.ent[,2],]
three.loc.ent <- three.loc.ent[three.loc.ent[,1]==three.loc.ent[,3],]
three.loc.stable.ent <- three.loc.ent
# we know that no samples.1 are stable across all locations
# so, get sample-enterotype table for samples with absolutely stable enterotype
abs.stable.ent <- stable.enterotypes.TR[match(rownames(three.loc.stable.ent),str_sub(stable.enterotypes.TR[,1],start = 4,end=-4)),]
# rename rows
rownames(abs.stable.ent) <- c(1:dim(abs.stable.ent)[1])
# add third row and remove .TR, add .all for the samples
abs.stable.ent[,3] <- abs.stable.ent[,2]
colnames(abs.stable.ent)[3] <- colnames(abs.stable.ent)[2]
abs.stable.ent[,1] <- gsub('\\.TR','\\.all',abs.stable.ent[,1])
write.table(abs.stable.ent,file = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype/Absolutely_stable_enterotypes_v12_v56_IL_TR_RE.csv')

