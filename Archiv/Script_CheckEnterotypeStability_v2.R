#  Script to check enterotype stability
# result: cond.probabilities - conditional probabilities of presence of enterotypes as well as
# overall stability
# StabStats (table of samples that have amplicon-switch-stable enterotype across each of locations) and result (random simulation)
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/CheckEnterotypeStability_v2.R')
# get enterotypes stable over different 16S locations
stable.enterotypes <- read.table('/media/lev-genetik/980E73270E72FD96/Liege/Results/Clustering/Sample_cluster/2018_09_06_enterotyping/Stable_enterotypes.csv')
# sample-enterotype table
stable.enterotypes <- subset(stable.enterotypes, select=c(1,3))
# get a table of samples that have amplicon-switch-stable enterotype across each of locations 
StabStats <- getStabStats(stable.enterotypes)
# random simulation
result <- getNoiseStabStats(stable.enterotypes,StabStats$three.loc.enterotypes,shuffle = 10000)

# get conditional probabilities of presence of enterotypes, e.g. Bacteroides enterotype probability in TR for the sample where the enterotype for IL is Bacteroides
cond.probabilities <- getCondProb.of.enterotypes(stable.enterotypes)
