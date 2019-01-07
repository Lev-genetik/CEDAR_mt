# Script to check enterotype stability
# only V12 and V56 amplicons taken into account in the input file  of stable enterotypes across amplicons
# result: cond.probabilities - conditional probabilities of presence of enterotypes as well as
# overall stability
# StabStats (table of samples that have amplicon-switch-stable enterotype across each of locations) and result (random simulation)
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/CheckEnterotypeStability_v3.R')
# get enterotypes stable over different 16S locations
stable.enterotypes <- read.table('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype/Stable_enterotypes_v12_v56.csv')
# sample-enterotype table
stable.enterotypes <- subset(stable.enterotypes, select=c(1,3))
# get a table of samples that have amplicon-switch-stable enterotype across each of locations 
StabStats <- getStabStats(stable.enterotypes)
# random simulation
result <- getNoiseStabStats(stable.enterotypes,StabStats$three.loc.enterotypes,shuffle = 10000)

# get conditional probabilities of presence of enterotypes, e.g. Bacteroides enterotype probability in TR for the sample where the enterotype for IL is Bacteroides
cond.probabilities <- getCondProb.of.enterotypes(stable.enterotypes)
