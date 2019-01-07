# Script to determine if IL enetrotypes are less stable than TR and RE (pairwise comparison for 
# different amplicons)
# by default, enterotyping done separately by location, 
# Input files can be found in /media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype
# Else, if want to check all locations enterotyped together, can use work_dir = 
# '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype/Split_by_locat'
# N - the number of simulations
whichLocMostStabl <- function(N=10000,work_dir='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/'){
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Enterotyping_compar_functions_v5.R')

# Step 1. Prepare what we need for calculations
# prepare matrix loc.concordance - number of samples similar 
# across amplicons (in row.names)
loc.concordance.seplocatent <- matrix(0,nrow=12,ncol=3) 
rownames(loc.concordance.seplocatent) <- c('v12-v34_2ent','v12-v56_2ent','v34-v56_2ent','v12-v34_3ent','v12-v56_3ent','v34-v56_3ent','v12-v34_4ent','v12-v56_4ent','v34-v56_4ent','v12-v34_5ent','v12-v56_5ent','v34-v56_5ent')
colnames(loc.concordance.seplocatent) <- c('IL','TR','RE')

loc.discordance.seplocatent <- matrix(0,nrow=12,ncol=3) 
rownames(loc.discordance.seplocatent) <- c('v12-v34_2ent','v12-v56_2ent','v34-v56_2ent','v12-v34_3ent','v12-v56_3ent','v34-v56_3ent','v12-v34_4ent','v12-v56_4ent','v34-v56_4ent','v12-v34_5ent','v12-v56_5ent','v34-v56_5ent')
colnames(loc.discordance.seplocatent) <- c('IL','TR','RE')

# prepare the table of concordance, all locations enterotyped separately
# first, read the sample-enterotype tables
setwd(work_dir)
# v1v2
ent.Raes12.2e.sep.by.loc.IL <- read.table('V1V2_2e_IL_sample-enterotype.csv', header = T)
ent.Raes12.2e.sep.by.loc.TR <- read.table('V1V2_2e_TR_sample-enterotype.csv', header = T)
ent.Raes12.2e.sep.by.loc.RE <- read.table('V1V2_2e_RE_sample-enterotype.csv', header = T)

ent.Raes12.3e.sep.by.loc.IL <- read.table('V1V2_3e_IL_sample-enterotype.csv', header = T)
ent.Raes12.3e.sep.by.loc.TR <- read.table('V1V2_3e_TR_sample-enterotype.csv', header = T)
ent.Raes12.3e.sep.by.loc.RE <- read.table('V1V2_3e_RE_sample-enterotype.csv', header = T)

ent.Raes12.4e.sep.by.loc.IL <- read.table('V1V2_4e_IL_sample-enterotype.csv', header = T)
ent.Raes12.4e.sep.by.loc.TR <- read.table('V1V2_4e_TR_sample-enterotype.csv', header = T)
ent.Raes12.4e.sep.by.loc.RE <- read.table('V1V2_4e_RE_sample-enterotype.csv', header = T)

ent.Raes12.5e.sep.by.loc.IL <- read.table('V1V2_5e_IL_sample-enterotype.csv', header = T)
ent.Raes12.5e.sep.by.loc.TR <- read.table('V1V2_5e_TR_sample-enterotype.csv', header = T)
ent.Raes12.5e.sep.by.loc.RE <- read.table('V1V2_5e_RE_sample-enterotype.csv', header = T)

# v3v4
ent.Raes34.2e.sep.by.loc.IL <- read.table('V3V4_2e_IL_sample-enterotype.csv', header = T)
ent.Raes34.2e.sep.by.loc.TR <- read.table('V3V4_2e_TR_sample-enterotype.csv', header = T)
ent.Raes34.2e.sep.by.loc.RE <- read.table('V3V4_2e_RE_sample-enterotype.csv', header = T)

ent.Raes34.3e.sep.by.loc.IL <- read.table('V3V4_3e_IL_sample-enterotype.csv', header = T)
ent.Raes34.3e.sep.by.loc.TR <- read.table('V3V4_3e_TR_sample-enterotype.csv', header = T)
ent.Raes34.3e.sep.by.loc.RE <- read.table('V3V4_3e_RE_sample-enterotype.csv', header = T)

ent.Raes34.4e.sep.by.loc.IL <- read.table('V3V4_4e_IL_sample-enterotype.csv', header = T)
ent.Raes34.4e.sep.by.loc.TR <- read.table('V3V4_4e_TR_sample-enterotype.csv', header = T)
ent.Raes34.4e.sep.by.loc.RE <- read.table('V3V4_4e_RE_sample-enterotype.csv', header = T)

ent.Raes34.5e.sep.by.loc.IL <- read.table('V3V4_5e_IL_sample-enterotype.csv', header = T)
ent.Raes34.5e.sep.by.loc.TR <- read.table('V3V4_5e_TR_sample-enterotype.csv', header = T)
ent.Raes34.5e.sep.by.loc.RE <- read.table('V3V4_5e_RE_sample-enterotype.csv', header = T)

# v5v6
ent.Raes56.2e.sep.by.loc.IL <- read.table('V5V6_2e_IL_sample-enterotype.csv', header = T)
ent.Raes56.2e.sep.by.loc.TR <- read.table('V5V6_2e_TR_sample-enterotype.csv', header = T)
ent.Raes56.2e.sep.by.loc.RE <- read.table('V5V6_2e_RE_sample-enterotype.csv', header = T)

ent.Raes56.3e.sep.by.loc.IL <- read.table('V5V6_3e_IL_sample-enterotype.csv', header = T)
ent.Raes56.3e.sep.by.loc.TR <- read.table('V5V6_3e_TR_sample-enterotype.csv', header = T)
ent.Raes56.3e.sep.by.loc.RE <- read.table('V5V6_3e_RE_sample-enterotype.csv', header = T)

ent.Raes56.4e.sep.by.loc.IL <- read.table('V5V6_4e_IL_sample-enterotype.csv', header = T)
ent.Raes56.4e.sep.by.loc.TR <- read.table('V5V6_4e_TR_sample-enterotype.csv', header = T)
ent.Raes56.4e.sep.by.loc.RE <- read.table('V5V6_4e_RE_sample-enterotype.csv', header = T)

ent.Raes56.5e.sep.by.loc.IL <- read.table('V5V6_5e_IL_sample-enterotype.csv', header = T)
ent.Raes56.5e.sep.by.loc.TR <- read.table('V5V6_5e_TR_sample-enterotype.csv', header = T)
ent.Raes56.5e.sep.by.loc.RE <- read.table('V5V6_5e_RE_sample-enterotype.csv', header = T)

# define the total number of samples in each location based on the tables read. We assume that number of samples for 2 enterotype groups
# =3 ent gr = 4 = 5
total.IL <- dim(ent.Raes12.2e.sep.by.loc.IL)[1]
total.TR <- dim(ent.Raes12.2e.sep.by.loc.TR)[1]
total.RE <- dim(ent.Raes12.2e.sep.by.loc.RE)[1]


# Step 2.
# Now calculate the number of samples the same [loc.concordance.seplocatent](and different [loc.discordance.seplocatent]) across locations
setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype')
# 1 row = 12vs34_2e comparison of clusterings for IL, TR and RE (the corresponding columns)
loc.concordance.seplocatent[1,1] <- as.numeric(CompEntPred.pair(ent.Raes12.2e.sep.by.loc.IL,ent.Raes34.2e.sep.by.loc.IL)$Maximum_number_samples_clustering_together[,1])
loc.discordance.seplocatent[1,1] <- total.IL - loc.concordance.seplocatent[1,1]
loc.concordance.seplocatent[1,2] <- as.numeric(CompEntPred.pair(ent.Raes12.2e.sep.by.loc.TR,ent.Raes34.2e.sep.by.loc.TR)$Maximum_number_samples_clustering_together[,1])
loc.discordance.seplocatent[1,2] <- total.TR - loc.concordance.seplocatent[1,2]
loc.concordance.seplocatent[1,3] <- as.numeric(CompEntPred.pair(ent.Raes12.2e.sep.by.loc.RE,ent.Raes34.2e.sep.by.loc.RE)$Maximum_number_samples_clustering_together[,1])
loc.discordance.seplocatent[1,3] <- total.RE - loc.concordance.seplocatent[1,3]
# 2 row = 12vs56_2e
loc.concordance.seplocatent[2,1] <- as.numeric(CompEntPred.pair(ent.Raes12.2e.sep.by.loc.IL,ent.Raes56.2e.sep.by.loc.IL)$Maximum_number_samples_clustering_together[,1])
loc.discordance.seplocatent[2,1] <- total.IL - loc.concordance.seplocatent[2,1]
loc.concordance.seplocatent[2,2] <- as.numeric(CompEntPred.pair(ent.Raes12.2e.sep.by.loc.TR,ent.Raes56.2e.sep.by.loc.TR)$Maximum_number_samples_clustering_together[,1])
loc.discordance.seplocatent[2,2] <- total.TR - loc.concordance.seplocatent[2,2]
loc.concordance.seplocatent[2,3] <- as.numeric(CompEntPred.pair(ent.Raes12.2e.sep.by.loc.RE,ent.Raes56.2e.sep.by.loc.RE)$Maximum_number_samples_clustering_together[,1])
loc.discordance.seplocatent[2,3] <- total.RE - loc.concordance.seplocatent[2,3]
# 3 row = 34vs56_2e 
loc.concordance.seplocatent[3,1] <- as.numeric(CompEntPred.pair(ent.Raes34.2e.sep.by.loc.IL,ent.Raes56.2e.sep.by.loc.IL)$Maximum_number_samples_clustering_together[,1])
loc.discordance.seplocatent[3,1] <- total.IL - loc.concordance.seplocatent[3,1]
loc.concordance.seplocatent[3,2] <- as.numeric(CompEntPred.pair(ent.Raes34.2e.sep.by.loc.TR,ent.Raes56.2e.sep.by.loc.TR)$Maximum_number_samples_clustering_together[,1])
loc.discordance.seplocatent[3,2] <- total.TR - loc.concordance.seplocatent[3,2]
loc.concordance.seplocatent[3,3] <- as.numeric(CompEntPred.pair(ent.Raes34.2e.sep.by.loc.RE,ent.Raes56.2e.sep.by.loc.RE)$Maximum_number_samples_clustering_together[,1])
loc.discordance.seplocatent[3,3] <- total.RE - loc.concordance.seplocatent[3,3]
# 4 row = 12vs34_3e 
loc.concordance.seplocatent[4,1] <- as.numeric(CompEntPred.pair(ent.Raes12.3e.sep.by.loc.IL,ent.Raes34.3e.sep.by.loc.IL)$Maximum_number_samples_clustering_together[,1])
loc.discordance.seplocatent[4,1] <- total.IL - loc.concordance.seplocatent[4,1]
loc.concordance.seplocatent[4,2] <- as.numeric(CompEntPred.pair(ent.Raes12.3e.sep.by.loc.TR,ent.Raes34.3e.sep.by.loc.TR)$Maximum_number_samples_clustering_together[,1])
loc.discordance.seplocatent[4,2] <- total.TR - loc.concordance.seplocatent[4,2]
loc.concordance.seplocatent[4,3] <- as.numeric(CompEntPred.pair(ent.Raes12.3e.sep.by.loc.RE,ent.Raes34.3e.sep.by.loc.RE)$Maximum_number_samples_clustering_together[,1])
loc.discordance.seplocatent[4,3] <- total.RE - loc.concordance.seplocatent[4,3]
# 5 row = 12vs56_3e
loc.concordance.seplocatent[5,1] <- as.numeric(CompEntPred.pair(ent.Raes12.3e.sep.by.loc.IL,ent.Raes56.3e.sep.by.loc.IL)$Maximum_number_samples_clustering_together[,1])
loc.discordance.seplocatent[5,1] <- total.IL - loc.concordance.seplocatent[5,1]
loc.concordance.seplocatent[5,2] <- as.numeric(CompEntPred.pair(ent.Raes12.3e.sep.by.loc.TR,ent.Raes56.3e.sep.by.loc.TR)$Maximum_number_samples_clustering_together[,1])
loc.discordance.seplocatent[5,2] <- total.TR - loc.concordance.seplocatent[5,2]
loc.concordance.seplocatent[5,3] <- as.numeric(CompEntPred.pair(ent.Raes12.3e.sep.by.loc.RE,ent.Raes56.3e.sep.by.loc.RE)$Maximum_number_samples_clustering_together[,1])
loc.discordance.seplocatent[5,3] <- total.RE - loc.concordance.seplocatent[5,3]
# 6 row = 34vs56_3e 
loc.concordance.seplocatent[6,1] <- as.numeric(CompEntPred.pair(ent.Raes34.3e.sep.by.loc.IL,ent.Raes56.3e.sep.by.loc.IL)$Maximum_number_samples_clustering_together[,1])
loc.discordance.seplocatent[6,1] <- total.IL - loc.concordance.seplocatent[6,1]
loc.concordance.seplocatent[6,2] <- as.numeric(CompEntPred.pair(ent.Raes34.3e.sep.by.loc.TR,ent.Raes56.3e.sep.by.loc.TR)$Maximum_number_samples_clustering_together[,1])
loc.discordance.seplocatent[6,2] <- total.TR - loc.concordance.seplocatent[6,2]
loc.concordance.seplocatent[6,3] <- as.numeric(CompEntPred.pair(ent.Raes34.3e.sep.by.loc.RE,ent.Raes56.3e.sep.by.loc.RE)$Maximum_number_samples_clustering_together[,1])
loc.discordance.seplocatent[6,3] <- total.RE - loc.concordance.seplocatent[6,3]
# 7 row = 12vs34_4e 
loc.concordance.seplocatent[7,1] <- as.numeric(CompEntPred.pair(ent.Raes12.4e.sep.by.loc.IL,ent.Raes34.4e.sep.by.loc.IL)$Maximum_number_samples_clustering_together[,1])
loc.discordance.seplocatent[7,1] <- total.IL - loc.concordance.seplocatent[7,1]
loc.concordance.seplocatent[7,2] <- as.numeric(CompEntPred.pair(ent.Raes12.4e.sep.by.loc.TR,ent.Raes34.4e.sep.by.loc.TR)$Maximum_number_samples_clustering_together[,1])
loc.discordance.seplocatent[7,2] <- total.TR - loc.concordance.seplocatent[7,2]
loc.concordance.seplocatent[7,3] <- as.numeric(CompEntPred.pair(ent.Raes12.4e.sep.by.loc.RE,ent.Raes34.4e.sep.by.loc.RE)$Maximum_number_samples_clustering_together[,1])
loc.discordance.seplocatent[7,3] <- total.RE - loc.concordance.seplocatent[7,3]
# 8 row = 12vs56_4e
loc.concordance.seplocatent[8,1] <- as.numeric(CompEntPred.pair(ent.Raes12.4e.sep.by.loc.IL,ent.Raes56.4e.sep.by.loc.IL)$Maximum_number_samples_clustering_together[,1])
loc.discordance.seplocatent[8,1] <- total.IL - loc.concordance.seplocatent[8,1]
loc.concordance.seplocatent[8,2] <- as.numeric(CompEntPred.pair(ent.Raes12.4e.sep.by.loc.TR,ent.Raes56.4e.sep.by.loc.TR)$Maximum_number_samples_clustering_together[,1])
loc.discordance.seplocatent[8,2] <- total.TR - loc.concordance.seplocatent[8,2]
loc.concordance.seplocatent[8,3] <- as.numeric(CompEntPred.pair(ent.Raes12.4e.sep.by.loc.RE,ent.Raes56.4e.sep.by.loc.RE)$Maximum_number_samples_clustering_together[,1])
loc.discordance.seplocatent[8,3] <- total.RE - loc.concordance.seplocatent[8,3]
# 9 row = 34vs56_4e 
loc.concordance.seplocatent[9,1] <- as.numeric(CompEntPred.pair(ent.Raes34.4e.sep.by.loc.IL,ent.Raes56.4e.sep.by.loc.IL)$Maximum_number_samples_clustering_together[,1])
loc.discordance.seplocatent[9,1] <- total.IL - loc.concordance.seplocatent[9,1]
loc.concordance.seplocatent[9,2] <- as.numeric(CompEntPred.pair(ent.Raes34.4e.sep.by.loc.TR,ent.Raes56.4e.sep.by.loc.TR)$Maximum_number_samples_clustering_together[,1])
loc.discordance.seplocatent[9,2] <- total.TR - loc.concordance.seplocatent[9,2]
loc.concordance.seplocatent[9,3] <- as.numeric(CompEntPred.pair(ent.Raes34.4e.sep.by.loc.RE,ent.Raes56.4e.sep.by.loc.RE)$Maximum_number_samples_clustering_together[,1])
loc.discordance.seplocatent[9,3] <- total.RE - loc.concordance.seplocatent[9,3]
# 10 row = 12vs34_5e 
loc.concordance.seplocatent[10,1] <- as.numeric(CompEntPred.pair(ent.Raes12.5e.sep.by.loc.IL,ent.Raes34.5e.sep.by.loc.IL)$Maximum_number_samples_clustering_together[,1])
loc.discordance.seplocatent[10,1] <- total.IL - loc.concordance.seplocatent[10,1]
loc.concordance.seplocatent[10,2] <- as.numeric(CompEntPred.pair(ent.Raes12.5e.sep.by.loc.TR,ent.Raes34.5e.sep.by.loc.TR)$Maximum_number_samples_clustering_together[,1])
loc.discordance.seplocatent[10,2] <- total.TR - loc.concordance.seplocatent[10,2]
loc.concordance.seplocatent[10,3] <- as.numeric(CompEntPred.pair(ent.Raes12.5e.sep.by.loc.RE,ent.Raes34.5e.sep.by.loc.RE)$Maximum_number_samples_clustering_together[,1])
loc.discordance.seplocatent[10,3] <- total.RE - loc.concordance.seplocatent[10,3]
# 11 row = 12vs56_5e
loc.concordance.seplocatent[11,1] <- as.numeric(CompEntPred.pair(ent.Raes12.5e.sep.by.loc.IL,ent.Raes56.5e.sep.by.loc.IL)$Maximum_number_samples_clustering_together[,1])
loc.discordance.seplocatent[11,1] <- total.IL - loc.concordance.seplocatent[11,1]
loc.concordance.seplocatent[11,2] <- as.numeric(CompEntPred.pair(ent.Raes12.5e.sep.by.loc.TR,ent.Raes56.5e.sep.by.loc.TR)$Maximum_number_samples_clustering_together[,1])
loc.discordance.seplocatent[11,2] <- total.TR - loc.concordance.seplocatent[11,2]
loc.concordance.seplocatent[11,3] <- as.numeric(CompEntPred.pair(ent.Raes12.5e.sep.by.loc.RE,ent.Raes56.5e.sep.by.loc.RE)$Maximum_number_samples_clustering_together[,1])
loc.discordance.seplocatent[11,3] <- total.RE - loc.concordance.seplocatent[11,3]
# 12 row = 34vs56_5e 
loc.concordance.seplocatent[12,1] <- as.numeric(CompEntPred.pair(ent.Raes34.5e.sep.by.loc.IL,ent.Raes56.5e.sep.by.loc.IL)$Maximum_number_samples_clustering_together[,1])
loc.discordance.seplocatent[12,1] <- total.IL - loc.concordance.seplocatent[12,1]
loc.concordance.seplocatent[12,2] <- as.numeric(CompEntPred.pair(ent.Raes34.5e.sep.by.loc.TR,ent.Raes56.5e.sep.by.loc.TR)$Maximum_number_samples_clustering_together[,1])
loc.discordance.seplocatent[12,2] <- total.TR - loc.concordance.seplocatent[12,2]
loc.concordance.seplocatent[12,3] <- as.numeric(CompEntPred.pair(ent.Raes34.5e.sep.by.loc.RE,ent.Raes56.5e.sep.by.loc.RE)$Maximum_number_samples_clustering_together[,1])
loc.discordance.seplocatent[12,3] <- total.RE - loc.concordance.seplocatent[12,3]

# Now make a matrix of SUM of concordance and discordance. This is total number of samples per location
total.samples <- matrix()
total.samples <- loc.concordance.seplocatent+loc.discordance.seplocatent
# Step 3. Calculate the sum of corresponding samples in a certain location divided by sum of non-corresponding samples
result.IL <- sum(loc.concordance.seplocatent[,1])/sum(loc.discordance.seplocatent[,1])
result.TR <- sum(loc.concordance.seplocatent[,2])/sum(loc.discordance.seplocatent[,2])
result.RE <- sum(loc.concordance.seplocatent[,3])/sum(loc.discordance.seplocatent[,3])

# Step 4.
# Now, simulate the random process of distr between locations in case IL and TR have the same 
# probability of concordance equal to their average

# First, we need a matrix of average percentage of sample concordance across locations (i.e. (concordant_IL+concordant_TR)/total(IL+TR))
probab.matrix <- matrix(0,12,3)
colnames(probab.matrix) <- c('IL_TR','IL_RE','TR_RE')
rownames(probab.matrix) <- c('v12-v34_2ent','v12-v56_2ent','v34-v56_2ent','v12-v34_3ent','v12-v56_3ent','v34-v56_3ent','v12-v34_4ent','v12-v56_4ent','v34-v56_4ent','v12-v34_5ent','v12-v56_5ent','v34-v56_5ent')
probab.matrix[1,1] <- (loc.concordance.seplocatent[1,1]+loc.concordance.seplocatent[1,2])/(total.IL+total.TR)
probab.matrix[1,2] <- (loc.concordance.seplocatent[1,1]+loc.concordance.seplocatent[1,3])/(total.IL+total.RE)
probab.matrix[1,3] <- (loc.concordance.seplocatent[1,2]+loc.concordance.seplocatent[1,3])/(total.TR+total.RE)

probab.matrix[2,1] <- (loc.concordance.seplocatent[2,1]+loc.concordance.seplocatent[2,2])/(total.IL+total.TR)
probab.matrix[2,2] <- (loc.concordance.seplocatent[2,1]+loc.concordance.seplocatent[2,3])/(total.IL+total.RE)
probab.matrix[2,3] <- (loc.concordance.seplocatent[2,2]+loc.concordance.seplocatent[2,3])/(total.TR+total.RE)

probab.matrix[3,1] <- (loc.concordance.seplocatent[3,1]+loc.concordance.seplocatent[3,2])/(total.IL+total.TR)
probab.matrix[3,2] <- (loc.concordance.seplocatent[3,1]+loc.concordance.seplocatent[3,3])/(total.IL+total.RE)
probab.matrix[3,3] <- (loc.concordance.seplocatent[3,2]+loc.concordance.seplocatent[3,3])/(total.TR+total.RE)

probab.matrix[4,1] <- (loc.concordance.seplocatent[4,1]+loc.concordance.seplocatent[4,2])/(total.IL+total.TR)
probab.matrix[4,2] <- (loc.concordance.seplocatent[4,1]+loc.concordance.seplocatent[4,3])/(total.IL+total.RE)
probab.matrix[4,3] <- (loc.concordance.seplocatent[4,2]+loc.concordance.seplocatent[4,3])/(total.TR+total.RE)

probab.matrix[5,1] <- (loc.concordance.seplocatent[5,1]+loc.concordance.seplocatent[5,2])/(total.IL+total.TR)
probab.matrix[5,2] <- (loc.concordance.seplocatent[5,1]+loc.concordance.seplocatent[5,3])/(total.IL+total.RE)
probab.matrix[5,3] <- (loc.concordance.seplocatent[5,2]+loc.concordance.seplocatent[5,3])/(total.TR+total.RE)

probab.matrix[6,1] <- (loc.concordance.seplocatent[6,1]+loc.concordance.seplocatent[6,2])/(total.IL+total.TR)
probab.matrix[6,2] <- (loc.concordance.seplocatent[6,1]+loc.concordance.seplocatent[6,3])/(total.IL+total.RE)
probab.matrix[6,3] <- (loc.concordance.seplocatent[6,2]+loc.concordance.seplocatent[6,3])/(total.TR+total.RE)

probab.matrix[7,1] <- (loc.concordance.seplocatent[7,1]+loc.concordance.seplocatent[7,2])/(total.IL+total.TR)
probab.matrix[7,2] <- (loc.concordance.seplocatent[7,1]+loc.concordance.seplocatent[7,3])/(total.IL+total.RE)
probab.matrix[7,3] <- (loc.concordance.seplocatent[7,2]+loc.concordance.seplocatent[7,3])/(total.TR+total.RE)

probab.matrix[8,1] <- (loc.concordance.seplocatent[8,1]+loc.concordance.seplocatent[8,2])/(total.IL+total.TR)
probab.matrix[8,2] <- (loc.concordance.seplocatent[8,1]+loc.concordance.seplocatent[8,3])/(total.IL+total.RE)
probab.matrix[8,3] <- (loc.concordance.seplocatent[8,2]+loc.concordance.seplocatent[8,3])/(total.TR+total.RE)

probab.matrix[9,1] <- (loc.concordance.seplocatent[9,1]+loc.concordance.seplocatent[9,2])/(total.IL+total.TR)
probab.matrix[9,2] <- (loc.concordance.seplocatent[9,1]+loc.concordance.seplocatent[9,3])/(total.IL+total.RE)
probab.matrix[9,3] <- (loc.concordance.seplocatent[9,2]+loc.concordance.seplocatent[9,3])/(total.TR+total.RE)

probab.matrix[10,1] <- (loc.concordance.seplocatent[10,1]+loc.concordance.seplocatent[10,2])/(total.IL+total.TR)
probab.matrix[10,2] <- (loc.concordance.seplocatent[10,1]+loc.concordance.seplocatent[10,3])/(total.IL+total.RE)
probab.matrix[10,3] <- (loc.concordance.seplocatent[10,2]+loc.concordance.seplocatent[10,3])/(total.TR+total.RE)

probab.matrix[11,1] <- (loc.concordance.seplocatent[11,1]+loc.concordance.seplocatent[11,2])/(total.IL+total.TR)
probab.matrix[11,2] <- (loc.concordance.seplocatent[11,1]+loc.concordance.seplocatent[11,3])/(total.IL+total.RE)
probab.matrix[11,3] <- (loc.concordance.seplocatent[11,2]+loc.concordance.seplocatent[11,3])/(total.TR+total.RE)

probab.matrix[12,1] <- (loc.concordance.seplocatent[12,1]+loc.concordance.seplocatent[12,2])/(total.IL+total.TR)
probab.matrix[12,2] <- (loc.concordance.seplocatent[12,1]+loc.concordance.seplocatent[12,3])/(total.IL+total.RE)
probab.matrix[12,3] <- (loc.concordance.seplocatent[12,2]+loc.concordance.seplocatent[12,3])/(total.TR+total.RE)


# And now, simulate the similarity of samples between locations. So, that we assume e.g. that the probability of a sample being
# similar in v12_2ent IL and v34_2ent IL is equal to the probability of a sample being similar in v12_2ent TR and v34_2ent TR 
# and the numbers are equal to probab.matrix[1,1]

# stage a. We simulate IL vs TR based on probab.matrix[,1] probabilities
# we assume that the probability of concordance of v12_2e_IL with v34_2e_IL = probability of concordance of  v12_2e_TR with v34_2e_TR=
# = probab.matrix[1,1]
simul.loc.concordance.seplocatent <- matrix(0,12,2)
colnames(simul.loc.concordance.seplocatent) <- c('IL','TR')
rownames(simul.loc.concordance.seplocatent) <- c('v12-v34_2ent','v12-v56_2ent','v34-v56_2ent','v12-v34_3ent','v12-v56_3ent','v34-v56_3ent','v12-v34_4ent','v12-v56_4ent','v34-v56_4ent','v12-v34_5ent','v12-v56_5ent','v34-v56_5ent')
# Now, make N simulations
# metrika is a vector of length N with simulations of  
# number_sampl_the_same_IL_v12_with_v34_2ent/number_sampl_different_IL_v12_with_v34_2ent - number_sampl_the_same_TR_v12_with_v34_2ent/number_sampl_different_TR_v12_with_v34_2ent
metrika_IL_TR <- c(rep(0,N))
for (i in 1:N){
simul.loc.concordance.seplocatent[1,1] <- sum(sample(c(1,0),size = total.IL, replace=T, prob = c(probab.matrix[1,1],1-probab.matrix[1,1])))
simul.loc.concordance.seplocatent[2,1] <- sum(sample(c(1,0),size = total.IL, replace=T, prob = c(probab.matrix[2,1],1-probab.matrix[2,1])))
simul.loc.concordance.seplocatent[3,1] <- sum(sample(c(1,0),size = total.IL, replace=T, prob = c(probab.matrix[3,1],1-probab.matrix[3,1])))
simul.loc.concordance.seplocatent[4,1] <- sum(sample(c(1,0),size = total.IL, replace=T, prob = c(probab.matrix[4,1],1-probab.matrix[4,1])))
simul.loc.concordance.seplocatent[5,1] <- sum(sample(c(1,0),size = total.IL, replace=T, prob = c(probab.matrix[5,1],1-probab.matrix[5,1])))
simul.loc.concordance.seplocatent[6,1] <- sum(sample(c(1,0),size = total.IL, replace=T, prob = c(probab.matrix[6,1],1-probab.matrix[6,1])))
simul.loc.concordance.seplocatent[7,1] <- sum(sample(c(1,0),size = total.IL, replace=T, prob = c(probab.matrix[7,1],1-probab.matrix[7,1])))
simul.loc.concordance.seplocatent[8,1] <- sum(sample(c(1,0),size = total.IL, replace=T, prob = c(probab.matrix[8,1],1-probab.matrix[8,1])))
simul.loc.concordance.seplocatent[9,1] <- sum(sample(c(1,0),size = total.IL, replace=T, prob = c(probab.matrix[9,1],1-probab.matrix[9,1])))
simul.loc.concordance.seplocatent[10,1] <- sum(sample(c(1,0),size = total.IL, replace=T, prob = c(probab.matrix[10,1],1-probab.matrix[10,1])))
simul.loc.concordance.seplocatent[11,1] <- sum(sample(c(1,0),size = total.IL, replace=T, prob = c(probab.matrix[11,1],1-probab.matrix[11,1])))
simul.loc.concordance.seplocatent[12,1] <- sum(sample(c(1,0),size = total.IL, replace=T, prob = c(probab.matrix[12,1],1-probab.matrix[12,1])))

# than similarly make column 2
simul.loc.concordance.seplocatent[1,2] <- sum(sample(c(1,0),size = total.TR, replace=T, prob = c(probab.matrix[1,1],1-probab.matrix[1,1])))
simul.loc.concordance.seplocatent[2,2] <- sum(sample(c(1,0),size = total.TR, replace=T, prob = c(probab.matrix[2,1],1-probab.matrix[2,1])))
simul.loc.concordance.seplocatent[3,2] <- sum(sample(c(1,0),size = total.TR, replace=T, prob = c(probab.matrix[3,1],1-probab.matrix[3,1])))
simul.loc.concordance.seplocatent[4,2] <- sum(sample(c(1,0),size = total.TR, replace=T, prob = c(probab.matrix[4,1],1-probab.matrix[4,1])))
simul.loc.concordance.seplocatent[5,2] <- sum(sample(c(1,0),size = total.TR, replace=T, prob = c(probab.matrix[5,1],1-probab.matrix[5,1])))
simul.loc.concordance.seplocatent[6,2] <- sum(sample(c(1,0),size = total.TR, replace=T, prob = c(probab.matrix[6,1],1-probab.matrix[6,1])))
simul.loc.concordance.seplocatent[7,2] <- sum(sample(c(1,0),size = total.TR, replace=T, prob = c(probab.matrix[7,1],1-probab.matrix[7,1])))
simul.loc.concordance.seplocatent[8,2] <- sum(sample(c(1,0),size = total.TR, replace=T, prob = c(probab.matrix[8,1],1-probab.matrix[8,1])))
simul.loc.concordance.seplocatent[9,2] <- sum(sample(c(1,0),size = total.TR, replace=T, prob = c(probab.matrix[9,1],1-probab.matrix[9,1])))
simul.loc.concordance.seplocatent[10,2] <- sum(sample(c(1,0),size = total.TR, replace=T, prob = c(probab.matrix[10,1],1-probab.matrix[10,1])))
simul.loc.concordance.seplocatent[11,2] <- sum(sample(c(1,0),size = total.TR, replace=T, prob = c(probab.matrix[11,1],1-probab.matrix[11,1])))
simul.loc.concordance.seplocatent[12,2] <- sum(sample(c(1,0),size = total.TR, replace=T, prob = c(probab.matrix[12,1],1-probab.matrix[12,1])))

# Now prepare the discordance matrix based on the concordance and the total number of samples
simul.loc.discordance.seplocatent <- total.samples[,1:2]-simul.loc.concordance.seplocatent

# get the metrics
metrika_IL_TR[i] <- sum(simul.loc.concordance.seplocatent[,1])/sum(simul.loc.discordance.seplocatent[,1])-sum(simul.loc.concordance.seplocatent[,2])/sum(simul.loc.discordance.seplocatent[,2])
}

metrika_IL_TR <- sort(metrika_IL_TR)
# now, define the percentile to which our result corresponds in the random simulation
difference_IL_TR <- result.IL-result.TR
i=1
while (metrika_IL_TR[i]<difference_IL_TR){
 i = i+1
}
sim.diff.IL.TR <- i/N
paste('The p-value for null hypothesis that there is no difference between IL and TR enterotype stability is: ',signif(100*i/N,digits = 3),'%','(a total of ',N,' simulations done)')




# stage b. We simulate IL vs RE based on probab.matrix[,2] probabilities and check if RE enterotype is significantly more stable than IL
# we assume that the probability of concordance of v12_2e_IL with v34_2e_IL = probability of concordance of  v12_2e_RE with v34_2e_RE=
# = probab.matrix[1,2]
# simul.loc.concordance.seplocatent - we re-use this matrix for simulation in this location
simul.loc.concordance.seplocatent <- matrix(0,12,2)
colnames(simul.loc.concordance.seplocatent) <- c('IL','RE')
rownames(simul.loc.concordance.seplocatent) <- c('v12-v34_2ent','v12-v56_2ent','v34-v56_2ent','v12-v34_3ent','v12-v56_3ent','v34-v56_3ent','v12-v34_4ent','v12-v56_4ent','v34-v56_4ent','v12-v34_5ent','v12-v56_5ent','v34-v56_5ent')
# Now, make N simulations
# metrika is a vector of length N with simulations of  
# number_sampl_the_same_IL_v12_with_v34_2ent/number_sampl_different_IL_v12_with_v34_2ent - number_sampl_the_same_RE_v12_with_v34_2ent/number_sampl_different_RE_v12_with_v34_2ent
metrika_IL_RE <- c(rep(0,N))
for (i in 1:N){
  simul.loc.concordance.seplocatent[1,1] <- sum(sample(c(1,0),size = total.IL, replace=T, prob = c(probab.matrix[1,2],1-probab.matrix[1,2])))
  simul.loc.concordance.seplocatent[2,1] <- sum(sample(c(1,0),size = total.IL, replace=T, prob = c(probab.matrix[2,2],1-probab.matrix[2,2])))
  simul.loc.concordance.seplocatent[3,1] <- sum(sample(c(1,0),size = total.IL, replace=T, prob = c(probab.matrix[3,2],1-probab.matrix[3,2])))
  simul.loc.concordance.seplocatent[4,1] <- sum(sample(c(1,0),size = total.IL, replace=T, prob = c(probab.matrix[4,2],1-probab.matrix[4,2])))
  simul.loc.concordance.seplocatent[5,1] <- sum(sample(c(1,0),size = total.IL, replace=T, prob = c(probab.matrix[5,2],1-probab.matrix[5,2])))
  simul.loc.concordance.seplocatent[6,1] <- sum(sample(c(1,0),size = total.IL, replace=T, prob = c(probab.matrix[6,2],1-probab.matrix[6,2])))
  simul.loc.concordance.seplocatent[7,1] <- sum(sample(c(1,0),size = total.IL, replace=T, prob = c(probab.matrix[7,2],1-probab.matrix[7,2])))
  simul.loc.concordance.seplocatent[8,1] <- sum(sample(c(1,0),size = total.IL, replace=T, prob = c(probab.matrix[8,2],1-probab.matrix[8,2])))
  simul.loc.concordance.seplocatent[9,1] <- sum(sample(c(1,0),size = total.IL, replace=T, prob = c(probab.matrix[9,2],1-probab.matrix[9,2])))
  simul.loc.concordance.seplocatent[10,1] <- sum(sample(c(1,0),size = total.IL, replace=T, prob = c(probab.matrix[10,2],1-probab.matrix[10,2])))
  simul.loc.concordance.seplocatent[11,1] <- sum(sample(c(1,0),size = total.IL, replace=T, prob = c(probab.matrix[11,2],1-probab.matrix[11,2])))
  simul.loc.concordance.seplocatent[12,1] <- sum(sample(c(1,0),size = total.IL, replace=T, prob = c(probab.matrix[12,2],1-probab.matrix[12,2])))
  
  # than similarly make column 2
  simul.loc.concordance.seplocatent[1,2] <- sum(sample(c(1,0),size = total.RE, replace=T, prob = c(probab.matrix[1,2],1-probab.matrix[1,2])))
  simul.loc.concordance.seplocatent[2,2] <- sum(sample(c(1,0),size = total.RE, replace=T, prob = c(probab.matrix[2,2],1-probab.matrix[2,2])))
  simul.loc.concordance.seplocatent[3,2] <- sum(sample(c(1,0),size = total.RE, replace=T, prob = c(probab.matrix[3,2],1-probab.matrix[3,2])))
  simul.loc.concordance.seplocatent[4,2] <- sum(sample(c(1,0),size = total.RE, replace=T, prob = c(probab.matrix[4,2],1-probab.matrix[4,2])))
  simul.loc.concordance.seplocatent[5,2] <- sum(sample(c(1,0),size = total.RE, replace=T, prob = c(probab.matrix[5,2],1-probab.matrix[5,2])))
  simul.loc.concordance.seplocatent[6,2] <- sum(sample(c(1,0),size = total.RE, replace=T, prob = c(probab.matrix[6,2],1-probab.matrix[6,2])))
  simul.loc.concordance.seplocatent[7,2] <- sum(sample(c(1,0),size = total.RE, replace=T, prob = c(probab.matrix[7,2],1-probab.matrix[7,2])))
  simul.loc.concordance.seplocatent[8,2] <- sum(sample(c(1,0),size = total.RE, replace=T, prob = c(probab.matrix[8,2],1-probab.matrix[8,2])))
  simul.loc.concordance.seplocatent[9,2] <- sum(sample(c(1,0),size = total.RE, replace=T, prob = c(probab.matrix[9,2],1-probab.matrix[9,2])))
  simul.loc.concordance.seplocatent[10,2] <- sum(sample(c(1,0),size = total.RE, replace=T, prob = c(probab.matrix[10,2],1-probab.matrix[10,2])))
  simul.loc.concordance.seplocatent[11,2] <- sum(sample(c(1,0),size = total.RE, replace=T, prob = c(probab.matrix[11,2],1-probab.matrix[11,2])))
  simul.loc.concordance.seplocatent[12,2] <- sum(sample(c(1,0),size = total.RE, replace=T, prob = c(probab.matrix[12,2],1-probab.matrix[12,2])))
  
  # Now prepare the discordance matrix based on the concordance and the total number of samples
  simul.loc.discordance.seplocatent <- total.samples[,c(1,3)]-simul.loc.concordance.seplocatent
  
  # get the metrics
  metrika_IL_RE[i] <- sum(simul.loc.concordance.seplocatent[,1])/sum(simul.loc.discordance.seplocatent[,1])-sum(simul.loc.concordance.seplocatent[,2])/sum(simul.loc.discordance.seplocatent[,2])
}

metrika_IL_RE <- sort(metrika_IL_RE)
# now, define the percentile to which our result corresponds in the random simulation
difference_IL_RE <- result.IL-result.RE
i=1
while (metrika_IL_RE[i]<difference_IL_RE){
  i = i+1
}

sim.diff.IL.RE <- i/N
paste('The p-value for null hypothesis that there is no difference between IL and RE enterotype stability is: ',signif(100*i/N,digits = 3),'%','(a total of ',N,' simulations done)')


# stage c. We simulate TR vs RE based on probab.matrix[,3] probabilities and check if RE enterotype is significantly more stable than TR
# we assume that the probability of concordance of v12_2e_TR with v34_2e_TR = probability of concordance of  v12_2e_RE with v34_2e_RE=
# = probab.matrix[1,3]
# simul.loc.concordance.seplocatent - we re-use this matrix for simulation in this location
simul.loc.concordance.seplocatent <- matrix(0,12,2)
colnames(simul.loc.concordance.seplocatent) <- c('TR','RE')
rownames(simul.loc.concordance.seplocatent) <- c('v12-v34_2ent','v12-v56_2ent','v34-v56_2ent','v12-v34_3ent','v12-v56_3ent','v34-v56_3ent','v12-v34_4ent','v12-v56_4ent','v34-v56_4ent','v12-v34_5ent','v12-v56_5ent','v34-v56_5ent')
# Now, make N simulations
# metrika is a vector of length N with simulations of  
# number_sampl_the_same_TR_v12_with_v34_2ent/number_sampl_different_TR_v12_with_v34_2ent - number_sampl_the_same_RE_v12_with_v34_2ent/number_sampl_different_RE_v12_with_v34_2ent
metrika_TR_RE <- c(rep(0,N))
for (i in 1:N){
  simul.loc.concordance.seplocatent[1,1] <- sum(sample(c(1,0),size = total.TR, replace=T, prob = c(probab.matrix[1,3],1-probab.matrix[1,3])))
  simul.loc.concordance.seplocatent[2,1] <- sum(sample(c(1,0),size = total.TR, replace=T, prob = c(probab.matrix[2,3],1-probab.matrix[2,3])))
  simul.loc.concordance.seplocatent[3,1] <- sum(sample(c(1,0),size = total.TR, replace=T, prob = c(probab.matrix[3,3],1-probab.matrix[3,3])))
  simul.loc.concordance.seplocatent[4,1] <- sum(sample(c(1,0),size = total.TR, replace=T, prob = c(probab.matrix[4,3],1-probab.matrix[4,3])))
  simul.loc.concordance.seplocatent[5,1] <- sum(sample(c(1,0),size = total.TR, replace=T, prob = c(probab.matrix[5,3],1-probab.matrix[5,3])))
  simul.loc.concordance.seplocatent[6,1] <- sum(sample(c(1,0),size = total.TR, replace=T, prob = c(probab.matrix[6,3],1-probab.matrix[6,3])))
  simul.loc.concordance.seplocatent[7,1] <- sum(sample(c(1,0),size = total.TR, replace=T, prob = c(probab.matrix[7,3],1-probab.matrix[7,3])))
  simul.loc.concordance.seplocatent[8,1] <- sum(sample(c(1,0),size = total.TR, replace=T, prob = c(probab.matrix[8,3],1-probab.matrix[8,3])))
  simul.loc.concordance.seplocatent[9,1] <- sum(sample(c(1,0),size = total.TR, replace=T, prob = c(probab.matrix[9,3],1-probab.matrix[9,3])))
  simul.loc.concordance.seplocatent[10,1] <- sum(sample(c(1,0),size = total.TR, replace=T, prob = c(probab.matrix[10,3],1-probab.matrix[10,3])))
  simul.loc.concordance.seplocatent[11,1] <- sum(sample(c(1,0),size = total.TR, replace=T, prob = c(probab.matrix[11,3],1-probab.matrix[11,3])))
  simul.loc.concordance.seplocatent[12,1] <- sum(sample(c(1,0),size = total.TR, replace=T, prob = c(probab.matrix[12,3],1-probab.matrix[12,3])))
  
  # than similarly make column 2
  simul.loc.concordance.seplocatent[1,2] <- sum(sample(c(1,0),size = total.RE, replace=T, prob = c(probab.matrix[1,3],1-probab.matrix[1,3])))
  simul.loc.concordance.seplocatent[2,2] <- sum(sample(c(1,0),size = total.RE, replace=T, prob = c(probab.matrix[2,3],1-probab.matrix[2,3])))
  simul.loc.concordance.seplocatent[3,2] <- sum(sample(c(1,0),size = total.RE, replace=T, prob = c(probab.matrix[3,3],1-probab.matrix[3,3])))
  simul.loc.concordance.seplocatent[4,2] <- sum(sample(c(1,0),size = total.RE, replace=T, prob = c(probab.matrix[4,3],1-probab.matrix[4,3])))
  simul.loc.concordance.seplocatent[5,2] <- sum(sample(c(1,0),size = total.RE, replace=T, prob = c(probab.matrix[5,3],1-probab.matrix[5,3])))
  simul.loc.concordance.seplocatent[6,2] <- sum(sample(c(1,0),size = total.RE, replace=T, prob = c(probab.matrix[6,3],1-probab.matrix[6,3])))
  simul.loc.concordance.seplocatent[7,2] <- sum(sample(c(1,0),size = total.RE, replace=T, prob = c(probab.matrix[7,3],1-probab.matrix[7,3])))
  simul.loc.concordance.seplocatent[8,2] <- sum(sample(c(1,0),size = total.RE, replace=T, prob = c(probab.matrix[8,3],1-probab.matrix[8,3])))
  simul.loc.concordance.seplocatent[9,2] <- sum(sample(c(1,0),size = total.RE, replace=T, prob = c(probab.matrix[9,3],1-probab.matrix[9,3])))
  simul.loc.concordance.seplocatent[10,2] <- sum(sample(c(1,0),size = total.RE, replace=T, prob = c(probab.matrix[10,3],1-probab.matrix[10,3])))
  simul.loc.concordance.seplocatent[11,2] <- sum(sample(c(1,0),size = total.RE, replace=T, prob = c(probab.matrix[11,3],1-probab.matrix[11,3])))
  simul.loc.concordance.seplocatent[12,2] <- sum(sample(c(1,0),size = total.RE, replace=T, prob = c(probab.matrix[12,3],1-probab.matrix[12,3])))
  
  # Now prepare the discordance matrix based on the concordance and the total number of samples
  simul.loc.discordance.seplocatent <- total.samples[,c(2,3)]-simul.loc.concordance.seplocatent
  
  # get the metrics
  metrika_TR_RE[i] <- sum(simul.loc.concordance.seplocatent[,1])/sum(simul.loc.discordance.seplocatent[,1])-sum(simul.loc.concordance.seplocatent[,2])/sum(simul.loc.discordance.seplocatent[,2])
}

metrika_TR_RE <- sort(metrika_TR_RE)
# now, define the percentile to which our result corresponds in the random simulation
difference_TR_RE <- result.TR-result.RE
i=1
while (metrika_TR_RE[i]<difference_TR_RE){
  i = i+1
}

sim.diff.TR.RE <- i/N
paste('The p-value for null hypothesis that there is no difference between TR and RE enterotype stability is: ',signif(100*i/N,digits = 3),'%','(a total of ',N,' simulations done)')
result <- c(sim.diff.IL.TR,sim.diff.IL.RE,sim.diff.TR.RE)
names(result) <- c('IL_TR','IL_RE','TR_RE')
return(result)
}