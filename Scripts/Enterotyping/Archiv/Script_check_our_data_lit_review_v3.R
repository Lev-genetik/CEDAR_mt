# This script is to check different things about our data that appear after literature review

# 1. CHECK FIRMICUTES According to Arumigam 2011, Costea 2018, the enterotypes they get are
# ET B lead by Bacteroides
# ET P lead by Prevotella
# ET F distinguished by overrepresentation of Firmicutes, most prominently Ruminococcus
# So, need to check if the 3rd enterotype is really distinguished by overrepresentation of Firmicutes
library(data.table)
library(stringr)
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Enterotyping_Raes_noise_removed_v9.R')

# make enterotyping
ent.Raes12.3e <- Enterotyping_Raes(file1, '/media/lev-genetik/980E73270E72FD96/Liege/Results/Clustering/V1V2', clust.number=3, noiserem=T,  locat = 'all', driv.number = 1000, plots=F)
ent.Raes34.3e <- Enterotyping_Raes(file2, '/media/lev-genetik/980E73270E72FD96/Liege/Results/Clustering/V3V4', clust.number=3, noiserem=T,  locat = 'all', driv.number = 1000, plots=F)
ent.Raes56.3e <- Enterotyping_Raes(file3, '/media/lev-genetik/980E73270E72FD96/Liege/Results/Clustering/V5V6', clust.number=3, noiserem=T,  locat = 'all', driv.number = 1000, plots=F)

# v1v2 enterotyping
# get sum of all taxa containing 'Firmicutes'
sum(ent.Raes12.3e$obs_bet_result$tab[1,colnames(ent.Raes12.3e$obs_bet_result$tab) %like% 'Firmicutes'])
sum(ent.Raes12.3e$obs_bet_result$tab[2,colnames(ent.Raes12.3e$obs_bet_result$tab) %like% 'Firmicutes'])
sum(ent.Raes12.3e$obs_bet_result$tab[3,colnames(ent.Raes12.3e$obs_bet_result$tab) %like% 'Firmicutes'])
# get value for Ruminococcus
ent.Raes12.3e$obs_bet_result$tab[1,colnames(ent.Raes12.3e$obs_bet_result$tab) %like% 'g__Ruminococcus']
ent.Raes12.3e$obs_bet_result$tab[2,colnames(ent.Raes12.3e$obs_bet_result$tab) %like% 'g__Ruminococcus']
ent.Raes12.3e$obs_bet_result$tab[3,colnames(ent.Raes12.3e$obs_bet_result$tab) %like% 'g__Ruminococcus']
# now see where is Bacteroides, Prevotella and the 3rd enterotype
which.max(ent.Raes12.3e$obs_bet_result$tab[1,])
which.max(ent.Raes12.3e$obs_bet_result$tab[2,])
which.max(ent.Raes12.3e$obs_bet_result$tab[3,])

# v3v4 enterotyping
# get sum of all taxa containing 'Firmicutes'
sum(ent.Raes34.3e$obs_bet_result$tab[1,colnames(ent.Raes34.3e$obs_bet_result$tab) %like% 'Firmicutes'])
sum(ent.Raes34.3e$obs_bet_result$tab[2,colnames(ent.Raes34.3e$obs_bet_result$tab) %like% 'Firmicutes'])
sum(ent.Raes34.3e$obs_bet_result$tab[3,colnames(ent.Raes34.3e$obs_bet_result$tab) %like% 'Firmicutes'])
# get value for Ruminococcus
ent.Raes34.3e$obs_bet_result$tab[1,colnames(ent.Raes34.3e$obs_bet_result$tab) %like% 'g__Ruminococcus']
ent.Raes34.3e$obs_bet_result$tab[2,colnames(ent.Raes34.3e$obs_bet_result$tab) %like% 'g__Ruminococcus']
ent.Raes34.3e$obs_bet_result$tab[3,colnames(ent.Raes34.3e$obs_bet_result$tab) %like% 'g__Ruminococcus']
# now see where is Bacteroides, Prevotella and the 3rd enterotype
which.max(ent.Raes34.3e$obs_bet_result$tab[1,])
which.max(ent.Raes34.3e$obs_bet_result$tab[2,])
which.max(ent.Raes34.3e$obs_bet_result$tab[3,])

# v5v6 enterotyping
# get sum of all taxa containing 'Firmicutes'
sum(ent.Raes56.3e$obs_bet_result$tab[1,colnames(ent.Raes56.3e$obs_bet_result$tab) %like% 'Firmicutes'])
sum(ent.Raes56.3e$obs_bet_result$tab[2,colnames(ent.Raes56.3e$obs_bet_result$tab) %like% 'Firmicutes'])
sum(ent.Raes56.3e$obs_bet_result$tab[3,colnames(ent.Raes56.3e$obs_bet_result$tab) %like% 'Firmicutes'])
# get value for Ruminococcus
ent.Raes56.3e$obs_bet_result$tab[1,colnames(ent.Raes56.3e$obs_bet_result$tab) %like% 'g__Ruminococcus']
ent.Raes56.3e$obs_bet_result$tab[2,colnames(ent.Raes56.3e$obs_bet_result$tab) %like% 'g__Ruminococcus']
ent.Raes56.3e$obs_bet_result$tab[3,colnames(ent.Raes56.3e$obs_bet_result$tab) %like% 'g__Ruminococcus']
# now see where is Bacteroides, Prevotella and the 3rd enterotype
which.max(ent.Raes56.3e$obs_bet_result$tab[1,])
which.max(ent.Raes56.3e$obs_bet_result$tab[2,])
which.max(ent.Raes56.3e$obs_bet_result$tab[3,])


# 2. WITHOUT BACTEROIDES AND PREVOTELLA ENTEROTYPING Check the hypothesis that without Bacteroides and Prevotella, samples would not cluster
# It seems they really cluster much worse
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Enterotyping_Raes_noise_removed_v9.R')
ent.Raes12.3e.noBactPrev <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6_g_Bact_Prev_removed_v2.txt', output_dir='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/g_Bact_Prev_removed/Graphs/V1V2', clust.number=3, noiserem=T, plots=T)
ent.Raes34.3e.noBactPrev <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6_g_Bact_Prev_removed_v2.txt', output_dir='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/g_Bact_Prev_removed/Graphs/V3V4', clust.number=3, noiserem=T, plots=T)
ent.Raes56.3e.noBactPrev <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6_g_Bact_Prev_removed_v2.txt', output_dir='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/g_Bact_Prev_removed/Graphs/V5V6', clust.number=3, noiserem=T, plots=T)

source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Enterotyping_compar_functions_v5.R')
# Now get percentage of samples with similar enterotypes between different amplicons
tmp12_34 <- CompEntPred.pair(ent.Raes34.3e.noBactPrev[[2]],ent.Raes12.3e.noBactPrev[[2]])
tmp34_56 <- CompEntPred.pair(ent.Raes34.3e.noBactPrev[[2]],ent.Raes56.3e.noBactPrev[[2]])
tmp12_56 <- CompEntPred.pair(ent.Raes12.3e.noBactPrev[[2]],ent.Raes56.3e.noBactPrev[[2]])
# and get cluster correspondences
tmp12_34$Optimal_cluster_correspondences
tmp12_56$Optimal_cluster_correspondences
tmp34_56$Optimal_cluster_correspondences


#3. BACTEROIDES TO PREVOTELLA RATIO 
# draw Bacteroides to Prevotella ratio as well as Bacteroides and Prevotella abundance
setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Check_liter_review/Bateroides_Prevotella')
locations <- c('','\\.IL','\\.TR','\\.RE')
# v1v2
v1v2_table <- read.table('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt', header=T, row.names=1, dec=".", sep="\t")
v1v2_col_sum <- colSums(v1v2_table)
v1v2_table_norm <- v1v2_table
for (i in 1:dim(v1v2_table)[2]){
  v1v2_table_norm[,i] <- v1v2_table[,i]/v1v2_col_sum[i]
}
BactPrev_v1v2b <- v1v2_table_norm[rownames(v1v2_table)%like%'g__Bacteroides',]
BactPrev_v1v2p <- v1v2_table_norm[rownames(v1v2_table)%like%'g__Prevotella',]
BactPrev_v1v2r <- as.numeric(BactPrev_v1v2b)/as.numeric(BactPrev_v1v2p)
BactPrev_v1v2r_log <- log(BactPrev_v1v2r)
names(BactPrev_v1v2r_log) <- names(BactPrev_v1v2p)
BactPrev_v1v2 <- rbind(BactPrev_v1v2b,BactPrev_v1v2p,BactPrev_v1v2r_log)
jpeg('v1v2_Bact_Prev_ratio_density.jpg',width = 1600,height = 1200)
par(mfrow = c(2, 2),mar = c (9,9,7,10),mgp = c(4,1,0),adj = 0.5)
for (i in locations){
  j = i
  if (i==''){j <- '__all loc'}
  plot(density(na.omit(BactPrev_v1v2r_log[names(BactPrev_v1v2r_log)%like%i])),
       main=paste0('Bacteroides/Prevotella rat. ',str_sub(j,start=3)), 
       xlab = 'log(Bacteroides/Prevotella)', cex.main = 5,cex.lab = 4, cex.axis = 2)
  title(sub = paste0('Number of samples: ',length(na.omit(BactPrev_v1v2r_log[names(BactPrev_v1v2r_log)%like%i]))),cex.sub = 3,line = -3,col.sub = 'grey')
  grid(nx=10,lty = 6)
}
dev.off()
# Draw Bacteroides abundance
jpeg('v1v2_Bact_abundance_density.jpg',width = 1600,height = 1200)
par(mfrow = c(2, 2),mar = c (9,9,7,10),mgp = c(4,1,0),adj = 0.5)
for (i in locations){
  j = i
  if (i==''){j <- '__all loc'}
  plot(density(as.numeric(BactPrev_v1v2b[names(BactPrev_v1v2b)%like%i])),
       main=paste0('Bacteroides abundance ',str_sub(j,start=3)),
       xlab = 'Bacteroides abundance',cex.main = 5,cex.lab = 4, cex.axis = 2)
  title(sub = paste0('Number of samples: ',length(na.omit(BactPrev_v1v2b[names(BactPrev_v1v2b)%like%i]))),cex.sub = 3,line = -3,col.sub = 'grey')
  grid(nx=10,lty = 6)
}
dev.off()
# Draw Prevotella abundance
jpeg('v1v2_Prev_abundance_density.jpg',width = 1600,height = 1200)
par(mfrow = c(2, 2),mar = c (5,9,7,10),mgp = c(4,1,0),adj = 1)
for (i in locations){
  j = i
  if (i==''){j <- '__all loc'}
  plot(density(as.numeric(BactPrev_v1v2p[names(BactPrev_v1v2p)%like%i])),
       main=paste0('Prevotella abundance ',str_sub(j,start=3)),
       xlab = 'Prevotella abundance',cex.main = 5,cex.lab = 4, cex.axis = 2)
  title(sub = paste0('Number of samples: ',length(na.omit(BactPrev_v1v2p[names(BactPrev_v1v2p)%like%i]))),cex.sub = 3,line = -13,col.sub = 'grey')
  grid(nx=10,lty = 6)
}
dev.off()


# v3v4
v3v4_table <- read.table('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt', header=T, row.names=1, dec=".", sep="\t")
v3v4_col_sum <- colSums(v3v4_table)
v3v4_table_norm <- v3v4_table
for (i in 1:dim(v3v4_table)[2]){
  v3v4_table_norm[,i] <- v3v4_table[,i]/v3v4_col_sum[i]
}
BactPrev_v3v4b <- v3v4_table_norm[rownames(v3v4_table)%like%'g__Bacteroides',]
BactPrev_v3v4p <- v3v4_table_norm[rownames(v3v4_table)%like%'g__Prevotella',]
BactPrev_v3v4r <- as.numeric(BactPrev_v3v4b)/as.numeric(BactPrev_v3v4p)
BactPrev_v3v4r_log <- log(BactPrev_v3v4r)
names(BactPrev_v3v4r_log) <- names(BactPrev_v3v4p)
BactPrev_v3v4 <- rbind(BactPrev_v3v4b,BactPrev_v3v4p,BactPrev_v3v4r_log)
jpeg('v3v4_Bact_Prev_ratio_density.jpg',width = 1600,height = 1200)
par(mfrow = c(2, 2),mar = c (5,9,7,10))
for (i in locations){
  j = i
  if (i==''){j <- '__all loc'}
  plot(density(na.omit(BactPrev_v3v4r_log[names(BactPrev_v3v4r_log)%like%i])),main=paste0('Bacteroides/Prevotella rat. ',str_sub(j,start=3)),xlab = 'log(Bacteroides/Prevotella)',cex.main = 5,cex.lab = 4, cex.axis = 2)
  grid(nx=10,lty = 6)
}
dev.off()
# Draw Bacteroides abundance
jpeg('v3v4_Bact_abundance_density.jpg',width = 1600,height = 1200)
par(mfrow = c(2, 2),mar = c (5,9,7,10))
for (i in locations){
  j = i
  if (i==''){j <- '__all loc'}
  plot(density(as.numeric(BactPrev_v3v4b[names(BactPrev_v3v4b)%like%i])),main=paste0('Bacteroides abundance ',str_sub(j,start=3)),xlab = 'Bacteroides abundance',cex.main = 5,cex.lab = 4, cex.axis = 2)
  grid(nx=10,lty = 6)
}
dev.off()
# Draw Prevotella abundance
jpeg('v3v4_Prev_abundance_density.jpg',width = 1600,height = 1200)
par(mfrow = c(2, 2),mar = c (5,9,7,10))
for (i in locations){
  j = i
  if (i==''){j <- '__all loc'}
  plot(density(as.numeric(BactPrev_v3v4p[names(BactPrev_v3v4p)%like%i])),main=paste0('Prevotella abundance ',str_sub(j,start=3)),xlab = 'Prevotella abundance',cex.main = 5,cex.lab = 4, cex.axis = 2)
  grid(nx=10,lty = 6)
}
dev.off()

# v5v6
v5v6_table <- read.table('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt', header=T, row.names=1, dec=".", sep="\t")
v5v6_col_sum <- colSums(v5v6_table)
v5v6_table_norm <- v5v6_table
for (i in 1:dim(v5v6_table)[2]){
  v5v6_table_norm[,i] <- v5v6_table[,i]/v5v6_col_sum[i]
}
BactPrev_v5v6b <- v5v6_table_norm[rownames(v5v6_table)%like%'g__Bacteroides',]
BactPrev_v5v6p <- v5v6_table_norm[rownames(v5v6_table)%like%'g__Prevotella',]
BactPrev_v5v6r <- as.numeric(BactPrev_v5v6b)/as.numeric(BactPrev_v5v6p)
BactPrev_v5v6r_log <- log(BactPrev_v5v6r)
names(BactPrev_v5v6r_log) <- names(BactPrev_v5v6p)
BactPrev_v5v6 <- rbind(BactPrev_v5v6b,BactPrev_v5v6p,BactPrev_v5v6r_log)
jpeg('v5v6_Bact_Prev_ratio_density.jpg',width = 1600,height = 1200)
par(mfrow = c(2, 2),mar = c (5,9,7,10))
for (i in locations){
  j = i
  if (i==''){j <- '__all loc'}
  plot(density(na.omit(BactPrev_v5v6r_log[names(BactPrev_v5v6r_log)%like%i])),main=paste0('Bacteroides/Prevotella rat. ',str_sub(j,start=3)),xlab = 'log(Bacteroides/Prevotella)',cex.main = 5,cex.lab = 4, cex.axis = 2)
  grid(nx=10,lty = 6)
}
dev.off()
# Draw Bacteroides abundance
jpeg('v5v6_Bact_abundance_density.jpg',width = 1600,height = 1200)
par(mfrow = c(2, 2),mar = c (5,9,7,10))
for (i in locations){
  j = i
  if (i==''){j <- '__all loc'}
  plot(density(as.numeric(BactPrev_v5v6b[names(BactPrev_v5v6b)%like%i])),main=paste0('Bacteroides abundance ',str_sub(j,start=3)),xlab = 'Bacteroides abundance',cex.main = 5,cex.lab = 4, cex.axis = 2)
  grid(nx=10,lty = 6)
}
dev.off()
# Draw Prevotella abundance
jpeg('v5v6_Prev_abundance_density.jpg',width = 1600,height = 1200)
par(mfrow = c(2, 2),mar = c (5,9,7,10))
for (i in locations){
  j = i
  if (i==''){j <- '__all loc'}
  plot(density(as.numeric(BactPrev_v5v6p[names(BactPrev_v5v6p)%like%i])),main=paste0('Prevotella abundance ',str_sub(j,start=3)),xlab = 'Prevotella abundance',cex.main = 5,cex.lab = 4, cex.axis = 2)
  grid(nx=10,lty = 6)
}
dev.off()

# SUPPLEMENT do some additional calculations based on v1v2
# calculate number of samples with high Bact/Prev ratio
RE <- na.omit(BactPrev_v1v2r_log[names(BactPrev_v1v2r_log)%like%'\\.RE'])
IL <- na.omit(BactPrev_v1v2r_log[names(BactPrev_v1v2r_log)%like%'\\.IL'])
TR <- na.omit(BactPrev_v1v2r_log[names(BactPrev_v1v2r_log)%like%'\\.TR'])
length(IL[IL>4])
length(TR[TR>4])
length(RE[RE>4])

# average Bacteroides prevalence
for (i in locations){
   print(paste0(i,mean(as.numeric(BactPrev_v1v2b[names(BactPrev_v1v2b)%like%i]))))
}

# average Prevotella prevalence
for (i in locations){
  print(paste0(i,mean(as.numeric(BactPrev_v1v2p[names(BactPrev_v1v2p)%like%i]))))
}


#4. COMPARE IL, TR and RE in alpha and beta diversity
file_input = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt'
file_input_v12 = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt'
# let us draw a graph of alpha diversity vs threshold
threshold <- seq(0,7,by=0.25)
species.IL <- threshold
for (i in 1:length(threshold)){
  species.IL[i] <- length(getSpec(file_input,loc = 'IL',percent_threshold = exp(-log(10)*threshold[i])))
}
species.IL.plot <- cbind(threshold,species.IL)
# TR
species.TR <- threshold
for (i in 1:length(threshold)){
  species.TR[i] <- length(getSpec(file_input,loc = 'TR',percent_threshold = exp(-log(10)*threshold[i])))
}
species.TR.plot <- cbind(threshold,species.TR)
plot(species.TR.plot,type = 'l',xlab='filtration threshold, -log(x) ,%',ylab='Transverse colon')
grid(nx=10,lty = 6)
# RE
species.RE <- threshold
for (i in 1:length(threshold)){
  species.RE[i] <- length(getSpec(file_input,loc = 'RE',percent_threshold = exp(-log(10)*threshold[i])))
}
species.RE.plot <- cbind(threshold,species.RE)
# make plots of species vs % threshold
setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Check_liter_review/Location_compar')
jpeg('IL_TR_RE_microbiota_richness.jpg',width = 800, height = 600)
par(mar = c(5,6,5,3))
plot(species.IL.plot,type = 'l',xlab='filtration threshold, -log scale,%', ylab = 'Number of genera detected', ylim = c(0,850), main ='Intestine location microbiota richness',col='red',cex.main = 3, cex.lab = 3, cex.axis = 2)
grid(nx=10,lty = 6)
lines(species.TR.plot, col='blue')
lines(species.RE.plot, col='black')
dev.off()
# this is atavistic a bit - better to use species.IL/TR/RE
# v5v6
species.v56.IL.01 <- getSpec(file_input,loc = 'IL',percent_threshold = 0.1)
species.v56.TR.01 <- getSpec(file_input,loc = 'TR',percent_threshold = 0.1)
species.v56.RE.01 <- getSpec(file_input,loc = 'RE',percent_threshold = 0.1)
species.v56.IL.001 <- getSpec(file_input,loc = 'IL',percent_threshold = 0.01)
species.v56.TR.001 <- getSpec(file_input,loc = 'TR',percent_threshold = 0.01)
species.v56.RE.001 <- getSpec(file_input,loc = 'RE',percent_threshold = 0.01)
species.v56.IL.002 <- getSpec(file_input,loc = 'IL',percent_threshold = 0.02)
species.v56.TR.002 <- getSpec(file_input,loc = 'TR',percent_threshold = 0.02)
species.v56.RE.002 <- getSpec(file_input,loc = 'RE',percent_threshold = 0.02)
species.v56.IL.0005 <- getSpec(file_input,loc = 'IL',percent_threshold = 0.005)
species.v56.TR.0005 <- getSpec(file_input,loc = 'TR',percent_threshold = 0.005)
species.v56.RE.0005 <- getSpec(file_input,loc = 'RE',percent_threshold = 0.005)
species.v56.IL.0001 <- getSpec(file_input,loc = 'IL',percent_threshold = 0.001)
species.v56.TR.0001 <- getSpec(file_input,loc = 'TR',percent_threshold = 0.001)
species.v56.RE.0001 <- getSpec(file_input,loc = 'RE',percent_threshold = 0.001)
# v1v2 atavistic
species.v12.IL.01 <- getSpec(file_input_v12,loc = 'IL',percent_threshold = 0.1)
species.v12.TR.01 <- getSpec(file_input_v12,loc = 'TR',percent_threshold = 0.1)
species.v12.RE.01 <- getSpec(file_input_v12,loc = 'RE',percent_threshold = 0.1)
species.v12.IL.0001 <- getSpec(file_input_v12,loc = 'IL',percent_threshold = 0.001)
species.v12.TR.0001 <- getSpec(file_input_v12,loc = 'TR',percent_threshold = 0.001)
species.v12.RE.0001 <- getSpec(file_input_v12,loc = 'RE',percent_threshold = 0.001)
# output location comparison as a table (use atavisms)
alpha_diversity <- matrix(0,5,3)
colnames(alpha_diversity) <- c('IL','TR','RE')
rownames(alpha_diversity) <- c('0.001%','0.005%','0.01%','0.02%','0.1%')
alpha_diversity['0.001%',] <- c(length(species.v56.IL.0001),length(species.v56.TR.0001),length(species.v56.RE.0001))
alpha_diversity['0.005%',] <- c(length(species.v56.IL.0005),length(species.v56.TR.0005),length(species.v56.RE.0005))
alpha_diversity['0.01%',] <- c(length(species.v56.IL.001),length(species.v56.TR.001),length(species.v56.RE.001))
alpha_diversity['0.02%',] <- c(length(species.v56.IL.002),length(species.v56.TR.002),length(species.v56.RE.002))
alpha_diversity['0.1%',] <- c(length(species.v56.IL.01),length(species.v56.TR.01),length(species.v56.RE.01))
alpha_diversity

# now see species different in different locations
# IL vs TR
IL.unique.comp.to.TR <- setdiff(species.v56.IL.002,species.v56.TR.0001)
TR.unique.comp.to.IL <- setdiff(species.v56.TR.002,species.v56.IL.0001)
# IL vs TR
IL.unique.comp.to.RE <- setdiff(species.v56.IL.01,species.v56.RE.0001)
RE.unique.comp.to.IL <- setdiff(species.v56.RE.01,species.v56.IL.0001)
# TR vs RE
TR.unique.comp.to.RE <- setdiff(species.v56.TR.01,species.v56.RE.0001)
RE.unique.comp.to.TR <- setdiff(species.v56.RE.01,species.v56.TR.0001)
IL.unique.comp.to.TR
IL.unique.comp.to.RE
TR.unique.comp.to.IL
TR.unique.comp.to.RE
RE.unique.comp.to.IL
RE.unique.comp.to.TR
# Now check for v1v2
# IL vs TR
IL.unique.comp.to.TR_v12 <- setdiff(species.v12.IL.01,species.v12.TR.0001)
TR.unique.comp.to.IL_v12 <- setdiff(species.v12.TR.01,species.v12.IL.0001)
# IL vs TR
IL.unique.comp.to.RE_v12 <- setdiff(species.v12.IL.01,species.v12.RE.0001)
RE.unique.comp.to.IL_v12 <- setdiff(species.v12.RE.01,species.v12.IL.0001)
# TR vs RE
TR.unique.comp.to.RE_v12 <- setdiff(species.v12.TR.01,species.v12.RE.0001)
RE.unique.comp.to.TR_v12 <- setdiff(species.v12.RE.01,species.v12.TR.0001)
IL.unique.comp.to.TR_v12
IL.unique.comp.to.RE_v12
TR.unique.comp.to.RE_v12
TR.unique.comp.to.IL_v12
RE.unique.comp.to.TR_v12
RE.unique.comp.to.IL_v12

