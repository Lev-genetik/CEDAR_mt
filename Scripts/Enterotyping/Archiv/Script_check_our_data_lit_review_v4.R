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
par(mfrow = c(2, 2),mar = c (10,9,7,10),mgp = c(5,1,0),adj = 0.5)
for (i in locations){
  j = i
  if (i==''){j <- '__all loc'}
  plot(density(na.omit(BactPrev_v1v2r_log[names(BactPrev_v1v2r_log)%like%i])),
       main=paste0('Bacteroides/Prevotella rat. ',str_sub(j,start=3)), 
       xlab = expression('log'[10]*'(Bacteroides/Prevotella)'), cex.main = 5,cex.lab = 4, cex.axis = 2)
  title(sub = paste0('Number of samples: ',length(na.omit(BactPrev_v1v2r_log[names(BactPrev_v1v2r_log)%like%i]))),cex.sub = 3,line = 8,col.sub = 'grey')
  grid(nx=10,lty = 6)
}
dev.off()
# Draw Bacteroides abundance
jpeg('v1v2_Bact_abundance_density.jpg',width = 1600,height = 1200)
par(mfrow = c(2, 2),mar = c (10,9,7,10),mgp = c(4,1,0),adj = 0.5)
for (i in locations){
  j = i
  if (i==''){j <- '__all loc'}
  plot(density(as.numeric(BactPrev_v1v2b[names(BactPrev_v1v2b)%like%i])),
       main=paste0('Bacteroides abundance ',str_sub(j,start=3)),
       xlab = 'Bacteroides fraction',cex.main = 5,cex.lab = 4, cex.axis = 2)
  title(sub = paste0('Number of samples: ',length(na.omit(BactPrev_v1v2b[names(BactPrev_v1v2b)%like%i]))),cex.sub = 3,line = 8,col.sub = 'grey')
  grid(nx=10,lty = 6)
}
dev.off()
# Draw Prevotella abundance
jpeg('v1v2_Prev_abundance_density.jpg',width = 1600,height = 1200)
par(mfrow = c(2, 2),mar = c (10,9,7,10),mgp = c(4,1,0))
for (i in locations){
  j = i
  if (i==''){j <- '__all loc'}
  plot(density(as.numeric(BactPrev_v1v2p[names(BactPrev_v1v2p)%like%i])),
       main=paste0('Prevotella abundance ',str_sub(j,start=3)),
       xlab = 'Prevotella fraction',cex.main = 5,cex.lab = 4, cex.axis = 2)
  title(sub = paste0('Number of samples: ',length(na.omit(BactPrev_v1v2p[names(BactPrev_v1v2p)%like%i]))),cex.sub = 3,line = 8,col.sub = 'grey')
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
par(mfrow = c(2, 2),mar = c (10,9,7,10),mgp = c(5,1,0),adj = 0.5)
for (i in locations){
  j = i
  if (i==''){j <- '__all loc'}
  plot(density(na.omit(BactPrev_v3v4r_log[names(BactPrev_v3v4r_log)%like%i])),
       main=paste0('Bacteroides/Prevotella rat. ',str_sub(j,start=3)), 
       xlab = expression('log'[10]*'(Bacteroides/Prevotella)'), cex.main = 5,cex.lab = 4, cex.axis = 2)
  title(sub = paste0('Number of samples: ',length(na.omit(BactPrev_v3v4r_log[names(BactPrev_v3v4r_log)%like%i]))),cex.sub = 3,line = 8,col.sub = 'grey')
  grid(nx=10,lty = 6)
}
dev.off()
# Draw Bacteroides abundance
jpeg('v3v4_Bact_abundance_density.jpg',width = 1600,height = 1200)
par(mfrow = c(2, 2),mar = c (10,9,7,10),mgp = c(4,1,0),adj = 0.5)
for (i in locations){
  j = i
  if (i==''){j <- '__all loc'}
  plot(density(as.numeric(BactPrev_v3v4b[names(BactPrev_v3v4b)%like%i])),
       main=paste0('Bacteroides abundance ',str_sub(j,start=3)),
       xlab = 'Bacteroides fraction',cex.main = 5,cex.lab = 4, cex.axis = 2)
  title(sub = paste0('Number of samples: ',length(na.omit(BactPrev_v3v4b[names(BactPrev_v3v4b)%like%i]))),cex.sub = 3,line = 8,col.sub = 'grey')
  grid(nx=10,lty = 6)
}
dev.off()
# Draw Prevotella abundance
jpeg('v3v4_Prev_abundance_density.jpg',width = 1600,height = 1200)
par(mfrow = c(2, 2),mar = c (10,9,7,10),mgp = c(4,1,0))
for (i in locations){
  j = i
  if (i==''){j <- '__all loc'}
  plot(density(as.numeric(BactPrev_v3v4p[names(BactPrev_v3v4p)%like%i])),
       main=paste0('Prevotella abundance ',str_sub(j,start=3)),
       xlab = 'Prevotella fraction',cex.main = 5,cex.lab = 4, cex.axis = 2)
  title(sub = paste0('Number of samples: ',length(na.omit(BactPrev_v3v4p[names(BactPrev_v3v4p)%like%i]))),cex.sub = 3,line = 8,col.sub = 'grey')
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
par(mfrow = c(2, 2),mar = c (10,9,7,10),mgp = c(5,1,0),adj = 0.5)
for (i in locations){
  j = i
  if (i==''){j <- '__all loc'}
  plot(density(na.omit(BactPrev_v5v6r_log[names(BactPrev_v5v6r_log)%like%i])),
       main=paste0('Bacteroides/Prevotella rat. ',str_sub(j,start=3)), 
       xlab = expression('log'[10]*'(Bacteroides/Prevotella)'), cex.main = 5,cex.lab = 4, cex.axis = 2)
  title(sub = paste0('Number of samples: ',length(na.omit(BactPrev_v5v6r_log[names(BactPrev_v5v6r_log)%like%i]))),cex.sub = 3,line = 8,col.sub = 'grey')
  grid(nx=10,lty = 6)
}
dev.off()
# Draw Bacteroides abundance
jpeg('v5v6_Bact_abundance_density.jpg',width = 1600,height = 1200)
par(mfrow = c(2, 2),mar = c (10,9,7,10),mgp = c(4,1,0),adj = 0.5)
for (i in locations){
  j = i
  if (i==''){j <- '__all loc'}
  plot(density(as.numeric(BactPrev_v5v6b[names(BactPrev_v5v6b)%like%i])),
       main=paste0('Bacteroides abundance ',str_sub(j,start=3)),
       xlab = 'Bacteroides fraction',cex.main = 5,cex.lab = 4, cex.axis = 2)
  title(sub = paste0('Number of samples: ',length(na.omit(BactPrev_v5v6b[names(BactPrev_v5v6b)%like%i]))),cex.sub = 3,line = 8,col.sub = 'grey')
  grid(nx=10,lty = 6)
}
dev.off()
# Draw Prevotella abundance
jpeg('v5v6_Prev_abundance_density.jpg',width = 1600,height = 1200)
par(mfrow = c(2, 2),mar = c (10,9,7,10),mgp = c(4,1,0))
for (i in locations){
  j = i
  if (i==''){j <- '__all loc'}
  plot(density(as.numeric(BactPrev_v5v6p[names(BactPrev_v5v6p)%like%i])),
       main=paste0('Prevotella abundance ',str_sub(j,start=3)),
       xlab = 'Prevotella fraction',cex.main = 5,cex.lab = 4, cex.axis = 2)
  title(sub = paste0('Number of samples: ',length(na.omit(BactPrev_v5v6p[names(BactPrev_v5v6p)%like%i]))),cex.sub = 3,line = 8,col.sub = 'grey')
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

# 4a. The total richness
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Check_data_lit_review_v3.R')
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

# V1V2
file_input = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt'
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
jpeg('v1v2_IL_TR_RE_microbiota_richness.jpg',width = 800, height = 600)
par(mar = c(5,6,5,3))
plot(species.IL.plot,type = 'l',xlab='filtration threshold, -log scale,%', ylab = 'Number of genera detected', ylim = c(0,850), main ='Intestine location microbiota richness',col='red',cex.main = 3, cex.lab = 3, cex.axis = 2)
grid(nx=10,lty = 6)
lines(species.TR.plot, col='blue')
lines(species.RE.plot, col='black')
dev.off()

# V3V4
file_input = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt'
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
jpeg('v3v4_IL_TR_RE_microbiota_richness.jpg',width = 800, height = 600)
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


# 4b. Shannon index and Hill number of order 1 comparison IL vs TR vs RE
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Check_data_lit_review_v3.R')
# V5V6
# get mean Shannon entropy and Hill numbers of order 1 per sample
threshold <- seq(0,7,by=0.25)
# IL
shannon.v56.IL <- threshold
Hill1.v56.IL <- threshold
for (i in 1:length(threshold)){
  enter.shannon <- entrShannon(input_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt',
                                    locat = '.IL', percent_threshold = exp(-log(10)*threshold[i]))
  shannon.v56.IL[i] <- mean(enter.shannon)
  Hill1.v56.IL[i] <- exp(mean(enter.shannon))
}
shannon.v56.IL.plot <- cbind(threshold,shannon.v56.IL)
Hill1.v56.IL.plot <- cbind(threshold,Hill1.v56.IL)
# TR
shannon.v56.TR <- threshold
Hill1.v56.TR <- threshold
for (i in 1:length(threshold)){
  enter.shannon <- entrShannon(input_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt',
                               locat = '.TR', percent_threshold = exp(-log(10)*threshold[i]))
  shannon.v56.TR[i] <- mean(enter.shannon)
  Hill1.v56.TR[i] <- exp(mean(enter.shannon))
}
shannon.v56.TR.plot <- cbind(threshold,shannon.v56.TR)
Hill1.v56.TR.plot <- cbind(threshold,Hill1.v56.TR)
# RE
shannon.v56.RE <- threshold
Hill1.v56.RE <- threshold
for (i in 1:length(threshold)){
  enter.shannon <- entrShannon(input_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt',
                               locat = '.RE', percent_threshold = exp(-log(10)*threshold[i]))
  shannon.v56.RE[i] <- mean(enter.shannon)
  Hill1.v56.RE[i] <- exp(mean(enter.shannon))
}
shannon.v56.RE.plot <- cbind(threshold,shannon.v56.RE)
Hill1.v56.RE.plot <- cbind(threshold,Hill1.v56.RE)

# V1V2
# get mean Shannon entropy and Hill numbers of order 1 per sample
threshold <- seq(0,7,by=0.25)
# IL
shannon.v12.IL <- threshold
Hill1.v12.IL <- threshold
for (i in 1:length(threshold)){
  enter.shannon <- entrShannon(input_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt',
                               locat = '.IL', percent_threshold = exp(-log(10)*threshold[i]))
  shannon.v12.IL[i] <- mean(enter.shannon)
  Hill1.v12.IL[i] <- exp(mean(enter.shannon))
}
shannon.v12.IL.plot <- cbind(threshold,shannon.v12.IL)
Hill1.v12.IL.plot <- cbind(threshold,Hill1.v12.IL)
# TR
shannon.v12.TR <- threshold
Hill1.v12.TR <- threshold
for (i in 1:length(threshold)){
  enter.shannon <- entrShannon(input_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt',
                               locat = '.TR', percent_threshold = exp(-log(10)*threshold[i]))
  shannon.v12.TR[i] <- mean(enter.shannon)
  Hill1.v12.TR[i] <- exp(mean(enter.shannon))
}
shannon.v12.TR.plot <- cbind(threshold,shannon.v12.TR)
Hill1.v12.TR.plot <- cbind(threshold,Hill1.v12.TR)
# RE
shannon.v12.RE <- threshold
Hill1.v12.RE <- threshold
for (i in 1:length(threshold)){
  enter.shannon <- entrShannon(input_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt',
                               locat = '.RE', percent_threshold = exp(-log(10)*threshold[i]))
  shannon.v12.RE[i] <- mean(enter.shannon)
  Hill1.v12.RE[i] <- exp(mean(enter.shannon))
}
shannon.v12.RE.plot <- cbind(threshold,shannon.v12.RE)
Hill1.v12.RE.plot <- cbind(threshold,Hill1.v12.RE)

# V3V4
# get mean Shannon entropy and Hill numbers of order 1 per sample
threshold <- seq(0,7,by=0.25)
# IL
shannon.v34.IL <- threshold
Hill1.v34.IL <- threshold
for (i in 1:length(threshold)){
  enter.shannon <- entrShannon(input_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt',
                               locat = '.IL', percent_threshold = exp(-log(10)*threshold[i]))
  shannon.v34.IL[i] <- mean(enter.shannon)
  Hill1.v34.IL[i] <- exp(mean(enter.shannon))
}
shannon.v34.IL.plot <- cbind(threshold,shannon.v34.IL)
Hill1.v34.IL.plot <- cbind(threshold,Hill1.v34.IL)
# TR
shannon.v34.TR <- threshold
Hill1.v34.TR <- threshold
for (i in 1:length(threshold)){
  enter.shannon <- entrShannon(input_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt',
                               locat = '.TR', percent_threshold = exp(-log(10)*threshold[i]))
  shannon.v34.TR[i] <- mean(enter.shannon)
  Hill1.v34.TR[i] <- exp(mean(enter.shannon))
}
shannon.v34.TR.plot <- cbind(threshold,shannon.v34.TR)
Hill1.v34.TR.plot <- cbind(threshold,Hill1.v34.TR)
# RE
shannon.v34.RE <- threshold
Hill1.v34.RE <- threshold
for (i in 1:length(threshold)){
  enter.shannon <- entrShannon(input_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt',
                               locat = '.RE', percent_threshold = exp(-log(10)*threshold[i]))
  shannon.v34.RE[i] <- mean(enter.shannon)
  Hill1.v34.RE[i] <- exp(mean(enter.shannon))
}
shannon.v34.RE.plot <- cbind(threshold,shannon.v34.RE)
Hill1.v34.RE.plot <- cbind(threshold,Hill1.v34.RE)

# V5V6
# make plot of Shannon entropy vs % threshold
setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Check_liter_review/Location_compar')
jpeg('v56_IL_TR_RE_microbiota_Shannon_entropy.jpg',width = 800, height = 600)
par(mar = c(5,6,5,3))
plot(shannon.v56.IL.plot,type = 'l',xlab='filtration threshold, -log scale,%', ylab = 'Shannon entropy', ylim = c(1.5,2.5), main ='Intestine location microbiota alpha diversity',col='red',cex.main = 3, cex.lab = 3, cex.axis = 2)
grid(nx=10,lty = 6)
lines(shannon.v56.TR.plot, col='blue')
lines(shannon.v56.RE.plot, col='black')
dev.off()
# make plot of Hill numbers of order 1 vs % threshold
setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Check_liter_review/Location_compar')
jpeg('v56_IL_TR_RE_microbiota_Hill_order1.jpg',width = 800, height = 600)
par(mar = c(5,6,5,3))
plot(Hill1.v56.IL.plot,type = 'l',xlab='filtration threshold, -log scale,%', ylab = expression(''^1*'D'), ylim = c(5,10), main ='Intestine location microbiota alpha diversity',col='red',cex.main = 3, cex.lab = 3, cex.axis = 2)
grid(nx=10,lty = 6)
lines(Hill1.v56.TR.plot, col='blue')
lines(Hill1.v56.RE.plot, col='black')
dev.off()

# V3V4
# make plot of Shannon entropy vs % threshold
setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Check_liter_review/Location_compar')
jpeg('v34_IL_TR_RE_microbiota_Shannon_entropy.jpg',width = 800, height = 600)
par(mar = c(5,6,5,3))
plot(shannon.v34.IL.plot,type = 'l',xlab='filtration threshold, -log scale,%', ylab = 'Shannon entropy', ylim = c(1.5,2.5), main ='Intestine location microbiota alpha diversity',col='red',cex.main = 3, cex.lab = 3, cex.axis = 2)
grid(nx=10,lty = 6)
lines(shannon.v34.TR.plot, col='blue')
lines(shannon.v34.RE.plot, col='black')
dev.off()
# make plot of Hill numbers of order 1 vs % threshold
setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Check_liter_review/Location_compar')
jpeg('v34_IL_TR_RE_microbiota_Hill_order1.jpg',width = 800, height = 600)
par(mar = c(5,6,5,3))
plot(Hill1.v34.IL.plot,type = 'l',xlab='filtration threshold, -log scale,%', ylab = expression(''^1*'D'), ylim = c(5,10), main ='Intestine location microbiota alpha diversity',col='red',cex.main = 3, cex.lab = 3, cex.axis = 2)
grid(nx=10,lty = 6)
lines(Hill1.v34.TR.plot, col='blue')
lines(Hill1.v34.RE.plot, col='black')
dev.off()

# V1V2
# make plot of Shannon entropy vs % threshold
setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Check_liter_review/Location_compar')
jpeg('v12_IL_TR_RE_microbiota_Shannon_entropy.jpg',width = 800, height = 600)
par(mar = c(5,6,5,3))
plot(shannon.v12.IL.plot,type = 'l',xlab='filtration threshold, -log scale,%', ylab = 'Shannon entropy', ylim = c(1.5,2.5), main ='Intestine location microbiota alpha diversity',col='red',cex.main = 3, cex.lab = 3, cex.axis = 2)
grid(nx=10,lty = 6)
lines(shannon.v12.TR.plot, col='blue')
lines(shannon.v12.RE.plot, col='black')
dev.off()
# make plot of Hill numbers of order 1 vs % threshold
setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Check_liter_review/Location_compar')
jpeg('v12_IL_TR_RE_microbiota_Hill_order1.jpg',width = 800, height = 600)
par(mar = c(5,6,5,3))
plot(Hill1.v12.IL.plot,type = 'l',xlab='filtration threshold, -log scale,%', ylab = expression(''^1*'D'), ylim = c(5,10), main ='Intestine location microbiota alpha diversity',col='red',cex.main = 3, cex.lab = 3, cex.axis = 2)
grid(nx=10,lty = 6)
lines(Hill1.v12.TR.plot, col='blue')
lines(Hill1.v12.RE.plot, col='black')
dev.off()


# 5. Compare taxa in different locations (and amplicons)
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Check_data_lit_review_v3.R')
setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Check_liter_review/Location_compar')
jpeg('v1v2_IL_TR_RE.jpg',width = 3800, height = 3800)
compTaxa(file1='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt',
         file2='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt',
         file3='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt',
         loc1='.IL',loc2 = '.TR',loc3 = '.RE')
dev.off()
jpeg('v3v4_IL_TR_RE.jpg',width = 3800, height = 3800)
compTaxa(file1='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt',
         file2='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt',
         file3='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt',
         loc1='.IL',loc2 = '.TR',loc3 = '.RE')
dev.off()
jpeg('v5v6_IL_TR_RE.jpg',width = 3800, height = 3800)
compTaxa(file1='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt',
         file2='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt',
         file3='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt',
         loc1='.IL',loc2 = '.TR',loc3 = '.RE')
dev.off()
jpeg('v1v2_v3v4_v5v6_IL',width = 3800, height = 3800)
compTaxa(file1='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt',
         file2='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt',
         file3='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt',
         loc1='.IL',loc2 = '.IL',loc3 = '.IL')
dev.off()
jpeg('v1v2_v3v4_v5v6_TR',width = 3800, height = 3800)
compTaxa(file1='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt',
         file2='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt',
         file3='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt',
         loc1='.TR',loc2 = '.TR',loc3 = '.TR')
dev.off()
jpeg('v1v2_v3v4_v5v6_RE',width = 3800, height = 3800)
compTaxa(file1='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt',
         file2='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt',
         file3='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt',
         loc1='.RE',loc2 = '.RE',loc3 = '.RE')
dev.off()


# Bacteria fraction agreement across amplicons
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Check_data_lit_review_v3.R')
setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Check_liter_review/Bacteria_fract_by_loc')
jpeg('v12_v34_IL_IPC293.jpg',width = 800,height = 800)
plotBactFract(file1 = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt',
              file2 = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt',
              sample_name = 'IPC293', locat='.IL', percent_threshold=0.01, remove.outliers=0.01)
dev.off()
jpeg('v12_v56_IL_IPC293.jpg',width = 800,height = 800)
plotBactFract(file1 = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt',
              file2 = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt',
              sample_name = 'IPC293', locat='.IL', percent_threshold=0.01, remove.outliers=0.01)
dev.off()  
jpeg('v34_v56_IL_IPC293.jpg',width = 800,height = 800)
plotBactFract(file1 = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt',
              file2 = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt',
              sample_name = 'IPC293', locat='.IL', percent_threshold=0.01, remove.outliers=0.01)
dev.off()



#6. Explore number of reads per sample
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Check_data_lit_review_v3.R')
setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Check_liter_review/Sample_coverage')
locations <- c('all','.IL','.TR','.RE')
# V1V2
jpeg('v1v2_samples_microbiota_coverage_stats_cumul.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2),mar = c(8,10,7,5), mgp=c(6,2,0))
for (i in locations){
exploreReadsPerSample(input_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt',
                       locat = i, percent_threshold = 0, remove.outliers = 0)
}
dev.off()
# V3V4
jpeg('v3v4_samples_microbiota_coverage_stats_cumul.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2),mar = c(8,10,7,5), mgp=c(6,2,0))
for (i in locations){
  exploreReadsPerSample(input_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt',
                        locat = i, percent_threshold = 0, remove.outliers = 0)
}
dev.off()
# V5V6
jpeg('v5v6_samples_microbiota_coverage_stats_cumul.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2),mar = c(8,10,7,5), mgp=c(6,2,0))
for (i in locations){
  exploreReadsPerSample(input_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt',
                        locat = i, percent_threshold = 0, remove.outliers = 0)
}
dev.off()

# non-cumulative graphs
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Check_data_lit_review_v3.R')
setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Check_liter_review/Sample_coverage')
locations <- c('all','.IL','.TR','.RE')
# V1V2
jpeg('v1v2_samples_microbiota_coverage_stats.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2),mar = c(8,10,7,5), mgp=c(6,2,0))
for (i in locations){
  exploreReadsPerSample(input_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt',
                        locat = i, percent_threshold = 0, remove.outliers = 0,cumulative=FALSE)
}
dev.off()
# V3V4
jpeg('v3v4_samples_microbiota_coverage_stats.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2),mar = c(8,10,7,5), mgp=c(6,2,0))
for (i in locations){
  exploreReadsPerSample(input_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt',
                        locat = i, percent_threshold = 0, remove.outliers = 0,cumulative=FALSE)
}
dev.off()
# V5V6
jpeg('v5v6_samples_microbiota_coverage_stats.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2),mar = c(8,10,7,5), mgp=c(6,2,0))
for (i in locations){
  exploreReadsPerSample(input_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt',
                        locat = i, percent_threshold = 0, remove.outliers = 0,cumulative=FALSE)
}
dev.off()



#7. Correlate number of reads per sample with enterotype stability
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Check_data_lit_review_v3.R')
library(data.table)
setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Check_liter_review/Enterotype_coverage')
locations234 <- c('.IL','.TR','.RE')
stable_enterotype_file <- c('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/Old-the-same/Stable_enterotypes_IL.csv',
                            '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/Old-the-same/Stable_enterotypes_TR.csv',
                            '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/Old-the-same/Stable_enterotypes_RE.csv')
# V1V2
v1v2_cov_vs_ent_stab_IL <- rep('',times=4)
v1v2_cov_vs_ent_stab_TR <- rep('',times=4)
v1v2_cov_vs_ent_stab_RE <- rep('',times=4)
v1v2_cov_vs_ent_stab <- rbind(v1v2_cov_vs_ent_stab_IL,v1v2_cov_vs_ent_stab_TR,v1v2_cov_vs_ent_stab_RE)
for(i in 1:dim(v1v2_cov_vs_ent_stab)[1]){
  v1v2_cov_vs_ent_stab[i,] <- corrCoverEnterotype(qiime_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt',
                        locat = locations234[i], stable_enterotype_file = stable_enterotype_file[i])
}
colnames(v1v2_cov_vs_ent_stab) <- c('Mann-Whitney_res','t.test_res', 'mean_cov_stable', 'mean_cov_unstable')
 write.table(v1v2_cov_vs_ent_stab,file = 'v1v2_coverage_vs_ent_stabil.csv',append = F,quote = F,sep = '\t',eol='\n') 
 # V3V4
 v3v4_cov_vs_ent_stab_IL <- rep('',times=4)
 v3v4_cov_vs_ent_stab_TR <- rep('',times=4)
 v3v4_cov_vs_ent_stab_RE <- rep('',times=4)
 v3v4_cov_vs_ent_stab <- rbind(v3v4_cov_vs_ent_stab_IL,v3v4_cov_vs_ent_stab_TR,v3v4_cov_vs_ent_stab_RE)
 for(i in 1:dim(v3v4_cov_vs_ent_stab)[1]){
  v3v4_cov_vs_ent_stab[i,] <- corrCoverEnterotype(qiime_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt',
                                                   locat = locations234[i], stable_enterotype_file = stable_enterotype_file[i])
 }
 colnames(v3v4_cov_vs_ent_stab) <- c('Mann-Whitney_res','t.test_res', 'mean_cov_stable', 'mean_cov_unstable')
 write.table(v3v4_cov_vs_ent_stab,file = 'v3v4_coverage_vs_ent_stabil.csv',append = F,quote = F,sep = '\t',eol='\n') 
 #v5v6
 v5v6_cov_vs_ent_stab_IL <- rep('',times=4)
 v5v6_cov_vs_ent_stab_TR <- rep('',times=4)
 v5v6_cov_vs_ent_stab_RE <- rep('',times=4)
 v5v6_cov_vs_ent_stab <- rbind(v5v6_cov_vs_ent_stab_IL,v5v6_cov_vs_ent_stab_TR,v5v6_cov_vs_ent_stab_RE)
 for(i in 1:dim(v5v6_cov_vs_ent_stab)[1]){
   v5v6_cov_vs_ent_stab[i,] <- corrCoverEnterotype(qiime_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt',
                                                   locat = locations234[i], stable_enterotype_file = stable_enterotype_file[i])
 }
 colnames(v5v6_cov_vs_ent_stab) <- c('Mann-Whitney_res','t.test_res', 'mean_cov_stable', 'mean_cov_unstable')
 write.table(v5v6_cov_vs_ent_stab,file = 'v5v6_coverage_vs_ent_stabil.csv',append = F,quote = F,sep = '\t',eol='\n') 
 
  
 #7b. Plot enterotype stability vs number of reads per sample 
 source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Check_data_lit_review_v3.R')
 setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Check_liter_review/Enterotype_coverage/Graphs')
 locations <- c('all','.IL','.TR','.RE')
 stable_enterotype_file <- c('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype/Stable_enterotypes.csv',
                             '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/Old-the-same/Stable_enterotypes_IL.csv',
                             '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/Old-the-same/Stable_enterotypes_TR.csv',
                             '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/Old-the-same/Stable_enterotypes_RE.csv')
 # V1V2
 jpeg('v1v2_enterot_stabil_vs_sample_coverage.jpg',width = 800,height = 800)
 par(mfrow=c(2,2), mar=c(7,7,7,4),mgp=c(4,2,0))
 for(i in 1:length(locations)){
 plotCoverEnterotype(qiime_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt', 
                     locat = locations[i], stable_enterotype_file = stable_enterotype_file[i])
 }
 dev.off()
 # V3V4
 jpeg('v3v4_enterot_stabil_vs_sample_coverage.jpg',width = 800,height = 800)
 par(mfrow=c(2,2), mar=c(7,7,7,4),mgp=c(4,2,0))
 for(i in 1:length(locations)){
   plotCoverEnterotype(qiime_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt', 
                       locat = locations[i], stable_enterotype_file = stable_enterotype_file[i])
 }
 dev.off()
 # V5V6
 jpeg('v5v6_enterot_stabil_vs_sample_coverage.jpg',width = 800,height = 800)
 par(mfrow=c(2,2), mar=c(7,7,7,4),mgp=c(4,2,0))
 for(i in 1:length(locations)){
   plotCoverEnterotype(qiime_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt', 
                       locat = locations[i], stable_enterotype_file = stable_enterotype_file[i])
 }
 dev.off()
 
 #8. How different are abilities of amplicons to catch taxa
 
 