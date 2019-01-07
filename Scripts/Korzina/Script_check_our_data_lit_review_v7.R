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
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Check_data_lit_review_v4.R')
setwd(paste0(home.dir,'/Results/Locations_compare/Genera_abundance_plot'))
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

# Make family-level boxplots
setwd(paste0(home.dir,'/Results/Locations_compare/Abundance_boxplots'))
jpeg('v3v4_location_compare_depth_13k_outl_rem',width = 3800, height = 3800)
compTaxa(file1='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt',
         file2='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt',
         file3='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt',
         loc1='.IL',loc2 = '.TR',loc3 = '.RE',sample.cov.filter=13000,draw.boxplots=T,draw.stackplots=F)
dev.off()


#6. Bacteria fraction agreement across amplicons
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



#7a. Explore number of reads per sample (total befor qual filtr and microbiota-assigned)
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Check_data_lit_review_v3.R')
setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Check_liter_review/Sample_coverage')
locations <- c('all','.IL','.TR','.RE')
# V1V2
jpeg('v1v2_samples_microbiota_coverage_stats_cumul.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2),mar = c(8,10,7,5), mgp=c(6,2,0))
for (i in locations){
exploreCoverage(input_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt',
                       locat = i, percent_threshold = 0, remove.outliers = 0)
}
dev.off()
# V3V4
jpeg('v3v4_samples_microbiota_coverage_stats_cumul.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2),mar = c(8,10,7,5), mgp=c(6,2,0))
for (i in locations){
  exploreCoverage(input_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt',
                        locat = i, percent_threshold = 0, remove.outliers = 0)
}
dev.off()
# V5V6
jpeg('v5v6_samples_microbiota_coverage_stats_cumul.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2),mar = c(8,10,7,5), mgp=c(6,2,0))
for (i in locations){
  exploreCoverage(input_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt',
                        locat = i, percent_threshold = 0, remove.outliers = 0)
}
dev.off()

#7b. non-cumulative graphs
# First for coverage window +-10k reads (e.g. 0 to 19,999) (usredneniye=20) - here we show both micobiota reads
# and reads before quality filtration and assignment to microbiota/human seq
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Check_data_lit_review_v3.R')
setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Check_liter_review/Sample_coverage')
locations <- c('all','.IL','.TR','.RE')
# V1V2
jpeg('v1v2_samples_microbiota_red_and_total_coverage_purple_stats_10k.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2),mar = c(8,10,7,5), mgp=c(6,2,0))
for (i in locations){
  exploreCoverage(input_type='QIIME', input_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt',
                        locat = i, percent_threshold = 0, remove.outliers = 0,cumulative=FALSE,usredneniye=20,y_axis_length_factor = 1.3,
                  color = 'red')
  exploreCoverage(input_type='total.coverage', input_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/stats_v12.txt',
                  locat = i, percent_threshold = 0, remove.outliers = 0,cumulative=FALSE,usredneniye=20,
                  to.append = T,color = 'purple')
  
}
dev.off()
# V3V4
jpeg('v3v4_samples_microbiota_red_and_total_coverage_purple_stats_10k.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2),mar = c(8,10,7,5), mgp=c(6,2,0))
for (i in locations){
  exploreCoverage(input_type='QIIME', input_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt',
                  locat = i, percent_threshold = 0, remove.outliers = 0,cumulative=FALSE,usredneniye=20,
                  y_axis_length_factor = 1.3,color = 'red')
  exploreCoverage(input_type='total.coverage', input_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/stats_v34.txt',
                  locat = i, percent_threshold = 0, remove.outliers = 0,cumulative=FALSE,usredneniye=20,
                  to.append = T,color = 'purple')
  
}
dev.off()
# V5V6
jpeg('v5v6_samples_microbiota_red_and_total_coverage_purple_stats_10k.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2),mar = c(8,10,7,5), mgp=c(6,2,0))
for (i in locations){
  exploreCoverage(input_type='QIIME', input_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt',
                  locat = i, percent_threshold = 0, remove.outliers = 0,cumulative=FALSE,usredneniye=20,y_axis_length_factor = 1.3,
                  color = 'red')
  exploreCoverage(input_type='total.coverage', input_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/stats_v56.txt',
                  locat = i, percent_threshold = 0, remove.outliers = 0,cumulative=FALSE,usredneniye=20,
                  to.append = T,color = 'purple')
  
}
dev.off()

# Than, for reads assigned to microbiota, coverage window +-1k reads (e.g. 0 to 1999) (usredneniye=2)
# V1V2
jpeg('v1v2_samples_microbiota_coverage_stats_1k.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2),mar = c(8,10,7,5), mgp=c(6,2,0))
for (i in locations){
  exploreCoverage(input_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt',
                        locat = i, percent_threshold = 0, remove.outliers = 0,cumulative=FALSE,usredneniye=2)
}
dev.off()
# V3V4
jpeg('v3v4_samples_microbiota_coverage_stats_1k.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2),mar = c(8,10,7,5), mgp=c(6,2,0))
for (i in locations){
  exploreCoverage(input_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt',
                        locat = i, percent_threshold = 0, remove.outliers = 0,cumulative=FALSE,usredneniye=2)
}
dev.off()
# V5V6
jpeg('v5v6_samples_microbiota_coverage_stats_1k.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2),mar = c(8,10,7,5), mgp=c(6,2,0))
for (i in locations){
  exploreCoverage(input_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt',
                        locat = i, percent_threshold = 0, remove.outliers = 0,cumulative=FALSE,usredneniye=2)
}
dev.off()


#8a. How different are abilities of amplicons to catch taxa
source(paste0(progr.dir,'/Check_data_lit_review_v3.R'))
# here are taxa that have significantly different representation across amplicons
diff12 <- checkTaxaAcrossAmpl(file1 = paste0(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt'), file2=paste0(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt'),sample.cov.filter=c(23000,13000))
diff13 <- checkTaxaAcrossAmpl(file1 = paste0(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt'), file2=paste0(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt'),sample.cov.filter=c(23000,21000))
diff23 <- checkTaxaAcrossAmpl(file1 = paste0(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt'), file2=paste0(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt'),sample.cov.filter=c(13000,21000))
# taxa where we have problems with in least 2 pairs (34 taxa)
union(union(rownames(diff13),rownames(diff12)),rownames(diff23))
# taxa where we have problems in all location pairs (5 taxa) dominated by o__Clostridiales
intersect(intersect(rownames(diff13),rownames(diff12)),rownames(diff23))

# log10 scale
# here are taxa that have significantly different representation across amplicons
diff12.lg <- checkTaxaAcrossAmpl(file1 = paste0(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt'), file2=paste0(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt'),log.scale = T,sample.cov.filter=c(23000,13000))
diff13.lg <- checkTaxaAcrossAmpl(file1 = paste0(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt'), file2=paste0(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt'),log.scale = T,sample.cov.filter=c(23000,21000))
diff23.lg <- checkTaxaAcrossAmpl(file1 = paste0(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt'), file2=paste0(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt'),log.scale = T,sample.cov.filter=c(13000,21000))
# taxa where we have problems with in least 2 pairs 
union(union(rownames(diff13.lg),rownames(diff12.lg)),rownames(diff23.lg))
# taxa where we have problems in all location pairs 
intersect(intersect(rownames(diff13.lg),rownames(diff12.lg)),rownames(diff23.lg))

# 8b plot all samples with signif difference across amplicons of >1 order
setwd(paste0(home.dir,'/Results/Check_liter_review/Amplicon_difference_detect_taxa'))

to.plot.diff12.lg <- diff12.lg[abs(diff12.lg[,3])>=1|diff12.lg[,1]>=-2|diff12.lg[,2]>=-2,1:2]
rownames(to.plot.diff12.lg) <- word(rownames(diff12.lg[abs(diff12.lg[,3])>=1|diff12.lg[,1]>=-2|diff12.lg[,2]>=-2,]),start = -2, end=-1, sep='__')
to.plot.diff13.lg <- diff13.lg[abs(diff13.lg[,3])>=1|diff13.lg[,1]>=-2|diff13.lg[,2]>=-2,1:2]
rownames(to.plot.diff13.lg) <- word(rownames(diff13.lg[abs(diff13.lg[,3])>=1|diff13.lg[,1]>=-2|diff13.lg[,2]>=-2,]),start = -2, end=-1, sep='__')
to.plot.diff23.lg <- diff23.lg[abs(diff23.lg[,3])>=1|diff23.lg[,1]>=-2|diff23.lg[,2]>=-2,1:2]
rownames(to.plot.diff23.lg) <- word(rownames(diff23.lg[abs(diff23.lg[,3])>=1|diff23.lg[,1]>=-2|diff23.lg[,2]>=-2,]),start = -2, end=-1, sep='__')
labels.12 <- seq(1,dim(to.plot.diff12.lg)[1],by=1)
labels.13 <- seq(1,dim(to.plot.diff13.lg)[1],by=1)
labels.23 <- seq(1,dim(to.plot.diff23.lg)[1],by=1)

# V12 vs V34
jpeg('v12-v34_log10_taxa_differ_detect.jpg',width = 1200,height = 1200)
par(mar=c(7,7,7,2))
plot(to.plot.diff12.lg,main='v1v2 vs v3v4 taxa average fraction',xlim = c(-8,0),ylim = c(-8,0),
     xlab = expression('v1v2 taxa average fraction across samples, log'[10]*' scale'), ylab = expression('v3v4 taxa average fraction across samples, log'[10]*' scale'),
     cex.main=5,cex.lab=3,cex.axis=2, col='red',pch=20)
lines(rbind(c(1,1),c(-13,-13)),col = 'green')
lines(rbind(c(0,1),c(-14,-13)),col = 'orange')
lines(rbind(c(1,0),c(-13,-14)),col = 'orange')
grid(lwd=3)
# identify(x=to.plot.diff12.lg[,1],y=to.plot.diff12.lg[,2],labels = rownames(to.plot.diff12.lg))
text(x=to.plot.diff12.lg[,1],y=to.plot.diff12.lg[,2],labels = labels.12,cex=2.5,col='black',pos = 3)
lines(diff12.lg[abs(diff12.lg[,3])<1,1:2],type='p',col='orange')
dev.off()
legend12 = paste0(labels.12,'=',rownames(to.plot.diff12.lg),'\n')
cat(legend12)
# V12 vs v56
jpeg('v12-v56_log10_taxa_differ_detect.jpg',width = 1200,height = 1200)
par(mar=c(7,7,7,2))
plot(to.plot.diff13.lg,main='v1v2 vs v5v6 taxa average fraction',xlim = c(-8,0),ylim = c(-8,0),
     xlab = expression('v1v2 taxa average fraction across samples, log'[10]*' scale'), ylab = expression('v5v6 taxa average fraction across samples, log'[10]*' scale'),
     cex.main=5,cex.lab=3,cex.axis=2, col='blue',pch=19)
lines(rbind(c(1,1),c(-13,-13)),col = 'green')
lines(rbind(c(0,1),c(-14,-13)),col = 'orange')
lines(rbind(c(1,0),c(-13,-14)),col = 'orange')
grid(lwd=3)
text(x=to.plot.diff13.lg[,1],y=to.plot.diff13.lg[,2],labels = labels.13,cex=2.5,col='black',pos = 3)
lines(diff13.lg[abs(diff13.lg[,3])<1,1:2],type='p',col='black',pch=19)
dev.off()
legend13 = paste0(labels.13,'=',rownames(to.plot.diff13.lg),'\n')
cat(legend13)
# v34 vs v56
jpeg('v34-v56_log10_taxa_differ_detect.jpg',width = 1200,height = 1200)
par(mar=c(7,7,7,2))
plot(to.plot.diff23.lg,main='v3v4 vs v5v6 taxa average fraction',xlim = c(-8,0),ylim = c(-8,0),
     xlab = expression('v3v4 taxa average fraction across samples, log'[10]*' scale'), ylab = expression('v5v6 taxa average fraction across samples, log'[10]*' scale'),
     cex.main=5,cex.lab=3,cex.axis=2, col='red')
lines(rbind(c(1,1),c(-13,-13)),col = 'green')
lines(rbind(c(0,1),c(-14,-13)),col = 'orange')
lines(rbind(c(1,0),c(-13,-14)),col = 'orange')
grid(lwd=3)
text(x=to.plot.diff23.lg[,1],y=to.plot.diff23.lg[,2],labels = labels.23,cex=2.5,col='black',pos = 3)
lines(diff23.lg[abs(diff23.lg[,3])<1,1:2],type='p',col='orange')
dev.off()
legend23 = paste0(labels.23,'=',rownames(to.plot.diff23.lg),'\n')
cat(legend23)

# Now see manually which bacteria are 'mis-estimated' by a particular amplicon
# V1V2 - and choose manually which have same-side difference for both pairs
setdiff(intersect(rownames(to.plot.diff12.lg),rownames(to.plot.diff13.lg)),rownames(to.plot.diff23.lg))
# V3V4 - and choose manually which have same-side difference for both pairs
setdiff(intersect(rownames(to.plot.diff12.lg),rownames(to.plot.diff23.lg)),rownames(to.plot.diff13.lg))
# V5V6 - and choose manually which have same-side difference for both pairs
setdiff(intersect(rownames(to.plot.diff13.lg),rownames(to.plot.diff23.lg)),rownames(to.plot.diff12.lg))


# 9. Descriptive statistics about the initial sample
# total number of samples per location
# number of replicates
input_file1 = file.path(home.dir,'Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6_g_Bact_Prev_removed_v2.txt')
input_file2 = file.path(home.dir,'Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6_g_Bact_Prev_removed_v2.txt')
input_file3 = file.path(home.dir,'Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6_g_Bact_Prev_removed_v2.txt')

data1=read.table(input_file1, header=T, row.names=1, dec=".", sep="\t")
data2=read.table(input_file2, header=T, row.names=1, dec=".", sep="\t")
data3=read.table(input_file3, header=T, row.names=1, dec=".", sep="\t")

# get the number of unique samples 
samples <- gsub('\\.IL|\\.TR|\\.RE|\\.V[0-9][0-9]|\\.1','',union(union(colnames(data1),colnames(data2)),colnames(data3)))
# samples <- samples[!samples%like%'\\.[0-9]']
samples <- unique(samples)
print(paste0('The total number of unique samples in V12, V34 and V56: ',length(samples)))
samples.IL <- gsub('\\.IL|\\.TR|\\.RE|\\.V[0-9][0-9]|\\.1','',union(union(colnames(data1)[colnames(data1)%like%'\\.IL'], colnames(data2)[colnames(data2)%like%'\\.IL']), colnames(data3)[colnames(data3)%like%'\\.IL']))
# samples.IL <- samples.IL[!samples.IL%like%'\\.[0-9]']
samples.IL <- unique(samples.IL)
print(paste0('The total number of unique samples in V12, V34 and V56 in IL: ',length(samples.IL)))
samples.TR <- gsub('\\.IL|\\.TR|\\.RE|\\.V[0-9][0-9]|\\.1','',union(union(colnames(data1)[colnames(data1)%like%'\\.TR'], colnames(data2)[colnames(data2)%like%'\\.TR']), colnames(data3)[colnames(data3)%like%'\\.TR']))
# samples.TR <- samples.TR[!samples.TR%like%'\\.[0-9]']
samples.TR <- unique(samples.TR)
print(paste0('The total number of unique samples in V12, V34 and V56 in TR: ',length(samples.TR)))
samples.RE <- gsub('\\.IL|\\.TR|\\.RE|\\.V[0-9][0-9]|\\.1','',union(union(colnames(data1)[colnames(data1)%like%'\\.RE'], colnames(data2)[colnames(data2)%like%'\\.RE']), colnames(data3)[colnames(data3)%like%'\\.RE']))
# samples.RE <- samples.RE[!samples.RE%like%'\\.[0-9]']
samples.RE <- unique(samples.RE)
print(paste0('The total number of unique samples in V12, V34 and V56: ',length(samples.RE)))


# get number of replicates as well as total number of samples per location
# for IL, 20 replicates
IL <- colnames(data1)[colnames(data1)%like%'\\.IL']
print(paste0('The total number of replicate samples in IL: ', length(IL[IL%like%'\\.[0-9]'])))
print(paste0('The total number of samples in IL including replicates: ', length(IL)))
print(paste0('The total number of samples in IL excluding replicates: ', length(IL[!IL%like%'\\.[0-9]'])))
# for TR, 0 replicates
TR <- colnames(data1)[colnames(data1)%like%'\\.TR']
print(paste0('The total number of replicate samples in TR: ', length(TR[TR%like%'\\.[0-9]'])))
print(paste0('The total number of samples in TR including replicates: ', length(TR)))
print(paste0('The total number of samples in TR excluding replicates: ', length(TR[!TR%like%'\\.[0-9]'])))
# for RE, 2 replicates
RE <- colnames(data1)[colnames(data1)%like%'\\.RE']
print(paste0('The total number of replicate samples in RE: ', length(RE[RE%like%'\\.[0-9]'])))
print(paste0('The total number of samples in RE including replicates: ', length(RE)))
print(paste0('The total number of samples in RE excluding replicates: ', length(RE[!RE%like%'\\.[0-9]'])))

