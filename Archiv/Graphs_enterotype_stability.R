# This program draws a scatter plot with points representing enterotypes. Set1 - stable across amplicons
# set 2 - stable across couple of amplicons: 12-34, 34-56, 12-56
# set 3 - not stable across amplicons at all
graphEntStabil <- function(file1 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/V1V2_3e_IL_sample-enterotype.csv', file2 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/V3V4_3e_IL_sample-enterotype.csv', file3 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/V5V6_3e_IL_sample-enterotype.csv'){
library(ade4)
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/GetStablyEnterotypedSamples_across_amplicons_smart.R')

# do the main job of enterotype classification. We get 5 data frames
classify.IL <- getStabAcrAmpl(file1,file2,file3)
  
# Get samples stable across all amplicons for IL

# all th

# Get samples stable across 2 of 3 locations for IL

# Get samples not


ent.Raes12.3e.sep.by.loc.IL <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt', '', clust.number=3, noiserem=T, locat='.IL', plots=F)
ent.Raes12.3e.sep.by.loc.IL$sample.enterotype[,1] <- gsub('\\.V[0-9][0-9]','',ent.Raes12.3e.sep.by.loc.IL$sample.enterotype[,1])

# replace enterotype numbers by B, P and 3
levels(ent.Raes12.3e.sep.by.loc.IL$sample.enterotype[,2])[levels(ent.Raes12.3e.sep.by.loc.IL$sample.enterotype[,2])%like%'g__Bacteroides'] <- 'B'
levels(ent.Raes12.3e.sep.by.loc.IL$sample.enterotype[,2])[levels(ent.Raes12.3e.sep.by.loc.IL$sample.enterotype[,2])%like%'g__Prevotella'] <- 'P'
levels(ent.Raes12.3e.sep.by.loc.IL$sample.enterotype[,2])[!levels(ent.Raes12.3e.sep.by.loc.IL$sample.enterotype[,2])%in%c('P','B')] <- 'e3'

# rows of samples that are stable for all amplicons in IL
stable.IL.12_34_56_numbers <- ent.Raes12.3e.sep.by.loc.IL$sample.enterotype[,1] %in% stable.ent.IL[,1]

# drawing scatter plot
# samples stable across 3 amplicons
s.class(ent.Raes12.3e.sep.by.loc.IL$obs_bet_result$ls[stable.IL.12_34_56_numbers,], fac=as.factor(ent.Raes12.3e.sep.by.loc.IL$sample.enterotype$enterotype[stable.IL.12_34_56_numbers]), col=c('red','red','red','blue'), cpoint = 1, grid=F, sub="Between-class analysis for the enterotypes",cellipse = 1)

# s.class(ent.Raes12.3e.sep.by.loc.IL$obs_bet_result$ls[51:150,], fac=as.factor(ent.Raes12.3e.sep.by.loc.IL$sample.enterotype$enterotype[51:150]), col=c('blue','blue','blue','blue'), cpoint = 1, grid=F, sub="Between-class analysis for the enterotypes",cellipse = 0, add.plot = T)

}