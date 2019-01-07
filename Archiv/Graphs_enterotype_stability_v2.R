# This program draws a scatter plot with points representing enterotypes. Set1 - stable across amplicons
# set 2 - stable across couple of amplicons: 12-34, 34-56, 12-56
# set 3 - not stable across amplicons at all
plotEntStabil <- function(file1 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/V1V2_3e_IL_sample-enterotype.csv', file2 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/V3V4_3e_IL_sample-enterotype.csv', file3 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/V5V6_3e_IL_sample-enterotype.csv',qiime_out = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt',outdir='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs/BCA',location='.IL'){
library(ade4)
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/GetStablyEnterotypedSamples_across_amplicons_smart_v2.R')

# do the main job of enterotype classification. We get 5 data frames
classify.IL <- getStabAcrAmpl(file1,file2,file3)
  
# Get samples stable across all amplicons for IL

# Get samples stable across 2 of 3 locations for IL

ent.Raes12.3e.sep.by.loc.IL <- Enterotyping_Raes(qiime_out, '', clust.number=3, noiserem=T, locat=location, plots=F)
sample.enterotype <- ent.Raes12.3e.sep.by.loc.IL$sample.enterotype
sample.enterotype[,1] <- gsub('\\.V[0-9][0-9]','',sample.enterotype[,1])

# replace enterotype numbers by B, P and 3
levels(sample.enterotype[,2])[levels(sample.enterotype[,2])%like%'g__Bacteroides'] <- 'B'
levels(sample.enterotype[,2])[levels(sample.enterotype[,2])%like%'g__Prevotella'] <- 'P'
levels(sample.enterotype[,2])[!levels(sample.enterotype[,2])%in%c('P','B')] <- 'e3'

# rows of samples that are stable for all amplicons in IL
stable.IL.12_34_56_numbers <- sample.enterotype[,1] %in% classify.IL$result.v12v34v56[,1]
stable.IL.12_34_numbers <- sample.enterotype[,1] %in% classify.IL$result.v12v34[,1]
stable.IL.12_56_numbers <- sample.enterotype[,1] %in% classify.IL$result.v12v56[,1]
stable.IL.34_56_numbers <- sample.enterotype[,1] %in% classify.IL$result.v34v56[,1]
unstable.IL_numbers <- sample.enterotype[,1] %in% classify.IL$result.nonconserved[,1]

# drawing scatter plot
# samples stable across 3 amplicons

vector_tmp <- word(word(file1,start = -1,sep = '/'),start=c(1,2,3),sep = '_')
vector_tmp1 <- paste(vector_tmp[1], vector_tmp[2], vector_tmp[3], sep ='_')
file_out <- paste0(vector_tmp1,'_bca_first_2_eignenvectors.jpg')

setwd(outdir)
jpeg(filename = file_out,width = 2000,height = 1600)
legend <- 'Red - stable all ampl., yellow - stable v12-v34, orange v12-v56, green v34-56, blue - unstable'
xlim = c(min(ent.Raes12.3e.sep.by.loc.IL$obs_bet_result$ls[,1])-0.5,max(ent.Raes12.3e.sep.by.loc.IL$obs_bet_result$ls[,1])+0.5)
ylim = c(min(ent.Raes12.3e.sep.by.loc.IL$obs_bet_result$ls[,2])-1,max(ent.Raes12.3e.sep.by.loc.IL$obs_bet_result$ls[,2])+1)
                          
s.class(ent.Raes12.3e.sep.by.loc.IL$obs_bet_result$ls[stable.IL.12_34_56_numbers,], fac=as.factor(sample.enterotype$enterotype[stable.IL.12_34_56_numbers]), col=c('red','red','red','red'), xlim=xlim, ylim=ylim, cpoint = 2, grid=F, sub=legend, cstar = 1, csub = 4, cellipse = 1)
s.class(ent.Raes12.3e.sep.by.loc.IL$obs_bet_result$ls[stable.IL.12_34_numbers,], fac=as.factor(sample.enterotype$enterotype[stable.IL.12_34_numbers]), col=c('orange','orange','orange','orange'), cpoint = 2, grid=F, cellipse = 0, cstar = 0, label = NULL, add.plot = T)
s.class(ent.Raes12.3e.sep.by.loc.IL$obs_bet_result$ls[stable.IL.12_56_numbers,], fac=as.factor(sample.enterotype$enterotype[stable.IL.12_56_numbers]), col=c('yellow','yellow','yellow','yellow'), cpoint = 2, grid=F, cellipse = 0, cstar = 0, label = NULL, add.plot = T)
s.class(ent.Raes12.3e.sep.by.loc.IL$obs_bet_result$ls[stable.IL.34_56_numbers,], fac=as.factor(sample.enterotype$enterotype[stable.IL.34_56_numbers]), col=c('green','green','green','green'), cpoint = 2, grid=F, cellipse = 0, cstar = 0, label = NULL, add.plot = T)
s.class(ent.Raes12.3e.sep.by.loc.IL$obs_bet_result$ls[unstable.IL_numbers,], fac=as.factor(sample.enterotype$enterotype[unstable.IL_numbers]), col=c('blue','blue','blue','blue'), cpoint = 2, grid=F, cellipse = 0, cstar = 0, label = NULL, add.plot = T)
dev.off()

# s.label(ent.Raes12.3e.sep.by.loc.IL$obs_bet_result$ls[stable.IL.12_34_56_numbers,],sub=legend, csub = 1.5, possub = "bottmoright")

}