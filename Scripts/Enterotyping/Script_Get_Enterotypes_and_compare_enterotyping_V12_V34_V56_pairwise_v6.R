# Make enterotyping for 3 locations for 2-5 enterotypes 
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Enterotyping_Raes_noise_removed_v8.R')
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Enterotyping_compar_functions_v5.R')
setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype')
# first, make enterotyping
# 2 enterotypes
ent.Raes12.2e <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=2, noiserem=T, plots=F)
ent.Raes34.2e <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V3V4', clust.number=2, noiserem=T, plots=F)
ent.Raes56.2e <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V5V6', clust.number=2, noiserem=T, plots=F)
write.table(ent.Raes12.2e$sample.enterotype, 'V1V2_2e_sample-enterotype.csv')
write.table(ent.Raes34.2e$sample.enterotype, 'V3V4_2e_sample-enterotype.csv')
write.table(ent.Raes56.2e$sample.enterotype, 'V5V6_2e_sample-enterotype.csv')
# 3 enterotypes
ent.Raes12.3e <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=3, noiserem=T, plots=F)
ent.Raes34.3e <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V3V4', clust.number=3, noiserem=T, plots=F)
ent.Raes56.3e <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V5V6', clust.number=3, noiserem=T, plots=F)
write.table(ent.Raes12.3e$sample.enterotype, 'V1V2_3e_sample-enterotype.csv')
write.table(ent.Raes34.3e$sample.enterotype, 'V3V4_3e_sample-enterotype.csv')
write.table(ent.Raes56.3e$sample.enterotype, 'V5V6_3e_sample-enterotype.csv')
# 4 enterotypes
ent.Raes12.4e <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=4, noiserem=T, plots=F)
ent.Raes34.4e <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V3V4', clust.number=4, noiserem=T, plots=F)
ent.Raes56.4e <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V5V6', clust.number=4, noiserem=T, plots=F)
write.table(ent.Raes12.4e$sample.enterotype, 'V1V2_4e_sample-enterotype.csv')
write.table(ent.Raes34.4e$sample.enterotype, 'V3V4_4e_sample-enterotype.csv')
write.table(ent.Raes56.4e$sample.enterotype, 'V5V6_4e_sample-enterotype.csv')
# 5 enterotypes
ent.Raes12.5e <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=5, noiserem=T, plots=F)
ent.Raes34.5e <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V3V4', clust.number=5, noiserem=T, plots=F)
ent.Raes56.5e <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V5V6', clust.number=5, noiserem=T, plots=F)
write.table(ent.Raes12.5e$sample.enterotype, 'V1V2_5e_sample-enterotype.csv')
write.table(ent.Raes34.5e$sample.enterotype, 'V3V4_5e_sample-enterotype.csv')
write.table(ent.Raes56.5e$sample.enterotype, 'V5V6_5e_sample-enterotype.csv')

#now get number of common elements between v12-v34-v56 clustering strategies and enterotype matches
# 2 enterotypes
compEnt_v12_v34.2e <- CompEntPred.pair(ent.Raes12.2e[[2]],ent.Raes34.2e[[2]],'CompEntPred_2e_V12_V34.jpg')
compEnt_v34_v56.2e <- CompEntPred.pair(ent.Raes34.2e[[2]],ent.Raes56.2e[[2]],'CompEntPred_2e_V34_V56.jpg')
compEnt_v12_v56.2e <- CompEntPred.pair(ent.Raes12.2e[[2]],ent.Raes56.2e[[2]],'CompEntPred_2e_V12_V56.jpg')

# 3 enterotypes
compEnt_v12_v34.3e <- CompEntPred.pair(ent.Raes12.3e[[2]],ent.Raes34.3e[[2]],'CompEntPred_3e_V12_V34.jpg')
compEnt_v34_v56.3e <- CompEntPred.pair(ent.Raes34.3e[[2]],ent.Raes56.3e[[2]],'CompEntPred_3e_V34_V56.jpg')
compEnt_v12_v56.3e <- CompEntPred.pair(ent.Raes12.3e[[2]],ent.Raes56.3e[[2]],'CompEntPred_3e_V12_V56.jpg')

# 4 enterotypes
compEnt_v12_v34.4e <- CompEntPred.pair(ent.Raes12.4e[[2]],ent.Raes34.4e[[2]],'CompEntPred_4e_V12_V34.jpg')
compEnt_v34_v56.4e <- CompEntPred.pair(ent.Raes34.4e[[2]],ent.Raes56.4e[[2]],'CompEntPred_4e_V34_V56.jpg')
compEnt_v12_v56.4e <- CompEntPred.pair(ent.Raes12.4e[[2]],ent.Raes56.4e[[2]],'CompEntPred_4e_V12_V56.jpg')

# 5 enterotypes
compEnt_v12_v34.5e <- CompEntPred.pair(ent.Raes12.5e$sample.enterotype,ent.Raes34.5e$sample.enterotype,'CompEntPred_5e_V12_V34.jpg')
compEnt_v34_v56.5e <- CompEntPred.pair(ent.Raes34.5e$sample.enterotype,ent.Raes56.5e$sample.enterotype,'CompEntPred_5e_V34_V56.jpg')
compEnt_v12_v56.5e <- CompEntPred.pair(ent.Raes12.5e$sample.enterotype,ent.Raes56.5e$sample.enterotype,'CompEntPred_5e_V12_V56.jpg')

# Part 2 (written 2018-10-09)
# AIM: Split enterotypes by locations
setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype')
# read the tables with enterotypes (obtained above)
v12_2e <- read.table('V1V2_2e_sample-enterotype.csv',header = T)
v12_3e <- read.table('V1V2_3e_sample-enterotype.csv',header = T)
v12_4e <- read.table('V1V2_4e_sample-enterotype.csv',header = T)
v12_5e <- read.table('V1V2_5e_sample-enterotype.csv',header = T)

v34_2e <- read.table('V3V4_2e_sample-enterotype.csv',header = T)
v34_3e <- read.table('V3V4_3e_sample-enterotype.csv',header = T)
v34_4e <- read.table('V3V4_4e_sample-enterotype.csv',header = T)
v34_5e <- read.table('V3V4_5e_sample-enterotype.csv',header = T)

v56_2e <- read.table('V5V6_2e_sample-enterotype.csv',header = T)
v56_3e <- read.table('V5V6_3e_sample-enterotype.csv',header = T)
v56_4e <- read.table('V5V6_4e_sample-enterotype.csv',header = T)
v56_5e <- read.table('V5V6_5e_sample-enterotype.csv',header = T)

library(data.table)
# 2 enterotypes
v12_2e.IL <- v12_2e[v12_2e[,1] %like% '\\.IL',]
v12_2e.TR <- v12_2e[v12_2e[,1] %like% '\\.TR',]
v12_2e.RE <- v12_2e[v12_2e[,1] %like% '\\.RE',]
v34_2e.IL <- v34_2e[v34_2e[,1] %like% '\\.IL',]
v34_2e.TR <- v34_2e[v34_2e[,1] %like% '\\.TR',]
v34_2e.RE <- v34_2e[v34_2e[,1] %like% '\\.RE',]
v56_2e.IL <- v56_2e[v56_2e[,1] %like% '\\.IL',]
v56_2e.TR <- v56_2e[v56_2e[,1] %like% '\\.TR',]
v56_2e.RE <- v56_2e[v56_2e[,1] %like% '\\.RE',]
# 3 enterotypes
v12_3e.IL <- v12_3e[v12_3e[,1] %like% '\\.IL',]
v12_3e.TR <- v12_3e[v12_3e[,1] %like% '\\.TR',]
v12_3e.RE <- v12_3e[v12_3e[,1] %like% '\\.RE',]
v34_3e.IL <- v34_3e[v34_3e[,1] %like% '\\.IL',]
v34_3e.TR <- v34_3e[v34_3e[,1] %like% '\\.TR',]
v34_3e.RE <- v34_3e[v34_3e[,1] %like% '\\.RE',]
v56_3e.IL <- v56_3e[v56_3e[,1] %like% '\\.IL',]
v56_3e.TR <- v56_3e[v56_3e[,1] %like% '\\.TR',]
v56_3e.RE <- v56_3e[v56_3e[,1] %like% '\\.RE',]
# 4 enterotypes
v12_4e.IL <- v12_4e[v12_4e[,1] %like% '\\.IL',]
v12_4e.TR <- v12_4e[v12_4e[,1] %like% '\\.TR',]
v12_4e.RE <- v12_4e[v12_4e[,1] %like% '\\.RE',]
v34_4e.IL <- v34_4e[v34_4e[,1] %like% '\\.IL',]
v34_4e.TR <- v34_4e[v34_4e[,1] %like% '\\.TR',]
v34_4e.RE <- v34_4e[v34_4e[,1] %like% '\\.RE',]
v56_4e.IL <- v56_4e[v56_4e[,1] %like% '\\.IL',]
v56_4e.TR <- v56_4e[v56_4e[,1] %like% '\\.TR',]
v56_4e.RE <- v56_4e[v56_4e[,1] %like% '\\.RE',]
# 5 enterotypes
v12_5e.IL <- v12_5e[v12_5e[,1] %like% '\\.IL',]
v12_5e.TR <- v12_5e[v12_5e[,1] %like% '\\.TR',]
v12_5e.RE <- v12_5e[v12_5e[,1] %like% '\\.RE',]
v34_5e.IL <- v34_5e[v34_5e[,1] %like% '\\.IL',]
v34_5e.TR <- v34_5e[v34_5e[,1] %like% '\\.TR',]
v34_5e.RE <- v34_5e[v34_5e[,1] %like% '\\.RE',]
v56_5e.IL <- v56_5e[v56_5e[,1] %like% '\\.IL',]
v56_5e.TR <- v56_5e[v56_5e[,1] %like% '\\.TR',]
v56_5e.RE <- v56_5e[v56_5e[,1] %like% '\\.RE',]

# now, write to files
setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype/Split_by_locat')
write.table(v12_2e.IL,'V1V2_2e_IL_sample-enterotype.csv')
write.table(v12_2e.TR,'V1V2_2e_TR_sample-enterotype.csv')
write.table(v12_2e.RE,'V1V2_2e_RE_sample-enterotype.csv')
write.table(v34_2e.IL,'V3V4_2e_IL_sample-enterotype.csv')
write.table(v34_2e.TR,'V3V4_2e_TR_sample-enterotype.csv')
write.table(v34_2e.RE,'V3V4_2e_RE_sample-enterotype.csv')
write.table(v56_2e.IL,'V5V6_2e_IL_sample-enterotype.csv')
write.table(v56_2e.TR,'V5V6_2e_TR_sample-enterotype.csv')
write.table(v56_2e.RE,'V5V6_2e_RE_sample-enterotype.csv')

write.table(v12_3e.IL,'V1V2_3e_IL_sample-enterotype.csv')
write.table(v12_3e.TR,'V1V2_3e_TR_sample-enterotype.csv')
write.table(v12_3e.RE,'V1V2_3e_RE_sample-enterotype.csv')
write.table(v34_3e.IL,'V3V4_3e_IL_sample-enterotype.csv')
write.table(v34_3e.TR,'V3V4_3e_TR_sample-enterotype.csv')
write.table(v34_3e.RE,'V3V4_3e_RE_sample-enterotype.csv')
write.table(v56_3e.IL,'V5V6_3e_IL_sample-enterotype.csv')
write.table(v56_3e.TR,'V5V6_3e_TR_sample-enterotype.csv')
write.table(v56_3e.RE,'V5V6_3e_RE_sample-enterotype.csv')

write.table(v12_4e.IL,'V1V2_4e_IL_sample-enterotype.csv')
write.table(v12_4e.TR,'V1V2_4e_TR_sample-enterotype.csv')
write.table(v12_4e.RE,'V1V2_4e_RE_sample-enterotype.csv')
write.table(v34_4e.IL,'V3V4_4e_IL_sample-enterotype.csv')
write.table(v34_4e.TR,'V3V4_4e_TR_sample-enterotype.csv')
write.table(v34_4e.RE,'V3V4_4e_RE_sample-enterotype.csv')
write.table(v56_4e.IL,'V5V6_4e_IL_sample-enterotype.csv')
write.table(v56_4e.TR,'V5V6_4e_TR_sample-enterotype.csv')
write.table(v56_4e.RE,'V5V6_4e_RE_sample-enterotype.csv')

write.table(v12_5e.IL,'V1V2_5e_IL_sample-enterotype.csv')
write.table(v12_5e.TR,'V1V2_5e_TR_sample-enterotype.csv')
write.table(v12_5e.RE,'V1V2_5e_RE_sample-enterotype.csv')
write.table(v34_5e.IL,'V3V4_5e_IL_sample-enterotype.csv')
write.table(v34_5e.TR,'V3V4_5e_TR_sample-enterotype.csv')
write.table(v34_5e.RE,'V3V4_5e_RE_sample-enterotype.csv')
write.table(v56_5e.IL,'V5V6_5e_IL_sample-enterotype.csv')
write.table(v56_5e.TR,'V5V6_5e_TR_sample-enterotype.csv')
write.table(v56_5e.RE,'V5V6_5e_RE_sample-enterotype.csv')
