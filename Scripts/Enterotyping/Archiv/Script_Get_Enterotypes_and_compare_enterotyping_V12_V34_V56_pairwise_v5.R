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


