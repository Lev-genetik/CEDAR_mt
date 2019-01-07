# script to compare enterotyping for IL|RE|TR separately and create 3 matrices for pairwise consistency between 16S amplicons (as in 
# report for 06.09.18 to Michel)
# rows = locations: All_loc, IL, RE, TR
# colons = cluster number: 2, 3 or 4

source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Enterotyping_Raes_noise_removed_v9.R')
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Enterotyping_compar_functions_v5.R')

# All locations
ent.Raes12.2e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=2, locat='all', noiserem=T, plots=F)
ent.Raes34.2e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=2, locat='all', noiserem=T, plots=F)
ent.Raes56.2e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=2, locat='all', noiserem=T, plots=F)

ent.Raes12.3e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=3, locat='all', noiserem=T, plots=F)
ent.Raes34.3e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=3, locat='all', noiserem=T, plots=F)
ent.Raes56.3e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=3, locat='all', noiserem=T, plots=F)

ent.Raes12.4e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=4, locat='all', noiserem=T, plots=F)
ent.Raes34.4e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=4, locat='all', noiserem=T, plots=F)
ent.Raes56.4e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=4, locat='all', noiserem=T, plots=F)

# Ileum
ent.Raes12.IL.2e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=2, locat='IL', noiserem=T, plots=F)
ent.Raes34.IL.2e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=2, locat='IL', noiserem=T, plots=F)
ent.Raes56.IL.2e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=2, locat='IL', noiserem=T, plots=F)

ent.Raes12.IL.3e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=3, locat='IL', noiserem=T, plots=F)
ent.Raes34.IL.3e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=3, locat='IL', noiserem=T, plots=F)
ent.Raes56.IL.3e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=3, locat='IL', noiserem=T, plots=F)

ent.Raes12.IL.4e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=4, locat='IL', noiserem=T, plots=F)
ent.Raes34.IL.4e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=4, locat='IL', noiserem=T, plots=F)
ent.Raes56.IL.4e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=4, locat='IL', noiserem=T, plots=F)

# Transverse colon
ent.Raes12.TR.2e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=2, locat='TR', noiserem=T, plots=F)
ent.Raes34.TR.2e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=2, locat='TR', noiserem=T, plots=F)
ent.Raes56.TR.2e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=2, locat='TR', noiserem=T, plots=F)

ent.Raes12.TR.3e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=3, locat='TR', noiserem=T, plots=F)
ent.Raes34.TR.3e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=3, locat='TR', noiserem=T, plots=F)
ent.Raes56.TR.3e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=3, locat='TR', noiserem=T, plots=F)

ent.Raes12.TR.4e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=4, locat='TR', noiserem=T, plots=F)
ent.Raes34.TR.4e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=4, locat='TR', noiserem=T, plots=F)
ent.Raes56.TR.4e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=4, locat='TR', noiserem=T, plots=F)

# Rectum
ent.Raes12.RE.2e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=2, locat='RE', noiserem=T, plots=F)
ent.Raes34.RE.2e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=2, locat='RE', noiserem=T, plots=F)
ent.Raes56.RE.2e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=2, locat='RE', noiserem=T, plots=F)

ent.Raes12.RE.3e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=3, locat='RE', noiserem=T, plots=F)
ent.Raes34.RE.3e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=3, locat='RE', noiserem=T, plots=F)
ent.Raes56.RE.3e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=3, locat='RE', noiserem=T, plots=F)

ent.Raes12.RE.4e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=4, locat='RE', noiserem=T, plots=F)
ent.Raes34.RE.4e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=4, locat='RE', noiserem=T, plots=F)
ent.Raes56.RE.4e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Graphs/V1V2', clust.number=4, locat='RE', noiserem=T, plots=F)

# Now compare clustering - IL
compEnt_v12_v34.IL.4e <- CompEntPred.pair(ent.Raes12.IL.4e$sample.enterotype,ent.Raes34.IL.4e$sample.enterotype,'CompEntPred_4e_V12_V34_IL.jpg')
compEnt_v34_v56.IL.4e <- CompEntPred.pair(ent.Raes34.IL.4e$sample.enterotype,ent.Raes56.IL.4e$sample.enterotype,'CompEntPred_4e_V34_V56_IL.jpg')
compEnt_v12_v56.IL.4e <- CompEntPred.pair(ent.Raes12.IL.4e$sample.enterotype,ent.Raes56.IL.4e$sample.enterotype,'CompEntPred_4e_V12_V56_IL.jpg')

compEnt_v12_v34.IL.3e <- CompEntPred.pair(ent.Raes12.IL.3e$sample.enterotype,ent.Raes34.IL.3e$sample.enterotype,'CompEntPred_3e_V12_V34_IL.jpg')
compEnt_v34_v56.IL.3e <- CompEntPred.pair(ent.Raes34.IL.3e$sample.enterotype,ent.Raes56.IL.3e$sample.enterotype,'CompEntPred_3e_V34_V56_IL.jpg')
compEnt_v12_v56.IL.3e <- CompEntPred.pair(ent.Raes12.IL.3e$sample.enterotype,ent.Raes56.IL.3e$sample.enterotype,'CompEntPred_3e_V12_V56_IL.jpg')

compEnt_v12_v34.IL.2e <- CompEntPred.pair(ent.Raes12.IL.2e$sample.enterotype,ent.Raes34.IL.2e$sample.enterotype,'CompEntPred_2e_V12_V34_IL.jpg')
compEnt_v34_v56.IL.2e <- CompEntPred.pair(ent.Raes34.IL.2e$sample.enterotype,ent.Raes56.IL.2e$sample.enterotype,'CompEntPred_2e_V34_V56_IL.jpg')
compEnt_v12_v56.IL.2e <- CompEntPred.pair(ent.Raes12.IL.2e$sample.enterotype,ent.Raes56.IL.2e$sample.enterotype,'CompEntPred_2e_V12_V56_IL.jpg')

# Now compare clustering - TR
compEnt_v12_v34.TR.4e <- CompEntPred.pair(ent.Raes12.TR.4e$sample.enterotype,ent.Raes34.TR.4e$sample.enterotype,'CompEntPred_4e_V12_V34_TR.jpg')
compEnt_v34_v56.TR.4e <- CompEntPred.pair(ent.Raes34.TR.4e$sample.enterotype,ent.Raes56.TR.4e$sample.enterotype,'CompEntPred_4e_V34_V56_TR.jpg')
compEnt_v12_v56.TR.4e <- CompEntPred.pair(ent.Raes12.TR.4e$sample.enterotype,ent.Raes56.TR.4e$sample.enterotype,'CompEntPred_4e_V12_V56_TR.jpg')

compEnt_v12_v34.TR.3e <- CompEntPred.pair(ent.Raes12.TR.3e$sample.enterotype,ent.Raes34.TR.3e$sample.enterotype,'CompEntPred_3e_V12_V34_TR.jpg')
compEnt_v34_v56.TR.3e <- CompEntPred.pair(ent.Raes34.TR.3e$sample.enterotype,ent.Raes56.TR.3e$sample.enterotype,'CompEntPred_3e_V34_V56_TR.jpg')
compEnt_v12_v56.TR.3e <- CompEntPred.pair(ent.Raes12.TR.3e$sample.enterotype,ent.Raes56.TR.3e$sample.enterotype,'CompEntPred_3e_V12_V56_TR.jpg')

compEnt_v12_v34.TR.2e <- CompEntPred.pair(ent.Raes12.TR.2e$sample.enterotype,ent.Raes34.TR.2e$sample.enterotype,'CompEntPred_2e_V12_V34_TR.jpg')
compEnt_v34_v56.TR.2e <- CompEntPred.pair(ent.Raes34.TR.2e$sample.enterotype,ent.Raes56.TR.2e$sample.enterotype,'CompEntPred_2e_V34_V56_TR.jpg')
compEnt_v12_v56.TR.2e <- CompEntPred.pair(ent.Raes12.TR.2e$sample.enterotype,ent.Raes56.TR.2e$sample.enterotype,'CompEntPred_2e_V12_V56_TR.jpg')

# Now compare clustering - RE
compEnt_v12_v34.RE.4e <- CompEntPred.pair(ent.Raes12.RE.4e$sample.enterotype,ent.Raes34.RE.4e$sample.enterotype,'CompEntPred_4e_V12_V34_RE.jpg')
compEnt_v34_v56.RE.4e <- CompEntPred.pair(ent.Raes34.RE.4e$sample.enterotype,ent.Raes56.RE.4e$sample.enterotype,'CompEntPred_4e_V34_V56_RE.jpg')
compEnt_v12_v56.RE.4e <- CompEntPred.pair(ent.Raes12.RE.4e$sample.enterotype,ent.Raes56.RE.4e$sample.enterotype,'CompEntPred_4e_V12_V56_RE.jpg')

compEnt_v12_v34.RE.3e <- CompEntPred.pair(ent.Raes12.RE.3e$sample.enterotype,ent.Raes34.RE.3e$sample.enterotype,'CompEntPred_3e_V12_V34_RE.jpg')
compEnt_v34_v56.RE.3e <- CompEntPred.pair(ent.Raes34.RE.3e$sample.enterotype,ent.Raes56.RE.3e$sample.enterotype,'CompEntPred_3e_V34_V56_RE.jpg')
compEnt_v12_v56.RE.3e <- CompEntPred.pair(ent.Raes12.RE.3e$sample.enterotype,ent.Raes56.RE.3e$sample.enterotype,'CompEntPred_3e_V12_V56_RE.jpg')

compEnt_v12_v34.RE.2e <- CompEntPred.pair(ent.Raes12.RE.2e$sample.enterotype,ent.Raes34.RE.2e$sample.enterotype,'CompEntPred_2e_V12_V34_RE.jpg')
compEnt_v34_v56.RE.2e <- CompEntPred.pair(ent.Raes34.RE.2e$sample.enterotype,ent.Raes56.RE.2e$sample.enterotype,'CompEntPred_2e_V34_V56_RE.jpg')
compEnt_v12_v56.RE.2e <- CompEntPred.pair(ent.Raes12.RE.2e$sample.enterotype,ent.Raes56.RE.2e$sample.enterotype,'CompEntPred_2e_V12_V56_RE.jpg')

#now get number of common elements between v12-v34-v56 clustering strategies and enterotype matches
compEnt_v12_v34.2e <- CompEntPred.pair(ent.Raes12.2e[[2]],ent.Raes34.2e[[2]],'CompEntPred_2e_V12_V34.jpg')
compEnt_v34_v56.2e <- CompEntPred.pair(ent.Raes34.2e[[2]],ent.Raes56.2e[[2]],'CompEntPred_2e_V34_V56.jpg')
compEnt_v12_v56.2e <- CompEntPred.pair(ent.Raes12.2e[[2]],ent.Raes56.2e[[2]],'CompEntPred_2e_V12_V56.jpg')

compEnt_v12_v34.3e <- CompEntPred.pair(ent.Raes12.3e[[2]],ent.Raes34.3e[[2]],'CompEntPred_3e_V12_V34.jpg')
compEnt_v34_v56.3e <- CompEntPred.pair(ent.Raes34.3e[[2]],ent.Raes56.3e[[2]],'CompEntPred_3e_V34_V56.jpg')
compEnt_v12_v56.3e <- CompEntPred.pair(ent.Raes12.3e[[2]],ent.Raes56.3e[[2]],'CompEntPred_3e_V12_V56.jpg')

compEnt_v12_v34.4e <- CompEntPred.pair(ent.Raes12.4e[[2]],ent.Raes34.4e[[2]],'CompEntPred_4e_V12_V34.jpg')
compEnt_v34_v56.4e <- CompEntPred.pair(ent.Raes34.4e[[2]],ent.Raes56.4e[[2]],'CompEntPred_4e_V34_V56.jpg')
compEnt_v12_v56.4e <- CompEntPred.pair(ent.Raes12.4e[[2]],ent.Raes56.4e[[2]],'CompEntPred_4e_V12_V56.jpg')

# Now compare clustering - all locations
compEnt_v12_v34.2e <- CompEntPred.pair(ent.Raes12.2e[[2]],ent.Raes34.2e[[2]],'CompEntPred_2e_V12_V34.jpg')
compEnt_v34_v56.2e <- CompEntPred.pair(ent.Raes34.2e[[2]],ent.Raes56.2e[[2]],'CompEntPred_2e_V34_V56.jpg')
compEnt_v12_v56.2e <- CompEntPred.pair(ent.Raes12.2e[[2]],ent.Raes56.2e[[2]],'CompEntPred_2e_V12_V56.jpg')

compEnt_v12_v34.3e <- CompEntPred.pair(ent.Raes12.3e[[2]],ent.Raes34.3e[[2]],'CompEntPred_3e_V12_V34.jpg')
compEnt_v34_v56.3e <- CompEntPred.pair(ent.Raes34.3e[[2]],ent.Raes56.3e[[2]],'CompEntPred_3e_V34_V56.jpg')
compEnt_v12_v56.3e <- CompEntPred.pair(ent.Raes12.3e[[2]],ent.Raes56.3e[[2]],'CompEntPred_3e_V12_V56.jpg')

compEnt_v12_v34.4e <- CompEntPred.pair(ent.Raes12.4e[[2]],ent.Raes34.4e[[2]],'CompEntPred_4e_V12_V34.jpg')
compEnt_v34_v56.4e <- CompEntPred.pair(ent.Raes34.4e[[2]],ent.Raes56.4e[[2]],'CompEntPred_4e_V34_V56.jpg')
compEnt_v12_v56.4e <- CompEntPred.pair(ent.Raes12.4e[[2]],ent.Raes56.4e[[2]],'CompEntPred_4e_V12_V56.jpg')

# prepare the final tables
# .number = number of samples the same across amplicons, .percent = percentage of the samples
# v12 vs v56
result_12_56.number <- matrix(0,5,3,dimnames = list(c('All_loc_togeth','IL','RE','TR','All_loc_sep'),c(2,3,4)))
result_12_56.number[1,1] <- as.numeric(compEnt_v12_v56.2e$Maximum_number_samples_clustering_together[,1])
result_12_56.number[2,1] <- as.numeric(compEnt_v12_v56.IL.2e$Maximum_number_samples_clustering_together[,1])
result_12_56.number[3,1] <- as.numeric(compEnt_v12_v56.RE.2e$Maximum_number_samples_clustering_together[,1])
result_12_56.number[4,1] <- as.numeric(compEnt_v12_v56.TR.2e$Maximum_number_samples_clustering_together[,1])
result_12_56.number[1,2] <- as.numeric(compEnt_v12_v56.3e$Maximum_number_samples_clustering_together[,1])
result_12_56.number[2,2] <- as.numeric(compEnt_v12_v56.IL.3e$Maximum_number_samples_clustering_together[,1])
result_12_56.number[3,2] <- as.numeric(compEnt_v12_v56.RE.3e$Maximum_number_samples_clustering_together[,1])
result_12_56.number[4,2] <- as.numeric(compEnt_v12_v56.TR.3e$Maximum_number_samples_clustering_together[,1])
result_12_56.number[1,3] <- as.numeric(compEnt_v12_v56.4e$Maximum_number_samples_clustering_together[,1])
result_12_56.number[2,3] <- as.numeric(compEnt_v12_v56.IL.4e$Maximum_number_samples_clustering_together[,1])
result_12_56.number[3,3] <- as.numeric(compEnt_v12_v56.RE.4e$Maximum_number_samples_clustering_together[,1])
result_12_56.number[4,3] <- as.numeric(compEnt_v12_v56.TR.4e$Maximum_number_samples_clustering_together[,1])
result_12_56.number[5,1] <- result_12_56.number[2,1]+result_12_56.number[3,1]+result_12_56.number[4,1]
result_12_56.number[5,2] <- result_12_56.number[2,2]+result_12_56.number[3,2]+result_12_56.number[4,2]
result_12_56.number[5,3] <- result_12_56.number[2,3]+result_12_56.number[3,3]+result_12_56.number[4,3]

result_12_56.percent <- matrix(0,5,3,dimnames = list(c('All_loc_togeth','IL','RE','TR','All_loc_sep'),c(2,3,4)))
result_12_56.percent[1,1] <- compEnt_v12_v56.2e$Maximum_number_samples_clustering_together[,2]
result_12_56.percent[2,1] <- compEnt_v12_v56.IL.2e$Maximum_number_samples_clustering_together[,2]
result_12_56.percent[3,1] <- compEnt_v12_v56.RE.2e$Maximum_number_samples_clustering_together[,2]
result_12_56.percent[4,1] <- compEnt_v12_v56.TR.2e$Maximum_number_samples_clustering_together[,2]
result_12_56.percent[1,2] <- compEnt_v12_v56.3e$Maximum_number_samples_clustering_together[,2]
result_12_56.percent[2,2] <- compEnt_v12_v56.IL.3e$Maximum_number_samples_clustering_together[,2]
result_12_56.percent[3,2] <- compEnt_v12_v56.RE.3e$Maximum_number_samples_clustering_together[,2]
result_12_56.percent[4,2] <- compEnt_v12_v56.TR.3e$Maximum_number_samples_clustering_together[,2]
result_12_56.percent[1,3] <- compEnt_v12_v56.4e$Maximum_number_samples_clustering_together[,2]
result_12_56.percent[2,3] <- compEnt_v12_v56.IL.4e$Maximum_number_samples_clustering_together[,2]
result_12_56.percent[3,3] <- compEnt_v12_v56.RE.4e$Maximum_number_samples_clustering_together[,2]
result_12_56.percent[4,3] <- compEnt_v12_v56.TR.4e$Maximum_number_samples_clustering_together[,2]
result_12_56.percent[5,1] <- paste(round(100*(result_12_56.number[2,1] + result_12_56.number[3,1] + result_12_56.number[4,1])/as.numeric(compEnt_v12_v56.2e$Maximum_number_samples_clustering_together[,3]),2),'%',sep='')
result_12_56.percent[5,2] <- paste(round(100*(result_12_56.number[2,2] + result_12_56.number[3,2] + result_12_56.number[4,2])/as.numeric(compEnt_v12_v56.2e$Maximum_number_samples_clustering_together[,3]),2),'%',sep='')
result_12_56.percent[5,3] <- paste(round(100*(result_12_56.number[2,3] + result_12_56.number[3,3] + result_12_56.number[4,3])/as.numeric(compEnt_v12_v56.2e$Maximum_number_samples_clustering_together[,3]),2),'%',sep='')

# v12 vs v34
result_12_34.number <- matrix(0,5,3,dimnames = list(c('All_loc_togeth','IL','RE','TR','All_loc_sep'),c(2,3,4)))
result_12_34.number[1,1] <- as.numeric(compEnt_v12_v34.2e$Maximum_number_samples_clustering_together[,1])
result_12_34.number[2,1] <- as.numeric(compEnt_v12_v34.IL.2e$Maximum_number_samples_clustering_together[,1])
result_12_34.number[3,1] <- as.numeric(compEnt_v12_v34.RE.2e$Maximum_number_samples_clustering_together[,1])
result_12_34.number[4,1] <- as.numeric(compEnt_v12_v34.TR.2e$Maximum_number_samples_clustering_together[,1])
result_12_34.number[1,2] <- as.numeric(compEnt_v12_v34.3e$Maximum_number_samples_clustering_together[,1])
result_12_34.number[2,2] <- as.numeric(compEnt_v12_v34.IL.3e$Maximum_number_samples_clustering_together[,1])
result_12_34.number[3,2] <- as.numeric(compEnt_v12_v34.RE.3e$Maximum_number_samples_clustering_together[,1])
result_12_34.number[4,2] <- as.numeric(compEnt_v12_v34.TR.3e$Maximum_number_samples_clustering_together[,1])
result_12_34.number[1,3] <- as.numeric(compEnt_v12_v34.4e$Maximum_number_samples_clustering_together[,1])
result_12_34.number[2,3] <- as.numeric(compEnt_v12_v34.IL.4e$Maximum_number_samples_clustering_together[,1])
result_12_34.number[3,3] <- as.numeric(compEnt_v12_v34.RE.4e$Maximum_number_samples_clustering_together[,1])
result_12_34.number[4,3] <- as.numeric(compEnt_v12_v34.TR.4e$Maximum_number_samples_clustering_together[,1])
result_12_34.number[5,1] <- result_12_34.number[2,1]+result_12_34.number[3,1]+result_12_34.number[4,1]
result_12_34.number[5,2] <- result_12_34.number[2,2]+result_12_34.number[3,2]+result_12_34.number[4,2]
result_12_34.number[5,3] <- result_12_34.number[2,3]+result_12_34.number[3,3]+result_12_34.number[4,3]

result_12_34.percent <- matrix(0,5,3,dimnames = list(c('All_loc_togeth','IL','RE','TR','All_loc_sep'),c(2,3,4)))
result_12_34.percent[1,1] <- compEnt_v12_v34.2e$Maximum_number_samples_clustering_together[,2]
result_12_34.percent[2,1] <- compEnt_v12_v34.IL.2e$Maximum_number_samples_clustering_together[,2]
result_12_34.percent[3,1] <- compEnt_v12_v34.RE.2e$Maximum_number_samples_clustering_together[,2]
result_12_34.percent[4,1] <- compEnt_v12_v34.TR.2e$Maximum_number_samples_clustering_together[,2]
result_12_34.percent[1,2] <- compEnt_v12_v34.3e$Maximum_number_samples_clustering_together[,2]
result_12_34.percent[2,2] <- compEnt_v12_v34.IL.3e$Maximum_number_samples_clustering_together[,2]
result_12_34.percent[3,2] <- compEnt_v12_v34.RE.3e$Maximum_number_samples_clustering_together[,2]
result_12_34.percent[4,2] <- compEnt_v12_v34.TR.3e$Maximum_number_samples_clustering_together[,2]
result_12_34.percent[1,3] <- compEnt_v12_v34.4e$Maximum_number_samples_clustering_together[,2]
result_12_34.percent[2,3] <- compEnt_v12_v34.IL.4e$Maximum_number_samples_clustering_together[,2]
result_12_34.percent[3,3] <- compEnt_v12_v34.RE.4e$Maximum_number_samples_clustering_together[,2]
result_12_34.percent[4,3] <- compEnt_v12_v34.TR.4e$Maximum_number_samples_clustering_together[,2]
result_12_34.percent[5,1] <- paste(round(100*(result_12_34.number[2,1] + result_12_34.number[3,1] + result_12_34.number[4,1])/as.numeric(compEnt_v12_v34.2e$Maximum_number_samples_clustering_together[,3]),2),'%',sep='')
result_12_34.percent[5,2] <- paste(round(100*(result_12_34.number[2,2] + result_12_34.number[3,2] + result_12_34.number[4,2])/as.numeric(compEnt_v12_v34.2e$Maximum_number_samples_clustering_together[,3]),2),'%',sep='')
result_12_34.percent[5,3] <- paste(round(100*(result_12_34.number[2,3] + result_12_34.number[3,3] + result_12_34.number[4,3])/as.numeric(compEnt_v12_v34.2e$Maximum_number_samples_clustering_together[,3]),2),'%',sep='')


# v34 vs v56
result_34_56.number <- matrix(0,5,3,dimnames = list(c('All_loc_togeth','IL','RE','TR','All_loc_sep'),c(2,3,4)))
result_34_56.number[1,1] <- as.numeric(compEnt_v34_v56.2e$Maximum_number_samples_clustering_together[,1])
result_34_56.number[2,1] <- as.numeric(compEnt_v34_v56.IL.2e$Maximum_number_samples_clustering_together[,1])
result_34_56.number[3,1] <- as.numeric(compEnt_v34_v56.RE.2e$Maximum_number_samples_clustering_together[,1])
result_34_56.number[4,1] <- as.numeric(compEnt_v34_v56.TR.2e$Maximum_number_samples_clustering_together[,1])
result_34_56.number[1,2] <- as.numeric(compEnt_v34_v56.3e$Maximum_number_samples_clustering_together[,1])
result_34_56.number[2,2] <- as.numeric(compEnt_v34_v56.IL.3e$Maximum_number_samples_clustering_together[,1])
result_34_56.number[3,2] <- as.numeric(compEnt_v34_v56.RE.3e$Maximum_number_samples_clustering_together[,1])
result_34_56.number[4,2] <- as.numeric(compEnt_v34_v56.TR.3e$Maximum_number_samples_clustering_together[,1])
result_34_56.number[1,3] <- as.numeric(compEnt_v34_v56.4e$Maximum_number_samples_clustering_together[,1])
result_34_56.number[2,3] <- as.numeric(compEnt_v34_v56.IL.4e$Maximum_number_samples_clustering_together[,1])
result_34_56.number[3,3] <- as.numeric(compEnt_v34_v56.RE.4e$Maximum_number_samples_clustering_together[,1])
result_34_56.number[4,3] <- as.numeric(compEnt_v34_v56.TR.4e$Maximum_number_samples_clustering_together[,1])
result_34_56.number[5,1] <- result_34_56.number[2,1]+result_34_56.number[3,1]+result_34_56.number[4,1]
result_34_56.number[5,2] <- result_34_56.number[2,2]+result_34_56.number[3,2]+result_34_56.number[4,2]
result_34_56.number[5,3] <- result_34_56.number[2,3]+result_34_56.number[3,3]+result_34_56.number[4,3]

result_34_56.percent <- matrix(0,5,3,dimnames = list(c('All_loc_togeth','IL','RE','TR','All_loc_sep'),c(2,3,4)))
result_34_56.percent[1,1] <- compEnt_v34_v56.2e$Maximum_number_samples_clustering_together[,2]
result_34_56.percent[2,1] <- compEnt_v34_v56.IL.2e$Maximum_number_samples_clustering_together[,2]
result_34_56.percent[3,1] <- compEnt_v34_v56.RE.2e$Maximum_number_samples_clustering_together[,2]
result_34_56.percent[4,1] <- compEnt_v34_v56.TR.2e$Maximum_number_samples_clustering_together[,2]
result_34_56.percent[1,2] <- compEnt_v34_v56.3e$Maximum_number_samples_clustering_together[,2]
result_34_56.percent[2,2] <- compEnt_v34_v56.IL.3e$Maximum_number_samples_clustering_together[,2]
result_34_56.percent[3,2] <- compEnt_v34_v56.RE.3e$Maximum_number_samples_clustering_together[,2]
result_34_56.percent[4,2] <- compEnt_v34_v56.TR.3e$Maximum_number_samples_clustering_together[,2]
result_34_56.percent[1,3] <- compEnt_v34_v56.4e$Maximum_number_samples_clustering_together[,2]
result_34_56.percent[2,3] <- compEnt_v34_v56.IL.4e$Maximum_number_samples_clustering_together[,2]
result_34_56.percent[3,3] <- compEnt_v34_v56.RE.4e$Maximum_number_samples_clustering_together[,2]
result_34_56.percent[4,3] <- compEnt_v34_v56.TR.4e$Maximum_number_samples_clustering_together[,2]
result_34_56.percent[5,1] <- paste(round(100*(result_34_56.number[2,1] + result_34_56.number[3,1] + result_34_56.number[4,1])/as.numeric(compEnt_v34_v56.2e$Maximum_number_samples_clustering_together[,3]),2),'%',sep='')
result_34_56.percent[5,2] <- paste(round(100*(result_34_56.number[2,2] + result_34_56.number[3,2] + result_34_56.number[4,2])/as.numeric(compEnt_v34_v56.2e$Maximum_number_samples_clustering_together[,3]),2),'%',sep='')
result_34_56.percent[5,3] <- paste(round(100*(result_34_56.number[2,3] + result_34_56.number[3,3] + result_34_56.number[4,3])/as.numeric(compEnt_v34_v56.2e$Maximum_number_samples_clustering_together[,3]),2),'%',sep='')

# now reorganize the tables so that each table corresponds to a number of enterotypes
# 2 enterotypes
result_2e.number <- cbind(result_12_34.number[,1],result_12_56.number[,1],result_34_56.number[,1])
colnames(result_2e.number) <- c('v12_v34','v12_v56','v34_v56')
result_2e.percent <- cbind(result_12_34.percent[,1],result_12_56.percent[,1],result_34_56.percent[,1])
colnames(result_2e.percent) <- c('v12_v34','v12_v56','v34_v56')

# 3 enterotypes
result_3e.number <- cbind(result_12_34.number[,2],result_12_56.number[,2],result_34_56.number[,2])
colnames(result_3e.number) <- c('v12_v34','v12_v56','v34_v56')
result_3e.percent <- cbind(result_12_34.percent[,2],result_12_56.percent[,2],result_34_56.percent[,2])
colnames(result_3e.percent) <- c('v12_v34','v12_v56','v34_v56')

# 4 enterotypes
result_4e.number <- cbind(result_12_34.number[,3],result_12_56.number[,3],result_34_56.number[,3])
colnames(result_4e.number) <- c('v12_v34','v12_v56','v34_v56')
result_4e.percent <- cbind(result_12_34.percent[,3],result_12_56.percent[,3],result_34_56.percent[,3])
colnames(result_4e.percent) <- c('v12_v34','v12_v56','v34_v56')


setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Compare_enterotypings')
write.table(result_12_34.number,'v12_v34_ent_comparison_incl_locat_ent_separ_number.csv')
write.table(result_12_34.percent,'v12_v34_ent_comparison_incl_locat_ent_separ_percent.csv')
write.table(result_12_56.number,'v12_v56_ent_comparison_incl_locat_ent_separ_number.csv')
write.table(result_12_56.percent,'v12_v56_ent_comparison_incl_locat_ent_separ_percent.csv')
write.table(result_34_56.number,'v34_v56_ent_comparison_incl_locat_ent_separ_number.csv')
write.table(result_34_56.percent,'v34_v56_ent_comparison_incl_locat_ent_separ_percent.csv')

