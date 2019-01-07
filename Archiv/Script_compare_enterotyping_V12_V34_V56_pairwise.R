source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Enterotyping_Raes_noise_removed_v8.R')
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Enterotyping_compar_functions_v3.R')
# first, make enterotyping
# 2 enterotypes
ent.Raes12.2e <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Clustering/V1V2', clust.number=2, noiserem=T, plots=F)
ent.Raes34.2e <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Clustering/V3V4', clust.number=2, noiserem=T, plots=F)
ent.Raes56.2e <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Clustering/V5V6', clust.number=2, noiserem=T, plots=F)
# 3 enterotypes
ent.Raes12.3e <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Clustering/V1V2', clust.number=3, noiserem=T, plots=F)
ent.Raes34.3e <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Clustering/V3V4', clust.number=3, noiserem=T, plots=F)
ent.Raes56.3e <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Clustering/V5V6', clust.number=3, noiserem=T, plots=F)
# 4 enterotypes
ent.Raes12.4e <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Clustering/V1V2', clust.number=4, noiserem=T, plots=F)
ent.Raes34.4e <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Clustering/V3V4', clust.number=4, noiserem=T, plots=F)
ent.Raes56.4e <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Clustering/V5V6', clust.number=4, noiserem=T, plots=F)
# 5 enterotypes
ent.Raes12.5e <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Clustering/V1V2', clust.number=5, noiserem=T, plots=F)
ent.Raes34.5e <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Clustering/V3V4', clust.number=5, noiserem=T, plots=F)
ent.Raes56.5e <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Clustering/V5V6', clust.number=5, noiserem=T, plots=F)

#now get number of common elements between v12-v34-v56 clustering strategies and enterotype matches
compEnt_v12_v34.2e <- CompEntPred.pair(ent.Raes12.2e[[2]],ent.Raes34.2e[[2]],as.vector(ent.Raes12.2e[[1]][,1]),as.vector(ent.Raes34.2e[[1]][,1]),'CompEntPred_2e_V12_V34.jpg')
compEnt_v34_v56.2e <- CompEntPred.pair(ent.Raes34.2e[[2]],ent.Raes56.2e[[2]],as.vector(ent.Raes34.2e[[1]][,1]),as.vector(ent.Raes56.2e[[1]][,1]),'CompEntPred_2e_V34_V56.jpg')
compEnt_v12_v56.2e <- CompEntPred.pair(ent.Raes12.2e[[2]],ent.Raes56.2e[[2]],as.vector(ent.Raes12.2e[[1]][,1]),as.vector(ent.Raes56.2e[[1]][,1]),'CompEntPred_2e_V12_V56.jpg')

compEnt_v12_v34.3e <- CompEntPred.pair(ent.Raes12.3e[[2]],ent.Raes34.3e[[2]],as.vector(ent.Raes12.3e[[1]][,1]),as.vector(ent.Raes34.3e[[1]][,1]),'CompEntPred_3e_V12_V34.jpg')
compEnt_v34_v56.3e <- CompEntPred.pair(ent.Raes34.3e[[2]],ent.Raes56.3e[[2]],as.vector(ent.Raes34.3e[[1]][,1]),as.vector(ent.Raes56.3e[[1]][,1]),'CompEntPred_3e_V34_V56.jpg')
compEnt_v12_v56.3e <- CompEntPred.pair(ent.Raes12.3e[[2]],ent.Raes56.3e[[2]],as.vector(ent.Raes12.3e[[1]][,1]),as.vector(ent.Raes56.3e[[1]][,1]),'CompEntPred_3e_V12_V56.jpg')

compEnt_v12_v34.4e <- CompEntPred.pair(ent.Raes12.4e[[2]],ent.Raes34.4e[[2]],as.vector(ent.Raes12.4e[[1]][,1]),as.vector(ent.Raes34.4e[[1]][,1]),'CompEntPred_4e_V12_V34.jpg')
compEnt_v34_v56.4e <- CompEntPred.pair(ent.Raes34.4e[[2]],ent.Raes56.4e[[2]],as.vector(ent.Raes34.4e[[1]][,1]),as.vector(ent.Raes56.4e[[1]][,1]),'CompEntPred_4e_V34_V56.jpg')
compEnt_v12_v56.4e <- CompEntPred.pair(ent.Raes12.4e[[2]],ent.Raes56.4e[[2]],as.vector(ent.Raes12.4e[[1]][,1]),as.vector(ent.Raes56.4e[[1]][,1]),'CompEntPred_4e_V12_V56.jpg')

# 5 enterotypes - problem solved! we added second and etc. drivers where required
compEnt_v12_v34.5e <- CompEntPred.pair(ent.Raes12.5e$sample.enterotype,ent.Raes34.5e$sample.enterotype,ent.Raes12.5e$drivers_unique,ent.Raes34.5e$drivers_unique,'CompEntPred_5e_V12_V34.jpg')
compEnt_v34_v56.5e <- CompEntPred.pair(ent.Raes34.5e$sample.enterotype,ent.Raes56.5e$sample.enterotype,ent.Raes34.5e$drivers_unique,ent.Raes56.5e$drivers_unique,'CompEntPred_5e_V34_V56.jpg')
compEnt_v12_v56.5e <- CompEntPred.pair(ent.Raes12.5e$sample.enterotype,ent.Raes56.5e$sample.enterotype,ent.Raes12.5e$drivers_unique,ent.Raes56.5e$drivers_unique,'CompEntPred_5e_V12_V56.jpg')

# now for 4 enterotypes for ileum only
ent.Raes12.IL.4e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Clustering/V1V2', clust.number=4, locat='IL', noiserem=T, plots=F)
ent.Raes34.IL.4e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Clustering/V1V2', clust.number=4, locat='IL', noiserem=T, plots=F)
ent.Raes56.IL.4e  <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Clustering/V1V2', clust.number=4, locat='IL', noiserem=T, plots=F)

compEnt_v12_v34.IL.4e <- CompEntPred.pair(ent.Raes12.IL.4e[[2]],ent.Raes34.IL.4e[[2]],ent.Raes12.IL.4e$drivers_unique,ent.Raes34.IL.4e$drivers_unique,'CompEntPred_4e_V12_V34_IL.jpg')
compEnt_v34_v56.IL.4e <- CompEntPred.pair(ent.Raes34.IL.4e[[2]],ent.Raes56.IL.4e[[2]],ent.Raes34.IL.4e$drivers_unique,ent.Raes56.IL.4e$drivers_unique,'CompEntPred_4e_V34_V56_IL.jpg')
compEnt_v12_v56.IL.4e <- CompEntPred.pair(ent.Raes12.IL.4e[[2]],ent.Raes56.IL.4e[[2]],ent.Raes12.IL.4e$drivers_unique,ent.Raes56.IL.4e$drivers_unique,'CompEntPred_4e_V12_V56_IL.jpg')
