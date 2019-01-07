# See the PCs that are associated with location difference - what genera do they correspond to?

# Make enterotyping for all locations V5V6
qiime_out = file.path(home.dir,'Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt')
sample.coverage.filter=21000
source(file.path(progr.dir,'Enterotyping_Raes_noise_removed_v12.R'))
tmp.dir <- file.path(home.dir,'tmp')
dir.create(tmp.dir,showWarnings = FALSE)
ent.Raes.2e.v5v6.all.locat <- Enterotyping_Raes(qiime_out, output_dir = tmp.dir, clust.number=2, noiserem=T, plots=T,
                                            sample.coverage.filter = sample.coverage.filter,
                                            remove.outliers = 0.01*outliers.removed,locat = 'all')
dim(ent.Raes.2e.v5v6.all.locat$obs_pca_result$c1)
eigen.vector1 <- ent.Raes.2e.v5v6.all.locat$obs_pca_result$c1[,1]
names(eigen.vector1) <- rownames(ent.Raes.2e.v5v6.all.locat$obs_pca_result$c1)
eigen.vector1[order(eigen.vector1,decreasing = TRUE)][1:5]
eigen.vector2 <- ent.Raes.2e.v5v6.all.locat$obs_pca_result$c1[,2]
names(eigen.vector2) <- rownames(ent.Raes.2e.v5v6.all.locat$obs_pca_result$c1)
eigen.vector2[order(eigen.vector2,decreasing = TRUE)][1:10]
# Enriched in "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria": 5/10 here and 20/97 in the whole table
eigen.vector3 <- ent.Raes.2e.v5v6.all.locat$obs_pca_result$c1[,3]
names(eigen.vector3) <- rownames(ent.Raes.2e.v5v6.all.locat$obs_pca_result$c1)
eigen.vector3[order(eigen.vector3,decreasing = TRUE)][1:5]
# This is most intr
eigen.vector4 <- ent.Raes.2e.v5v6.all.locat$obs_pca_result$c1[,4]
names(eigen.vector4) <- rownames(ent.Raes.2e.v5v6.all.locat$obs_pca_result$c1)
eigen.vector4[order(eigen.vector4,decreasing = TRUE)][1:10]
# The 4th eigenvector is enriched by f__Pasteurellaceae: 1,2,3 and 7th largest coordinate, no other 
# f__Pasteurellaceae among the 97 taxa!!!
rownames(ent.Raes.2e.v5v6.all.locat$obs_pca_result$c1)[rownames(ent.Raes.2e.v5v6.all.locat$obs_pca_result$c1) %like% 
                                                    'k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae']
# [1] "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.g__"               
# [2] "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.g__Actinobacillus" 
# [3] "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.g__Aggregatibacter"
# [4] "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.g__Haemophilus"    




# Make enterotyping for all locations v3v4
# PC2 is the best separator
qiime_out = file.path(home.dir,'Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt')
sample.coverage.filter=13000
source(file.path(progr.dir,'Enterotyping_Raes_noise_removed_v12.R'))
tmp.dir <- file.path(home.dir,'tmp')
dir.create(tmp.dir,showWarnings = FALSE)
ent.Raes.2e.v3v4.all.locat <- Enterotyping_Raes(qiime_out, output_dir = tmp.dir, clust.number=2, noiserem=T, plots=T,
                                                sample.coverage.filter = sample.coverage.filter,
                                                remove.outliers = 0.01*outliers.removed,locat = 'all')
dim(ent.Raes.2e.v3v4.all.locat$obs_pca_result$c1)
# eigen.v3v4.vector1 <- ent.Raes.2e.v3v4.all.locat$obs_pca_result$c1[,1]
# names(eigen.v3v4.vector1) <- rownames(ent.Raes.2e.v3v4.all.locat$obs_pca_result$c1)
# eigen.v3v4.vector1[order(eigen.v3v4.vector1,decreasing = TRUE)][1:10]
# do not know what it is

# This is the most intr
eigen.v3v4.vector2 <- ent.Raes.2e.v3v4.all.locat$obs_pca_result$c1[,2]
names(eigen.v3v4.vector2) <- rownames(ent.Raes.2e.v3v4.all.locat$obs_pca_result$c1)
eigen.v3v4.vector2[order(eigen.v3v4.vector2,decreasing = TRUE)][1:10]
# Here we again see 
# Enriched in "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria": 4/4 leader species (1-3 have bif difference
# with others).

eigen.v3v4.vector3 <- ent.Raes.2e.v3v4.all.locat$obs_pca_result$c1[,3]
names(eigen.v3v4.vector3) <- rownames(ent.Raes.2e.v3v4.all.locat$obs_pca_result$c1)
eigen.v3v4.vector3[order(eigen.v3v4.vector3,decreasing = TRUE)][1:10]

eigen.v3v4.vector4 <- ent.Raes.2e.v3v4.all.locat$obs_pca_result$c1[,4]
names(eigen.v3v4.vector4) <- rownames(ent.Raes.2e.v3v4.all.locat$obs_pca_result$c1)
eigen.v3v4.vector4[order(eigen.v3v4.vector4,decreasing = TRUE)][1:10]

# For v1V2




# Get correlations between PCs
