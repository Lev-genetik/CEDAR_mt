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
eigen.v5v6.vector2 <- ent.Raes.2e.v5v6.all.locat$obs_pca_result$c1[,2]
names(eigen.v5v6.vector2) <- rownames(ent.Raes.2e.v5v6.all.locat$obs_pca_result$c1)
eigen.v5v6.vector2[order(eigen.v5v6.vector2,decreasing = TRUE)][1:10]
# Enriched in "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria": 5/10 here and 20/97 in the whole table
# This is most intr #4
eigen.v5v6.vector4 <- ent.Raes.2e.v5v6.all.locat$obs_pca_result$c1[,4]
names(eigen.v5v6.vector4) <- rownames(ent.Raes.2e.v5v6.all.locat$obs_pca_result$c1)
eigen.v5v6.vector4[order(eigen.v5v6.vector4,decreasing = TRUE)][1:10]
# The 4th eigenvector is enriched by f__Pasteurellaceae: 1,2,3 and 7th largest coordinate, no other 
# f__Pasteurellaceae among the 97 taxa!!!
rownames(ent.Raes.2e.v5v6.all.locat$obs_pca_result$c1)[rownames(ent.Raes.2e.v5v6.all.locat$obs_pca_result$c1) %like% 
                                                    'k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae']
# [1] "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.g__"               
# [2] "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.g__Actinobacillus" 
# [3] "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.g__Aggregatibacter"
# [4] "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.g__Haemophilus"    

# 5
eigen.v5v6.vector5 <- ent.Raes.2e.v5v6.all.locat$obs_pca_result$c1[,5]
names(eigen.v5v6.vector5) <- rownames(ent.Raes.2e.v5v6.all.locat$obs_pca_result$c1)
# eigen.v5v6.vector5[order(eigen.v5v6.vector5,decreasing = TRUE)][1:10]
# 13
eigen.v5v6.vector13 <- ent.Raes.2e.v5v6.all.locat$obs_pca_result$c1[,13]
names(eigen.v5v6.vector13) <- rownames(ent.Raes.2e.v5v6.all.locat$obs_pca_result$c1)
# eigen.v5v6.vector13[order(eigen.v5v6.vector13,decreasing = TRUE)][1:10]

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
eigen.v3v4.vector1 <- ent.Raes.2e.v3v4.all.locat$obs_pca_result$c1[,1]
names(eigen.v3v4.vector1) <- rownames(ent.Raes.2e.v3v4.all.locat$obs_pca_result$c1)
eigen.v3v4.vector1[order(eigen.v3v4.vector1,decreasing = TRUE)][1:10]

# 2. This is the most intr
eigen.v3v4.vector2 <- ent.Raes.2e.v3v4.all.locat$obs_pca_result$c1[,2]
names(eigen.v3v4.vector2) <- rownames(ent.Raes.2e.v3v4.all.locat$obs_pca_result$c1)
eigen.v3v4.vector2[order(eigen.v3v4.vector2,decreasing = TRUE)][1:10]
# Here we again see 
# Enriched in "k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria": 4/4 leader species (1-3 have bif difference
# with others).
# 3
eigen.v3v4.vector3 <- ent.Raes.2e.v3v4.all.locat$obs_pca_result$c1[,3]
names(eigen.v3v4.vector3) <- rownames(ent.Raes.2e.v3v4.all.locat$obs_pca_result$c1)
# eigen.v3v4.vector3[order(eigen.v3v4.vector3,decreasing = TRUE)][1:10]
# 4.
eigen.v3v4.vector4 <- ent.Raes.2e.v3v4.all.locat$obs_pca_result$c1[,4]
names(eigen.v3v4.vector4) <- rownames(ent.Raes.2e.v3v4.all.locat$obs_pca_result$c1)
# eigen.v3v4.vector4[order(eigen.v3v4.vector4,decreasing = TRUE)][1:10]
# 8.
eigen.v3v4.vector8 <- ent.Raes.2e.v3v4.all.locat$obs_pca_result$c1[,8]
names(eigen.v3v4.vector8) <- rownames(ent.Raes.2e.v3v4.all.locat$obs_pca_result$c1)
# eigen.v3v4.vector8[order(eigen.v3v4.vector8,decreasing = TRUE)][1:10]

# For v1V2
sample.coverage.filter=23000
qiime_out = file.path(home.dir,'Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt')
ent.Raes.2e.v1v2.all.locat <- Enterotyping_Raes(qiime_out, output_dir = tmp.dir, clust.number=2, noiserem=T, plots=T,
                                                sample.coverage.filter = sample.coverage.filter,
                                                remove.outliers = 0.01*outliers.removed,locat = 'all')
eigen.v1v2.vector6 <- ent.Raes.2e.v1v2.all.locat$obs_pca_result$c1[,6]
names(eigen.v1v2.vector6) <- rownames(ent.Raes.2e.v1v2.all.locat$obs_pca_result$c1)
eigen.v1v2.vector18 <- ent.Raes.2e.v1v2.all.locat$obs_pca_result$c1[,18]
names(eigen.v1v2.vector18) <- rownames(ent.Raes.2e.v1v2.all.locat$obs_pca_result$c1)



# Get correlations between PCs
common.taxa <- intersect(intersect(rownames(ent.Raes.2e.v1v2.all.locat$obs_pca_result$c1),
                                   rownames(ent.Raes.2e.v3v4.all.locat$obs_pca_result$c1)),
                         rownames(ent.Raes.2e.v5v6.all.locat$obs_pca_result$c1))
# prepare shortened PCs on the coordinates of taxa (only including the taxa common across ampl.)
v1v2.6 <- eigen.v1v2.vector6[common.taxa]
v1v2.18 <- eigen.v1v2.vector18[common.taxa]

v3v4.1 <- eigen.v3v4.vector1[common.taxa]
v3v4.2 <- eigen.v3v4.vector2[common.taxa]
v3v4.3 <- eigen.v3v4.vector3[common.taxa]
v3v4.4 <- eigen.v3v4.vector4[common.taxa]
v3v4.8 <- eigen.v3v4.vector8[common.taxa]

v5v6.2 <- eigen.v5v6.vector2[common.taxa]
v5v6.4 <- eigen.v5v6.vector4[common.taxa]
v5v6.5 <- eigen.v5v6.vector5[common.taxa]
v5v6.13 <- eigen.v5v6.vector13[common.taxa]

PCs.df <- cbind(v1v2.6,v1v2.18,v3v4.1,v3v4.2,v3v4.3,v3v4.4,v3v4.8,v5v6.2,v5v6.4,v5v6.5,v5v6.13)
PCA.cor.table <- round(cor(PCs.df),2)
setwd(file.path(home.dir,'Results/Locations_compare'))
write.table(PCA.cor.table,'PCA.important.PCs.cor.table.csv',quote=F,sep='\t')
# The best correlation
pair.1 <- as.data.frame(cbind(v3v4.1,v5v6.2))
regr.1 <- lm(v3v4.1~v5v6.2,data=pair.1)
good.taxa <- pair.1[abs(pair.1[,1]-pair.1[,2]*regr.1$coefficients[2]-regr.1$coefficients[1])<0.1,]
good.taxa <- good.taxa[order(good.taxa[,1],decreasing = T),]
plot(pair.1)
good.taxa[rownames(good.taxa)%like% 'k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria',]
good.taxa[rownames(good.taxa)%like% 'k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.',]
# We can speculate that TR is characterized by high presence of class Gammaproteobacteria and low presence of  o__Clostridiales