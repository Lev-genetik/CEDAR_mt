# Compare intestine locations
# 1. plot IL, TR and RE locations on first 2 PCs of the PCA analysis
source(file.path(progr.dir,'CompareLocations_or_amplicons.R'))
plotIntestineLocations(qiime_out = file.path(home.dir,'Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt'),
                                             outdir=file.path(home.dir,'Results/Locations_compare/Graphs'),
                                             sample.coverage.filter=23000)
plotIntestineLocations(qiime_out = file.path(home.dir,'Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt'),
                       outdir=file.path(home.dir,'Results/Locations_compare/Graphs'),
                       sample.coverage.filter=23000,draw.PCoA.sep=seq(10,18), draw.PCA.sep=seq(10,18))

plotIntestineLocations(qiime_out = file.path(home.dir,'Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt'),
                       outdir=file.path(home.dir,'Results/Locations_compare/Graphs'),
                       sample.coverage.filter=13000)
plotIntestineLocations(qiime_out = file.path(home.dir,'Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt'),
                       outdir=file.path(home.dir,'Results/Locations_compare/Graphs'),
                       sample.coverage.filter=13000,draw.PCoA.sep=seq(10,18), draw.PCA.sep=seq(10,18))
plotIntestineLocations(qiime_out = file.path(home.dir,'Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt'),
                       outdir=file.path(home.dir,'Results/Locations_compare/Graphs'),
                       sample.coverage.filter=21000)
plotIntestineLocations(qiime_out = file.path(home.dir,'Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt'),
                       outdir=file.path(home.dir,'Results/Locations_compare/Graphs'),
                       sample.coverage.filter=21000,draw.PCoA.sep=seq(10,18), draw.PCA.sep=seq(10,18))

# Now, draw pairwise plots for the PCs that  show difference across locations
# V1V2
plotIntestineLocations(qiime_out = file.path(home.dir,'Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt'),
                       outdir=file.path(home.dir,'Results/Locations_compare/Graphs'),
                       sample.coverage.filter=23000,PCoA.individual=FALSE, PCA.individual=FALSE, 
                       draw.PCoA.pairwise = c(1,2), draw.PCA.pairwise = c(1,2))

plotIntestineLocations(qiime_out = file.path(home.dir,'Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt'),
                       outdir=file.path(home.dir,'Results/Locations_compare/Graphs'),
                       sample.coverage.filter=23000,PCoA.individual=FALSE, PCA.individual=FALSE, 
                       draw.PCoA.pairwise = c(6,7), draw.PCA.pairwise = c(6,18))
plotIntestineLocations(qiime_out = file.path(home.dir,'Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt'),
                       outdir=file.path(home.dir,'Results/Locations_compare/Graphs'),
                       sample.coverage.filter=23000,PCoA.individual=FALSE, PCA.individual=FALSE, PCA.pair=FALSE,
                       draw.PCoA.pairwise = c(8,9))
plotIntestineLocations(qiime_out = file.path(home.dir,'Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt'),
                       outdir=file.path(home.dir,'Results/Locations_compare/Graphs'),
                       sample.coverage.filter=23000,PCoA.individual=FALSE, PCA.individual=FALSE, PCA.pair=FALSE,
                       draw.PCoA.pairwise = c(6,10))
# V3V4
plotIntestineLocations(qiime_out = file.path(home.dir,'Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt'),
                       outdir=file.path(home.dir,'Results/Locations_compare/Graphs'),
                       sample.coverage.filter=13000,PCoA.individual=FALSE, PCA.individual=FALSE, 
                       draw.PCoA.pairwise = c(1,2), draw.PCA.pairwise = c(1,2))
plotIntestineLocations(qiime_out = file.path(home.dir,'Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt'),
                       outdir=file.path(home.dir,'Results/Locations_compare/Graphs'),
                       sample.coverage.filter=13000,PCoA.individual=FALSE, PCA.individual=FALSE, 
                       draw.PCoA.pairwise = c(1,9), draw.PCA.pairwise = c(1,2))
plotIntestineLocations(qiime_out = file.path(home.dir,'Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt'),
                       outdir=file.path(home.dir,'Results/Locations_compare/Graphs'),
                       sample.coverage.filter=13000,PCoA.individual=FALSE, PCA.individual=FALSE, 
                       draw.PCoA.pairwise = c(9,10), draw.PCA.pairwise = c(3,4))
plotIntestineLocations(qiime_out = file.path(home.dir,'Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt'),
                       outdir=file.path(home.dir,'Results/Locations_compare/Graphs'),
                       sample.coverage.filter=13000,PCoA.individual=FALSE, PCA.individual=FALSE, 
                       draw.PCoA.pairwise = c(9,10), draw.PCA.pairwise = c(3,8))

# V5V6
plotIntestineLocations(qiime_out = file.path(home.dir,'Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt'),
                       outdir=file.path(home.dir,'Results/Locations_compare/Graphs'),
                       sample.coverage.filter=21000,PCoA.individual=FALSE, PCA.individual=FALSE, 
                       draw.PCoA.pairwise = c(1,2), draw.PCA.pairwise = c(1,2))

plotIntestineLocations(qiime_out = file.path(home.dir,'Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt'),
                       outdir=file.path(home.dir,'Results/Locations_compare/Graphs'),
                       sample.coverage.filter=21000,PCoA.individual=FALSE, PCA.individual=FALSE, 
                       draw.PCoA.pairwise = c(2,6), draw.PCA.pairwise = c(2,4))
plotIntestineLocations(qiime_out = file.path(home.dir,'Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt'),
                       outdir=file.path(home.dir,'Results/Locations_compare/Graphs'),
                       sample.coverage.filter=21000,PCoA.individual=FALSE, PCA.individual=FALSE, 
                       draw.PCoA.pairwise = c(6,7), draw.PCA.pairwise = c(5,13))


# 2. BoxPlots: Plot abundance of bacterial families across locations
source(file.path(progr.dir,'CompareLocations_or_amplicons.R'))
compTaxaBoxplots(file1 = file.path(home.dir,'Summarize_taxa_counts_original/summarize_taxa_counts_v12/otu_table_L5_v2.txt'),
                 file2 = file.path(home.dir,'Summarize_taxa_counts_original/summarize_taxa_counts_v12/otu_table_L5_v2.txt'),
                 file3 = file.path(home.dir,'Summarize_taxa_counts_original/summarize_taxa_counts_v34/otu_table_L5_v2.txt'),
                 loc1 = '.IL',loc2 = '.TR', loc3 = '.RE', 
                 file_out = file.path(home.dir,'Results/Locations_compare/Abundance_boxplots/V1V2_cov23k_outlRem_IL_TR_RE_comp_families.jpg'), 
                 sample.cov.filter = 23000,percent_threshold=0,min.abundance.any.sample = 0.001)
compTaxaBoxplots(file1 = file.path(home.dir,'Summarize_taxa_counts_original/summarize_taxa_counts_v12/otu_table_L5_v2.txt'),
                 file2 = file.path(home.dir,'Summarize_taxa_counts_original/summarize_taxa_counts_v12/otu_table_L5_v2.txt'),
                 file3 = file.path(home.dir,'Summarize_taxa_counts_original/summarize_taxa_counts_v34/otu_table_L5_v2.txt'),
                 loc1 = '.IL',loc2 = '.TR', loc3 = '.RE',taxon.abundance.threshold.plot = 0,
                 file_out = file.path(home.dir,'Results/Locations_compare/Abundance_boxplots/V1V2_cov23k_outlRem_IL_TR_RE_comp_families_abund_thr_0.jpg'),
                 sample.cov.filter = 23000,percent_threshold=0,min.abundance.any.sample = 0)
compTaxaBoxplots(file1 = file.path(home.dir,'Summarize_taxa_counts_original/summarize_taxa_counts_v34/otu_table_L5_v2.txt'),
                 file2 = file.path(home.dir,'Summarize_taxa_counts_original/summarize_taxa_counts_v34/otu_table_L5_v2.txt'),
                 file3 = file.path(home.dir,'Summarize_taxa_counts_original/summarize_taxa_counts_v34/otu_table_L5_v2.txt'),
                 loc1 = '.IL',loc2 = '.TR', loc3 = '.RE', 
                 file_out = file.path(home.dir,'Results/Locations_compare/Abundance_boxplots/V3V4_cov13k_outlRem_IL_TR_RE_comp_families.jpg'), 
                 sample.cov.filter = 13000,percent_threshold=0,min.abundance.any.sample = 0.001)
# compTaxaBoxplots(file1 = file.path(home.dir,'Summarize_taxa_counts_original/summarize_taxa_counts_v34/otu_table_L5_v2.txt'),
#                  file2 = file.path(home.dir,'Summarize_taxa_counts_original/summarize_taxa_counts_v34/otu_table_L5_v2.txt'),
#                  file3 = file.path(home.dir,'Summarize_taxa_counts_original/summarize_taxa_counts_v34/otu_table_L5_v2.txt'),
#                  loc1 = '.IL',loc2 = '.TR', loc3 = '.RE', 
#                  file_out = file.path(home.dir,'Results/Locations_compare/Abundance_boxplots/V3V4_cov13k_outlRem_IL_TR_RE_comp_families_abund_thr_0_00001.jpg'), 
#                  sample.cov.filter = 13000,percent_threshold=0,min.abundance.any.sample = 0.0001)
compTaxaBoxplots(file1 = file.path(home.dir,'Summarize_taxa_counts_original/summarize_taxa_counts_v56/otu_table_L5_v2.txt'),
                 file2 = file.path(home.dir,'Summarize_taxa_counts_original/summarize_taxa_counts_v56/otu_table_L5_v2.txt'),
                 file3 = file.path(home.dir,'Summarize_taxa_counts_original/summarize_taxa_counts_v56/otu_table_L5_v2.txt'),
                 loc1 = '.IL',loc2 = '.TR', loc3 = '.RE', 
                 file_out = file.path(home.dir,'Results/Locations_compare/Abundance_boxplots/V5V6_cov21k_outlRem_IL_TR_RE_comp_families.jpg'), 
                 sample.cov.filter = 21000,percent_threshold=0,min.abundance.any.sample = 0.001)
# compTaxaBoxplots(file1 = file.path(home.dir,'Summarize_taxa_counts_original/summarize_taxa_counts_v56/otu_table_L5_v2.txt'),
#                  file2 = file.path(home.dir,'Summarize_taxa_counts_original/summarize_taxa_counts_v56/otu_table_L5_v2.txt'),
#                  file3 = file.path(home.dir,'Summarize_taxa_counts_original/summarize_taxa_counts_v56/otu_table_L5_v2.txt'),
#                  loc1 = '.IL',loc2 = '.TR', loc3 = '.RE', 
#                  file_out = file.path(home.dir,'Results/Locations_compare/Abundance_boxplots/V5V6_cov21k_outlRem_IL_TR_RE_comp_families_abund_thr_0_00001.jpg'), 
#                  sample.cov.filter = 21000,percent_threshold=0,min.abundance.any.sample = 0.0001)