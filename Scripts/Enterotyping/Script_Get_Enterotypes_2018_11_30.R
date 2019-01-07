# Perform enterotyping: 
# each location separately, filter out samples with < threshold coverage (13-25k depending on ampl-loc), filter rare taxa <1%, 
# remove outliers <1% p-value of distance to other samples,
# 2-5 enterotypes,
source(file.path(progr.dir,'/Enterotyping_Raes_noise_removed_v12.R'))
setwd(file.path(home.dir,'Results/Enterotyping/2018_11_30_enterotyping/Sample-enterotype'))
output_dir = file.path(home.dir,'Results/Enterotyping/2018_11_30_enterotyping/Graphs')
input_file_v12 <- file.path(home.dir,'Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt')
input_file_v34 <- file.path(home.dir,'Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt')
input_file_v56 <- file.path(home.dir,'Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt')

coverage.thersholds <- rep(0,times=9)
names(coverage.thersholds) <- c('v12.IL','v12.TR','v12.RE','v34.IL','v34.TR','v34.RE','v56.IL','v56.TR','v56.RE')
coverage.thersholds['v12.IL'] <- 17000
coverage.thersholds['v12.TR'] <- 24000
coverage.thersholds['v12.RE'] <- 22000
coverage.thersholds['v34.IL'] <- 12000
coverage.thersholds['v34.TR'] <- 15000
coverage.thersholds['v34.RE'] <- 13000
coverage.thersholds['v56.IL'] <- 25000
coverage.thersholds['v56.TR'] <- 17000
coverage.thersholds['v56.RE'] <- 12000

# 2 enterotypes
ent.Raes12.2e.sep.by.loc.IL <- Enterotyping_Raes(input_file_v12, output_dir=output_dir, clust.number=2, noiserem=T, remove.outliers = 0.01, locat='.IL', sample.coverage.filter=coverage.thersholds['v12.IL'], plots=T)
ent.Raes12.2e.sep.by.loc.TR <- Enterotyping_Raes(input_file_v12, output_dir=output_dir, clust.number=2, noiserem=T, remove.outliers = 0.01, locat='.TR', sample.coverage.filter=coverage.thersholds['v12.TR'], plots=T)
ent.Raes12.2e.sep.by.loc.RE <- Enterotyping_Raes(input_file_v12, output_dir=output_dir, clust.number=2, noiserem=T, remove.outliers = 0.01, locat='.RE', sample.coverage.filter=coverage.thersholds['v12.RE'], plots=T)
ent.Raes34.2e.sep.by.loc.IL <- Enterotyping_Raes(input_file_v34, output_dir=output_dir, clust.number=2, noiserem=T, remove.outliers = 0.01, locat='.IL', sample.coverage.filter=coverage.thersholds['v34.IL'], plots=T)
ent.Raes34.2e.sep.by.loc.TR <- Enterotyping_Raes(input_file_v34, output_dir=output_dir, clust.number=2, noiserem=T, remove.outliers = 0.01, locat='.TR', sample.coverage.filter=coverage.thersholds['v34.TR'], plots=T)
ent.Raes34.2e.sep.by.loc.RE <- Enterotyping_Raes(input_file_v34, output_dir=output_dir, clust.number=2, noiserem=T, remove.outliers = 0.01, locat='.RE', sample.coverage.filter=coverage.thersholds['v34.RE'], plots=T)
ent.Raes56.2e.sep.by.loc.IL <- Enterotyping_Raes(input_file_v56, output_dir=output_dir, clust.number=2, noiserem=T, remove.outliers = 0.01, locat='.IL', sample.coverage.filter=coverage.thersholds['v56.IL'], plots=T)
ent.Raes56.2e.sep.by.loc.TR <- Enterotyping_Raes(input_file_v56, output_dir=output_dir, clust.number=2, noiserem=T, remove.outliers = 0.01, locat='.TR', sample.coverage.filter=coverage.thersholds['v56.TR'], plots=T)
ent.Raes56.2e.sep.by.loc.RE <- Enterotyping_Raes(input_file_v56, output_dir=output_dir, clust.number=2, noiserem=T, remove.outliers = 0.01, locat='.RE', sample.coverage.filter=coverage.thersholds['v56.RE'], plots=T)

setwd(file.path(home.dir,'Results/Enterotyping/2018_11_30_enterotyping/Sample-enterotype'))
write.table(ent.Raes12.2e.sep.by.loc.IL$sample.enterotype, 'V1V2_2e_IL_sample-enterotype.csv')
write.table(ent.Raes12.2e.sep.by.loc.TR$sample.enterotype, 'V1V2_2e_TR_sample-enterotype.csv')
write.table(ent.Raes12.2e.sep.by.loc.RE$sample.enterotype, 'V1V2_2e_RE_sample-enterotype.csv')

write.table(ent.Raes34.2e.sep.by.loc.IL$sample.enterotype, 'V3V4_2e_IL_sample-enterotype.csv')
write.table(ent.Raes34.2e.sep.by.loc.TR$sample.enterotype, 'V3V4_2e_TR_sample-enterotype.csv')
write.table(ent.Raes34.2e.sep.by.loc.RE$sample.enterotype, 'V3V4_2e_RE_sample-enterotype.csv')

write.table(ent.Raes56.2e.sep.by.loc.IL$sample.enterotype, 'V5V6_2e_IL_sample-enterotype.csv')
write.table(ent.Raes56.2e.sep.by.loc.TR$sample.enterotype, 'V5V6_2e_TR_sample-enterotype.csv')
write.table(ent.Raes56.2e.sep.by.loc.RE$sample.enterotype, 'V5V6_2e_RE_sample-enterotype.csv')

# 3 enterotypes
ent.Raes12.3e.sep.by.loc.IL <- Enterotyping_Raes(input_file_v12, output_dir=output_dir, clust.number=3, noiserem=T, remove.outliers = 0.01, locat='.IL', sample.coverage.filter=coverage.thersholds['v12.IL'], plots=T)
ent.Raes12.3e.sep.by.loc.TR <- Enterotyping_Raes(input_file_v12, output_dir=output_dir, clust.number=3, noiserem=T, remove.outliers = 0.01, locat='.TR', sample.coverage.filter=coverage.thersholds['v12.TR'], plots=T)
ent.Raes12.3e.sep.by.loc.RE <- Enterotyping_Raes(input_file_v12, output_dir=output_dir, clust.number=3, noiserem=T, remove.outliers = 0.01, locat='.RE', sample.coverage.filter=coverage.thersholds['v12.RE'], plots=T)
ent.Raes34.3e.sep.by.loc.IL <- Enterotyping_Raes(input_file_v34, output_dir=output_dir, clust.number=3, noiserem=T, remove.outliers = 0.01, locat='.IL', sample.coverage.filter=coverage.thersholds['v34.IL'], plots=T)
ent.Raes34.3e.sep.by.loc.TR <- Enterotyping_Raes(input_file_v34, output_dir=output_dir, clust.number=3, noiserem=T, remove.outliers = 0.01, locat='.TR', sample.coverage.filter=coverage.thersholds['v34.TR'], plots=T)
ent.Raes34.3e.sep.by.loc.RE <- Enterotyping_Raes(input_file_v34, output_dir=output_dir, clust.number=3, noiserem=T, remove.outliers = 0.01, locat='.RE', sample.coverage.filter=coverage.thersholds['v34.RE'], plots=T)
ent.Raes56.3e.sep.by.loc.IL <- Enterotyping_Raes(input_file_v56, output_dir=output_dir, clust.number=3, noiserem=T, remove.outliers = 0.01, locat='.IL', sample.coverage.filter=coverage.thersholds['v56.IL'], plots=T)
ent.Raes56.3e.sep.by.loc.TR <- Enterotyping_Raes(input_file_v56, output_dir=output_dir, clust.number=3, noiserem=T, remove.outliers = 0.01, locat='.TR', sample.coverage.filter=coverage.thersholds['v56.TR'], plots=T)
ent.Raes56.3e.sep.by.loc.RE <- Enterotyping_Raes(input_file_v56, output_dir=output_dir, clust.number=3, noiserem=T, remove.outliers = 0.01, locat='.RE', sample.coverage.filter=coverage.thersholds['v56.RE'], plots=T)

setwd(file.path(home.dir,'Results/Enterotyping/2018_11_30_enterotyping/Sample-enterotype'))

write.table(ent.Raes12.3e.sep.by.loc.IL$sample.enterotype, 'V1V2_3e_IL_sample-enterotype.csv')
write.table(ent.Raes12.3e.sep.by.loc.TR$sample.enterotype, 'V1V2_3e_TR_sample-enterotype.csv')
write.table(ent.Raes12.3e.sep.by.loc.RE$sample.enterotype, 'V1V2_3e_RE_sample-enterotype.csv')

write.table(ent.Raes34.3e.sep.by.loc.IL$sample.enterotype, 'V3V4_3e_IL_sample-enterotype.csv')
write.table(ent.Raes34.3e.sep.by.loc.TR$sample.enterotype, 'V3V4_3e_TR_sample-enterotype.csv')
write.table(ent.Raes34.3e.sep.by.loc.RE$sample.enterotype, 'V3V4_3e_RE_sample-enterotype.csv')

write.table(ent.Raes56.3e.sep.by.loc.IL$sample.enterotype, 'V5V6_3e_IL_sample-enterotype.csv')
write.table(ent.Raes56.3e.sep.by.loc.TR$sample.enterotype, 'V5V6_3e_TR_sample-enterotype.csv')
write.table(ent.Raes56.3e.sep.by.loc.RE$sample.enterotype, 'V5V6_3e_RE_sample-enterotype.csv')


# 4 enterotypes
ent.Raes12.4e.sep.by.loc.IL <- Enterotyping_Raes(input_file_v12, output_dir=output_dir, clust.number=4, noiserem=T, remove.outliers = 0.01, locat='.IL', sample.coverage.filter=coverage.thersholds['v12.IL'], plots=T)
ent.Raes12.4e.sep.by.loc.TR <- Enterotyping_Raes(input_file_v12, output_dir=output_dir, clust.number=4, noiserem=T, remove.outliers = 0.01, locat='.TR', sample.coverage.filter=coverage.thersholds['v12.TR'], plots=T)
ent.Raes12.4e.sep.by.loc.RE <- Enterotyping_Raes(input_file_v12, output_dir=output_dir, clust.number=4, noiserem=T, remove.outliers = 0.01, locat='.RE', sample.coverage.filter=coverage.thersholds['v12.RE'], plots=T)
ent.Raes34.4e.sep.by.loc.IL <- Enterotyping_Raes(input_file_v34, output_dir=output_dir, clust.number=4, noiserem=T, remove.outliers = 0.01, locat='.IL', sample.coverage.filter=coverage.thersholds['v34.IL'], plots=T)
ent.Raes34.4e.sep.by.loc.TR <- Enterotyping_Raes(input_file_v34, output_dir=output_dir, clust.number=4, noiserem=T, remove.outliers = 0.01, locat='.TR', sample.coverage.filter=coverage.thersholds['v34.TR'], plots=T)
ent.Raes34.4e.sep.by.loc.RE <- Enterotyping_Raes(input_file_v34, output_dir=output_dir, clust.number=4, noiserem=T, remove.outliers = 0.01, locat='.RE', sample.coverage.filter=coverage.thersholds['v34.RE'], plots=T)
ent.Raes56.4e.sep.by.loc.IL <- Enterotyping_Raes(input_file_v56, output_dir=output_dir, clust.number=4, noiserem=T, remove.outliers = 0.01, locat='.IL', sample.coverage.filter=coverage.thersholds['v56.IL'], plots=T)
ent.Raes56.4e.sep.by.loc.TR <- Enterotyping_Raes(input_file_v56, output_dir=output_dir, clust.number=4, noiserem=T, remove.outliers = 0.01, locat='.TR', sample.coverage.filter=coverage.thersholds['v56.TR'], plots=T)
ent.Raes56.4e.sep.by.loc.RE <- Enterotyping_Raes(input_file_v56, output_dir=output_dir, clust.number=4, noiserem=T, remove.outliers = 0.01, locat='.RE', sample.coverage.filter=coverage.thersholds['v56.RE'], plots=T)

setwd(file.path(home.dir,'Results/Enterotyping/2018_11_30_enterotyping/Sample-enterotype'))

write.table(ent.Raes12.4e.sep.by.loc.IL$sample.enterotype, 'V1V2_4e_IL_sample-enterotype.csv')
write.table(ent.Raes12.4e.sep.by.loc.TR$sample.enterotype, 'V1V2_4e_TR_sample-enterotype.csv')
write.table(ent.Raes12.4e.sep.by.loc.RE$sample.enterotype, 'V1V2_4e_RE_sample-enterotype.csv')

write.table(ent.Raes34.4e.sep.by.loc.IL$sample.enterotype, 'V3V4_4e_IL_sample-enterotype.csv')
write.table(ent.Raes34.4e.sep.by.loc.TR$sample.enterotype, 'V3V4_4e_TR_sample-enterotype.csv')
write.table(ent.Raes34.4e.sep.by.loc.RE$sample.enterotype, 'V3V4_4e_RE_sample-enterotype.csv')

write.table(ent.Raes56.4e.sep.by.loc.IL$sample.enterotype, 'V5V6_4e_IL_sample-enterotype.csv')
write.table(ent.Raes56.4e.sep.by.loc.TR$sample.enterotype, 'V5V6_4e_TR_sample-enterotype.csv')
write.table(ent.Raes56.4e.sep.by.loc.RE$sample.enterotype, 'V5V6_4e_RE_sample-enterotype.csv')


# 5 enterotypes
ent.Raes12.5e.sep.by.loc.IL <- Enterotyping_Raes(input_file_v12, output_dir=output_dir, clust.number=5, noiserem=T, remove.outliers = 0.01, locat='.IL', sample.coverage.filter=coverage.thersholds['v12.IL'], plots=T)
ent.Raes12.5e.sep.by.loc.TR <- Enterotyping_Raes(input_file_v12, output_dir=output_dir, clust.number=5, noiserem=T, remove.outliers = 0.01, locat='.TR', sample.coverage.filter=coverage.thersholds['v12.TR'], plots=T)
ent.Raes12.5e.sep.by.loc.RE <- Enterotyping_Raes(input_file_v12, output_dir=output_dir, clust.number=5, noiserem=T, remove.outliers = 0.01, locat='.RE', sample.coverage.filter=coverage.thersholds['v12.RE'], plots=T)
ent.Raes34.5e.sep.by.loc.IL <- Enterotyping_Raes(input_file_v34, output_dir=output_dir, clust.number=5, noiserem=T, remove.outliers = 0.01, locat='.IL', sample.coverage.filter=coverage.thersholds['v34.IL'], plots=T)
ent.Raes34.5e.sep.by.loc.TR <- Enterotyping_Raes(input_file_v34, output_dir=output_dir, clust.number=5, noiserem=T, remove.outliers = 0.01, locat='.TR', sample.coverage.filter=coverage.thersholds['v34.TR'], plots=T)
ent.Raes34.5e.sep.by.loc.RE <- Enterotyping_Raes(input_file_v34, output_dir=output_dir, clust.number=5, noiserem=T, remove.outliers = 0.01, locat='.RE', sample.coverage.filter=coverage.thersholds['v34.RE'], plots=T)
ent.Raes56.5e.sep.by.loc.IL <- Enterotyping_Raes(input_file_v56, output_dir=output_dir, clust.number=5, noiserem=T, remove.outliers = 0.01, locat='.IL', sample.coverage.filter=coverage.thersholds['v56.IL'], plots=T)
ent.Raes56.5e.sep.by.loc.TR <- Enterotyping_Raes(input_file_v56, output_dir=output_dir, clust.number=5, noiserem=T, remove.outliers = 0.01, locat='.TR', sample.coverage.filter=coverage.thersholds['v56.TR'], plots=T)
ent.Raes56.5e.sep.by.loc.RE <- Enterotyping_Raes(input_file_v56, output_dir=output_dir, clust.number=5, noiserem=T, remove.outliers = 0.01, locat='.RE', sample.coverage.filter=coverage.thersholds['v56.RE'], plots=T)

setwd(file.path(home.dir,'Results/Enterotyping/2018_11_30_enterotyping/Sample-enterotype'))

write.table(ent.Raes12.5e.sep.by.loc.IL$sample.enterotype, 'V1V2_5e_IL_sample-enterotype.csv')
write.table(ent.Raes12.5e.sep.by.loc.TR$sample.enterotype, 'V1V2_5e_TR_sample-enterotype.csv')
write.table(ent.Raes12.5e.sep.by.loc.RE$sample.enterotype, 'V1V2_5e_RE_sample-enterotype.csv')

write.table(ent.Raes34.5e.sep.by.loc.IL$sample.enterotype, 'V3V4_5e_IL_sample-enterotype.csv')
write.table(ent.Raes34.5e.sep.by.loc.TR$sample.enterotype, 'V3V4_5e_TR_sample-enterotype.csv')
write.table(ent.Raes34.5e.sep.by.loc.RE$sample.enterotype, 'V3V4_5e_RE_sample-enterotype.csv')

write.table(ent.Raes56.5e.sep.by.loc.IL$sample.enterotype, 'V5V6_5e_IL_sample-enterotype.csv')
write.table(ent.Raes56.5e.sep.by.loc.TR$sample.enterotype, 'V5V6_5e_TR_sample-enterotype.csv')
write.table(ent.Raes56.5e.sep.by.loc.RE$sample.enterotype, 'V5V6_5e_RE_sample-enterotype.csv')

