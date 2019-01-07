# v2 the only difference with v1 is that we draw plots
# Perform enterotyping: each location separately, 3 enterotypes
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Enterotyping_Raes_noise_removed_v9a.R')
setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype')


# 2 enterotypes
ent.Raes12.2e.sep.by.loc.IL <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs', clust.number=2, noiserem=T, locat='.IL', plots=T)
ent.Raes12.2e.sep.by.loc.TR <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs', clust.number=2, noiserem=T, locat='.TR', plots=T)
ent.Raes12.2e.sep.by.loc.RE <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs', clust.number=2, noiserem=T, locat='.RE', plots=T)
ent.Raes34.2e.sep.by.loc.IL <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs', clust.number=2, noiserem=T, locat='.IL', plots=T)
ent.Raes34.2e.sep.by.loc.TR <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs', clust.number=2, noiserem=T, locat='.TR', plots=T)
ent.Raes34.2e.sep.by.loc.RE <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs', clust.number=2, noiserem=T, locat='.RE', plots=T)
ent.Raes56.2e.sep.by.loc.IL <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs', clust.number=2, noiserem=T, locat='.IL', plots=T)
ent.Raes56.2e.sep.by.loc.TR <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs', clust.number=2, noiserem=T, locat='.TR', plots=T)
ent.Raes56.2e.sep.by.loc.RE <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs', clust.number=2, noiserem=T, locat='.RE', plots=T)

setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype')
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
ent.Raes12.3e.sep.by.loc.IL <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs', clust.number=3, noiserem=T, locat='.IL', plots=T)
ent.Raes12.3e.sep.by.loc.TR <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs', clust.number=3, noiserem=T, locat='.TR', plots=T)
ent.Raes12.3e.sep.by.loc.RE <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs', clust.number=3, noiserem=T, locat='.RE', plots=T)
ent.Raes34.3e.sep.by.loc.IL <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs', clust.number=3, noiserem=T, locat='.IL', plots=T)
ent.Raes34.3e.sep.by.loc.TR <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs', clust.number=3, noiserem=T, locat='.TR', plots=T)
ent.Raes34.3e.sep.by.loc.RE <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs', clust.number=3, noiserem=T, locat='.RE', plots=T)
ent.Raes56.3e.sep.by.loc.IL <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs', clust.number=3, noiserem=T, locat='.IL', plots=T)
ent.Raes56.3e.sep.by.loc.TR <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs', clust.number=3, noiserem=T, locat='.TR', plots=T)
ent.Raes56.3e.sep.by.loc.RE <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs', clust.number=3, noiserem=T, locat='.RE', plots=T)

setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype')
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
ent.Raes12.4e.sep.by.loc.IL <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs', clust.number=4, noiserem=T, locat='.IL', plots=T)
ent.Raes12.4e.sep.by.loc.TR <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs', clust.number=4, noiserem=T, locat='.TR', plots=T)
ent.Raes12.4e.sep.by.loc.RE <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs', clust.number=4, noiserem=T, locat='.RE', plots=T)
ent.Raes34.4e.sep.by.loc.IL <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs', clust.number=4, noiserem=T, locat='.IL', plots=T)
ent.Raes34.4e.sep.by.loc.TR <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs', clust.number=4, noiserem=T, locat='.TR', plots=T)
ent.Raes34.4e.sep.by.loc.RE <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs', clust.number=4, noiserem=T, locat='.RE', plots=T)
ent.Raes56.4e.sep.by.loc.IL <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs', clust.number=4, noiserem=T, locat='.IL', plots=T)
ent.Raes56.4e.sep.by.loc.TR <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs', clust.number=4, noiserem=T, locat='.TR', plots=T)
ent.Raes56.4e.sep.by.loc.RE <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs', clust.number=4, noiserem=T, locat='.RE', plots=T)

setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype')
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
ent.Raes12.5e.sep.by.loc.IL <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs', clust.number=5, noiserem=T, locat='.IL', plots=T)
ent.Raes12.5e.sep.by.loc.TR <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs', clust.number=5, noiserem=T, locat='.TR', plots=T)
ent.Raes12.5e.sep.by.loc.RE <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs', clust.number=5, noiserem=T, locat='.RE', plots=T)
ent.Raes34.5e.sep.by.loc.IL <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs', clust.number=5, noiserem=T, locat='.IL', plots=T)
ent.Raes34.5e.sep.by.loc.TR <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs', clust.number=5, noiserem=T, locat='.TR', plots=T)
ent.Raes34.5e.sep.by.loc.RE <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs', clust.number=5, noiserem=T, locat='.RE', plots=T)
ent.Raes56.5e.sep.by.loc.IL <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs', clust.number=5, noiserem=T, locat='.IL', plots=T)
ent.Raes56.5e.sep.by.loc.TR <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs', clust.number=5, noiserem=T, locat='.TR', plots=T)
ent.Raes56.5e.sep.by.loc.RE <- Enterotyping_Raes('/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt', '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs', clust.number=5, noiserem=T, locat='.RE', plots=T)

setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype')
write.table(ent.Raes12.5e.sep.by.loc.IL$sample.enterotype, 'V1V2_5e_IL_sample-enterotype.csv')
write.table(ent.Raes12.5e.sep.by.loc.TR$sample.enterotype, 'V1V2_5e_TR_sample-enterotype.csv')
write.table(ent.Raes12.5e.sep.by.loc.RE$sample.enterotype, 'V1V2_5e_RE_sample-enterotype.csv')

write.table(ent.Raes34.5e.sep.by.loc.IL$sample.enterotype, 'V3V4_5e_IL_sample-enterotype.csv')
write.table(ent.Raes34.5e.sep.by.loc.TR$sample.enterotype, 'V3V4_5e_TR_sample-enterotype.csv')
write.table(ent.Raes34.5e.sep.by.loc.RE$sample.enterotype, 'V3V4_5e_RE_sample-enterotype.csv')

write.table(ent.Raes56.5e.sep.by.loc.IL$sample.enterotype, 'V5V6_5e_IL_sample-enterotype.csv')
write.table(ent.Raes56.5e.sep.by.loc.TR$sample.enterotype, 'V5V6_5e_TR_sample-enterotype.csv')
write.table(ent.Raes56.5e.sep.by.loc.RE$sample.enterotype, 'V5V6_5e_RE_sample-enterotype.csv')

