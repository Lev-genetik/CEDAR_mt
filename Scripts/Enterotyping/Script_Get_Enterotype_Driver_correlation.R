# Get correlations of drivers of all enterotypes and write them to a file
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Enterotyping_compar_functions_v5.R')
file1='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt'
file2='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt'
file3='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt'
# we need ent.driver.corr.matrix[[1]] in case check=TRUE
setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping')
ent.driver.corr.matrix <- getEntCorr(file1,file2,file3,correl.method = 'pearson',rnd=3,check=T)
write.table(as.data.frame(ent.driver.corr.matrix$`correlation matrix for all enterotypes`),'Enterotype_driver_abundance_correlation_Pearson.csv',sep='\t')
ent.driver.corr.matrix <- getEntCorr(file1,file2,file3,correl.method = 'spearman',rnd=3,check=T)
write.table(as.data.frame(ent.driver.corr.matrix$`correlation matrix for all enterotypes`),'Enterotype_driver_abundance_correlation_Spearman.csv',sep='\t')
# and now Spearman correlation where Bacteroides and Prevotella are removed
ent.driver.corr.matrix.trunc <- getEntCorr(file1,file2,file3,correl.method = 'spearman',rnd=3,check=T,remove.driver=c('g__Bacteroides','g__Prevotella'))
write.table(as.data.frame(ent.driver.corr.matrix.trunc$`correlation matrix for all enterotypes`),'Enterotype_driver_abundance_correlation_Spearman_g_Bact_g_Prev_dropped.csv',sep='\t')


# Now, do the same for separate enterotyping in IL, TR and RE
# Checked that actually the driver vector elements do not have normal distribution
# (e.g. e12_1i,e12_2i) => can not use Pearson test!!!
setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Driver_correlat')
ent.driver.corr.matrix.IL <- getEntCorr(file1,file2,file3,locat='.IL',correl.method = 'spearman',rnd=3,check=T)
write.table(as.data.frame(ent.driver.corr.matrix.IL$`correlation matrix for all enterotypes`),'Enterotype_driver_abundance_correlation_Spearman_IL.csv')
ent.driver.corr.matrix.TR <- getEntCorr(file1,file2,file3,locat='.TR',correl.method = 'spearman',rnd=3,check=T)
write.table(as.data.frame(ent.driver.corr.matrix.TR$`correlation matrix for all enterotypes`),'Enterotype_driver_abundance_correlation_Spearman_TR.csv')
ent.driver.corr.matrix.RE <- getEntCorr(file1,file2,file3,locat='.RE',correl.method = 'spearman',rnd=3,check=T)
write.table(as.data.frame(ent.driver.corr.matrix.RE$`correlation matrix for all enterotypes`),'Enterotype_driver_abundance_correlation_Spearman_RE.csv')