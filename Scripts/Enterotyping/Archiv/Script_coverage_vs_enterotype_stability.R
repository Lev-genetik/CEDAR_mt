# We have found that sample coverage influences enterotyping stability.
# ent stability calculated both for locations together ('all') and separately (IL,TR,RE)
# 1. Here, we determine correlation of cov. with enterotype stability
# 2. and make plots of coverage of stable vs unstable samples(2)
# 3. Plot enterotype stability vs number of reads per sample + sample coverage for all samples with this amplicon
# Than, the thresholds are determined (iteration 1) and
# 4. Plot SMART enterotype stability (low-covered samples from other amplicons excluded) vs number of reads per sample
# 5. Stable and unstable samples ae calculated for v1v2 versus only well-covered v3v4 and v5v6.
# Second iteration thresholds determined visually from these graphs

#1. Correlate number of reads per sample with enterotype stability
home.dir <- '/media/lev-genetik/980E73270E72FD96/Liege'
progr.dir <- '/home/lev-genetik/Desktop/Projects/liege/src/Lev'
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Check_data_lit_review_v3.R')
library(data.table)
setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Check_liter_review/Enterotype_coverage')
locations234 <- c('.IL','.TR','.RE')
stable_enterotype_file <- c('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/Old-the-same/Stable_enterotypes_IL.csv',
                            '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/Old-the-same/Stable_enterotypes_TR.csv',
                            '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/Old-the-same/Stable_enterotypes_RE.csv')
# V1V2
v1v2_cov_vs_ent_stab_IL <- rep('',times=4)
v1v2_cov_vs_ent_stab_TR <- rep('',times=4)
v1v2_cov_vs_ent_stab_RE <- rep('',times=4)
v1v2_cov_vs_ent_stab <- rbind(v1v2_cov_vs_ent_stab_IL,v1v2_cov_vs_ent_stab_TR,v1v2_cov_vs_ent_stab_RE)
for(i in 1:dim(v1v2_cov_vs_ent_stab)[1]){
  v1v2_cov_vs_ent_stab[i,] <- corrCoverEnterotype(qiime_file = file.path(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt'),
                                                  locat = locations234[i], stable_enterotype_file = stable_enterotype_file[i])
}
colnames(v1v2_cov_vs_ent_stab) <- c('Mann-Whitney_res','t.test_res', 'mean_cov_stable', 'mean_cov_unstable')
write.table(v1v2_cov_vs_ent_stab,file = 'v1v2_coverage_vs_ent_stabil.csv',append = F,quote = F,sep = '\t',eol='\n') 
# V3V4
v3v4_cov_vs_ent_stab_IL <- rep('',times=4)
v3v4_cov_vs_ent_stab_TR <- rep('',times=4)
v3v4_cov_vs_ent_stab_RE <- rep('',times=4)
v3v4_cov_vs_ent_stab <- rbind(v3v4_cov_vs_ent_stab_IL,v3v4_cov_vs_ent_stab_TR,v3v4_cov_vs_ent_stab_RE)
for(i in 1:dim(v3v4_cov_vs_ent_stab)[1]){
  v3v4_cov_vs_ent_stab[i,] <- corrCoverEnterotype(qiime_file = file.path(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt'),
                                                  locat = locations234[i], stable_enterotype_file = stable_enterotype_file[i])
}
colnames(v3v4_cov_vs_ent_stab) <- c('Mann-Whitney_res','t.test_res', 'mean_cov_stable', 'mean_cov_unstable')
write.table(v3v4_cov_vs_ent_stab,file = 'v3v4_coverage_vs_ent_stabil.csv',append = F,quote = F,sep = '\t',eol='\n') 
#v5v6
v5v6_cov_vs_ent_stab_IL <- rep('',times=4)
v5v6_cov_vs_ent_stab_TR <- rep('',times=4)
v5v6_cov_vs_ent_stab_RE <- rep('',times=4)
v5v6_cov_vs_ent_stab <- rbind(v5v6_cov_vs_ent_stab_IL,v5v6_cov_vs_ent_stab_TR,v5v6_cov_vs_ent_stab_RE)
for(i in 1:dim(v5v6_cov_vs_ent_stab)[1]){
  v5v6_cov_vs_ent_stab[i,] <- corrCoverEnterotype(qiime_file = file.path(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt'),
                                                  locat = locations234[i], stable_enterotype_file = stable_enterotype_file[i])
}
colnames(v5v6_cov_vs_ent_stab) <- c('Mann-Whitney_res','t.test_res', 'mean_cov_stable', 'mean_cov_unstable')
write.table(v5v6_cov_vs_ent_stab,file = 'v5v6_coverage_vs_ent_stabil.csv',append = F,quote = F,sep = '\t',eol='\n') 


#2. Plot enterotype stability vs number of reads per sample 
# the classic version with kernel density, done by 
home.dir <- '/media/lev-genetik/980E73270E72FD96/Liege'
progr.dir <- '/home/lev-genetik/Desktop/Projects/liege/src/Lev'
source(file.path(progr.dir,'Check_data_lit_review_v3.R'))
setwd(file.path(home.dir,'/Results/Check_liter_review/Enterotype_coverage/Graphs'))
locations <- c('all','.IL','.TR','.RE')
stable_enterotype_file <- c(file.path(home.dir,'/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype/Stable_enterotypes.csv'),
                            file.path(home.dir,'/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/Old-the-same/Stable_enterotypes_IL.csv'),
                            file.path(home.dir,'/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/Old-the-same/Stable_enterotypes_TR.csv'),
                            file.path(home.dir,'/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/Old-the-same/Stable_enterotypes_RE.csv'))
# V1V2
jpeg('v1v2_enterot_stabil_vs_sample_coverage.jpg',width = 800,height = 800)
par(mfrow=c(2,2), mar=c(7,7,7,4),mgp=c(4,2,0))
for(i in 1:length(locations)){
  plotCoverEnterotype(qiime_file = file.path(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt'), 
                      locat = locations[i], stable_enterotype_file = stable_enterotype_file[i])
}
dev.off()
# V3V4
jpeg('v3v4_enterot_stabil_vs_sample_coverage.jpg',width = 800,height = 800)
par(mfrow=c(2,2), mar=c(7,7,7,4),mgp=c(4,2,0))
for(i in 1:length(locations)){
  plotCoverEnterotype(qiime_file = file.path(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt'), 
                      locat = locations[i], stable_enterotype_file = stable_enterotype_file[i])
}
dev.off()
# V5V6
jpeg('v5v6_enterot_stabil_vs_sample_coverage.jpg',width = 800,height = 800)
par(mfrow=c(2,2), mar=c(7,7,7,4),mgp=c(4,2,0))
for(i in 1:length(locations)){
  plotCoverEnterotype(qiime_file = file.path(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt'), 
                      locat = locations[i], stable_enterotype_file = stable_enterotype_file[i])
}
dev.off()

# the new version with checking directly how many samples have coverage n+-10k reads
home.dir <- '/media/lev-genetik/980E73270E72FD96/Liege'
progr.dir <- '/home/lev-genetik/Desktop/Projects/liege/src/Lev'
source(file.path(progr.dir,'Check_data_lit_review_v3.R'))
setwd(file.path(home.dir,'/Results/Check_liter_review/Enterotype_coverage/Graphs'))
locations <- c('all','.IL','.TR','.RE')
stable_enterotype_file <- c(file.path(home.dir,'/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype/Stable_enterotypes.csv'),
                            file.path(home.dir,'/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/Old-the-same/Stable_enterotypes_IL.csv'),
                            file.path(home.dir,'/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/Old-the-same/Stable_enterotypes_TR.csv'),
                            file.path(home.dir,'/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/Old-the-same/Stable_enterotypes_RE.csv'))
# V1V2
jpeg('v1v2_enterot_stabil_vs_sample_coverage_window_10k.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2), mar=c(8,10,7,5), mgp=c(6,2,0))
qiime_file = file.path(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt')
for(i in 1:length(stable_enterotype_file)){
  # From these, we need colnames(data) = samples. A bit long time, but do not care
  data <- readQiimeSmart(input_file=qiime_file, locat = locations[i], percent_threshold = 0,remove.outliers = 0,to.scale=F)
  colnames(data) <- gsub('.V[0-9][0-9]','',colnames(data))
  # get stably enterotyped samples   
  stable.enterotypes.current <- stable_enterotype_file[i] 
  stable.enterotypes.table <- read.table(stable.enterotypes.current)
  stable.samples <- as.vector(stable.enterotypes.table[,1])
  # get unstable samples
  unstable.samples <- colnames(data)[!colnames(data) %in% stable.samples]
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=20,to.append=FALSE,color='black',y_axis_length_factor = 1.5, x_axis_length_factor = 1.03,
                  my.subset=stable.samples,adjust.factor=1000/length(stable.samples))
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=20,to.append=TRUE,color='red',my.subset=unstable.samples,adjust.factor=1000/length(unstable.samples))
}
dev.off()
# V3V4
jpeg('v3v4_enterot_stabil_vs_sample_coverage_window_10k.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2), mar=c(8,10,7,5), mgp=c(6,2,0))
qiime_file = file.path(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt')
for(i in 1:length(stable_enterotype_file)){
  # From these, we need colnames(data) = samples. A bit long time, but do not care
  data <- readQiimeSmart(input_file=qiime_file, locat = locations[i], percent_threshold = 0,remove.outliers = 0,to.scale=F)
  colnames(data) <- gsub('.V[0-9][0-9]','',colnames(data))
  # get stably enterotyped samples   
  stable.enterotypes.current <- stable_enterotype_file[i] 
  stable.enterotypes.table <- read.table(stable.enterotypes.current)
  stable.samples <- as.vector(stable.enterotypes.table[,1])
  # get unstable samples
  unstable.samples <- colnames(data)[!colnames(data) %in% stable.samples]
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=20,to.append=FALSE,color='black',y_axis_length_factor = 1.5, x_axis_length_factor = 1.03,
                  my.subset=stable.samples,adjust.factor=1000/length(stable.samples))
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=20,to.append=TRUE,color='red',my.subset=unstable.samples,adjust.factor=1000/length(unstable.samples))
}
dev.off()
# V5V6
jpeg('v5v6_enterot_stabil_vs_sample_coverage_window_10k.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2), mar=c(8,10,7,5), mgp=c(6,2,0))
qiime_file = file.path(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt')
for(i in 1:length(stable_enterotype_file)){
  # From these, we need colnames(data) = samples. A bit long time, but do not care
  data <- readQiimeSmart(input_file=qiime_file, locat = locations[i], percent_threshold = 0,remove.outliers = 0,to.scale=F)
  colnames(data) <- gsub('.V[0-9][0-9]','',colnames(data))
  # get stably enterotyped samples   
  stable.enterotypes.current <- stable_enterotype_file[i] 
  stable.enterotypes.table <- read.table(stable.enterotypes.current)
  stable.samples <- as.vector(stable.enterotypes.table[,1])
  # get unstable samples
  unstable.samples <- colnames(data)[!colnames(data) %in% stable.samples]
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=20,to.append=FALSE,color='black',y_axis_length_factor = 1.5, x_axis_length_factor = 1.03,
                  my.subset=stable.samples,adjust.factor=1000/length(stable.samples))
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=20,to.append=TRUE,color='red',my.subset=unstable.samples,adjust.factor=1000/length(unstable.samples))
}
dev.off()

# the new version with checking directly how many samples have coverage n+-5k reads
# V1V2
jpeg('v1v2_enterot_stabil_vs_sample_coverage_window_5k.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2), mar=c(8,10,7,5), mgp=c(6,2,0))
qiime_file = file.path(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt')
for(i in 1:length(stable_enterotype_file)){
  # From these, we need colnames(data) = samples. A bit long time, but do not care
  data <- readQiimeSmart(input_file=qiime_file, locat = locations[i], percent_threshold = 0,remove.outliers = 0,to.scale=F)
  colnames(data) <- gsub('.V[0-9][0-9]','',colnames(data))
  # get stably enterotyped samples   
  stable.enterotypes.current <- stable_enterotype_file[i] 
  stable.enterotypes.table <- read.table(stable.enterotypes.current)
  stable.samples <- as.vector(stable.enterotypes.table[,1])
  # get unstable samples
  unstable.samples <- colnames(data)[!colnames(data) %in% stable.samples]
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=10,to.append=FALSE,color='black',y_axis_length_factor = 1.5, x_axis_length_factor = 1.03,
                  my.subset=stable.samples,adjust.factor=1000/length(stable.samples))
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=10,to.append=TRUE,color='red',my.subset=unstable.samples,adjust.factor=1000/length(unstable.samples))
}
dev.off()
# V3V4
jpeg('v3v4_enterot_stabil_vs_sample_coverage_window_5k.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2), mar=c(8,10,7,5), mgp=c(6,2,0))
qiime_file = file.path(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt')
for(i in 1:length(stable_enterotype_file)){
  # From these, we need colnames(data) = samples. A bit long time, but do not care
  data <- readQiimeSmart(input_file=qiime_file, locat = locations[i], percent_threshold = 0,remove.outliers = 0,to.scale=F)
  colnames(data) <- gsub('.V[0-9][0-9]','',colnames(data))
  # get stably enterotyped samples   
  stable.enterotypes.current <- stable_enterotype_file[i] 
  stable.enterotypes.table <- read.table(stable.enterotypes.current)
  stable.samples <- as.vector(stable.enterotypes.table[,1])
  # get unstable samples
  unstable.samples <- colnames(data)[!colnames(data) %in% stable.samples]
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=10,to.append=FALSE,color='black',y_axis_length_factor = 1.5, x_axis_length_factor = 1.03,
                  my.subset=stable.samples,adjust.factor=1000/length(stable.samples))
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=10,to.append=TRUE,color='red',my.subset=unstable.samples,adjust.factor=1000/length(unstable.samples))
}
dev.off()
# V5V6
jpeg('v5v6_enterot_stabil_vs_sample_coverage_window_5k.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2), mar=c(8,10,7,5), mgp=c(6,2,0))
qiime_file = file.path(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt')
for(i in 1:length(stable_enterotype_file)){
  # From these, we need colnames(data) = samples. A bit long time, but do not care
  data <- readQiimeSmart(input_file=qiime_file, locat = locations[i], percent_threshold = 0,remove.outliers = 0,to.scale=F)
  colnames(data) <- gsub('.V[0-9][0-9]','',colnames(data))
  # get stably enterotyped samples   
  stable.enterotypes.current <- stable_enterotype_file[i] 
  stable.enterotypes.table <- read.table(stable.enterotypes.current)
  stable.samples <- as.vector(stable.enterotypes.table[,1])
  # get unstable samples
  unstable.samples <- colnames(data)[!colnames(data) %in% stable.samples]
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=10,to.append=FALSE,color='black',y_axis_length_factor = 1.5, x_axis_length_factor = 1.03,
                  my.subset=stable.samples,adjust.factor=1000/length(stable.samples))
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=10,to.append=TRUE,color='red',my.subset=unstable.samples,adjust.factor=1000/length(unstable.samples))
}
dev.off()

#3. Plot enterotype stability vs number of reads per sample + sample coverage for all samples with this amplicon
home.dir <- '/media/lev-genetik/980E73270E72FD96/Liege'
progr.dir <- '/home/lev-genetik/Desktop/Projects/liege/src/Lev'
source(file.path(progr.dir,'Check_data_lit_review_v3.R'))
setwd(file.path(home.dir,'/Results/Check_liter_review/Enterotype_coverage/Graphs'))
locations <- c('all','.IL','.TR','.RE')
stable_enterotype_file <- c(file.path(home.dir,'/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype/Stable_enterotypes.csv'),
                            file.path(home.dir,'/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/Old-the-same/Stable_enterotypes_IL.csv'),
                            file.path(home.dir,'/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/Old-the-same/Stable_enterotypes_TR.csv'),
                            file.path(home.dir,'/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/Old-the-same/Stable_enterotypes_RE.csv'))
# V1V2
jpeg('v1v2_enterot_stabil_vs_sample_coverage_window_10k_and_microbiota_cov.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2), mar=c(8,10,7,5), mgp=c(6,2,0))
qiime_file = file.path(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt')
for(i in 1:length(locations)){
  # From these, we need colnames(data) = samples. A bit long time, but do not care
  data <- readQiimeSmart(input_file=qiime_file, locat = locations[i], percent_threshold = 0,remove.outliers = 0,to.scale=F)
  colnames(data) <- gsub('.V[0-9][0-9]','',colnames(data))
  # get stably enterotyped samples   
  stable.enterotypes.current <- stable_enterotype_file[i] 
  stable.enterotypes.table <- read.table(stable.enterotypes.current)
  stable.samples <- as.vector(stable.enterotypes.table[,1])
  # get unstable samples
  unstable.samples <- colnames(data)[!colnames(data) %in% stable.samples]
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=20,to.append=FALSE,color='black',y_axis_length_factor = 1.5, x_axis_length_factor = 1.03,
                  my.subset=stable.samples,adjust.factor=(length(stable.samples)+length(unstable.samples))/length(stable.samples))
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=20,to.append=TRUE,color='red',my.subset=unstable.samples,adjust.factor=(length(stable.samples)+length(unstable.samples))/length(unstable.samples))
  exploreCoverage(input_type='QIIME', input_file = qiime_file,
                  locat = locations[i], percent_threshold = 0, remove.outliers = 0,cumulative=FALSE,usredneniye=20,y_axis_length_factor = 1,
                  to.append = T, color = 'green')
  # exploreCoverage(input_type='total.coverage', input_file = qiime_file,
  #                 locat = locations[i], percent_threshold = 0, remove.outliers = 0,cumulative=FALSE,usredneniye=20,
  #                 to.append = T,color = 'purple')
  
}
dev.off()
# v3v4
jpeg('v3v4_enterot_stabil_vs_sample_coverage_window_10k_and_microbiota_cov.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2), mar=c(8,10,7,5), mgp=c(6,2,0))
qiime_file = file.path(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt')
for(i in 1:length(locations)){
  # From these, we need colnames(data) = samples. A bit long time, but do not care
  data <- readQiimeSmart(input_file=qiime_file, locat = locations[i], percent_threshold = 0,remove.outliers = 0,to.scale=F)
  colnames(data) <- gsub('.V[0-9][0-9]','',colnames(data))
  # get stably enterotyped samples   
  stable.enterotypes.current <- stable_enterotype_file[i] 
  stable.enterotypes.table <- read.table(stable.enterotypes.current)
  stable.samples <- as.vector(stable.enterotypes.table[,1])
  # get unstable samples
  unstable.samples <- colnames(data)[!colnames(data) %in% stable.samples]
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=20,to.append=FALSE,color='black',y_axis_length_factor = 1.5, x_axis_length_factor = 1.03,
                  my.subset=stable.samples,adjust.factor=(length(stable.samples)+length(unstable.samples))/length(stable.samples))
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=20,to.append=TRUE,color='red',my.subset=unstable.samples,adjust.factor=(length(stable.samples)+length(unstable.samples))/length(unstable.samples))
  exploreCoverage(input_type='QIIME', input_file = qiime_file,
                  locat = locations[i], percent_threshold = 0, remove.outliers = 0,cumulative=FALSE,usredneniye=20,y_axis_length_factor = 1,
                  to.append = T, color = 'green')
  # exploreCoverage(input_type='total.coverage', input_file = qiime_file,
  #                 locat = locations[i], percent_threshold = 0, remove.outliers = 0,cumulative=FALSE,usredneniye=20,
  #                 to.append = T,color = 'purple')
  
}
dev.off()
# V5V6
jpeg('v5v6_enterot_stabil_vs_sample_coverage_window_10k_and_microbiota_cov.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2), mar=c(8,10,7,5), mgp=c(6,2,0))
qiime_file = file.path(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt')
for(i in 1:length(locations)){
  # From these, we need colnames(data) = samples. A bit long time, but do not care
  data <- readQiimeSmart(input_file=qiime_file, locat = locations[i], percent_threshold = 0,remove.outliers = 0,to.scale=F)
  colnames(data) <- gsub('.V[0-9][0-9]','',colnames(data))
  # get stably enterotyped samples   
  stable.enterotypes.current <- stable_enterotype_file[i] 
  stable.enterotypes.table <- read.table(stable.enterotypes.current)
  stable.samples <- as.vector(stable.enterotypes.table[,1])
  # get unstable samples
  unstable.samples <- colnames(data)[!colnames(data) %in% stable.samples]
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=20,to.append=FALSE,color='black',y_axis_length_factor = 1.5, x_axis_length_factor = 1.03,
                  my.subset=stable.samples,adjust.factor=(length(stable.samples)+length(unstable.samples))/length(stable.samples))
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=20,to.append=TRUE,color='red',my.subset=unstable.samples,adjust.factor=(length(stable.samples)+length(unstable.samples))/length(unstable.samples))
  exploreCoverage(input_type='QIIME', input_file = qiime_file,
                  locat = locations[i], percent_threshold = 0, remove.outliers = 0,cumulative=FALSE,usredneniye=20,y_axis_length_factor = 1,
                  to.append = T, color = 'green')
  # exploreCoverage(input_type='total.coverage', input_file = qiime_file,
  #                 locat = locations[i], percent_threshold = 0, remove.outliers = 0,cumulative=FALSE,usredneniye=20,
  #                 to.append = T,color = 'purple')
  
}
dev.off()

#4. Plot SMART enterotype stability (low-covered samples from other amplicons excluded) vs number of reads per sample
home.dir <- '/media/lev-genetik/980E73270E72FD96/Liege'
progr.dir <- '/home/lev-genetik/Desktop/Projects/liege/src/Lev'
source(file.path(progr.dir,'Check_data_lit_review_v3.R'))
setwd(file.path(home.dir,'/Results/Check_liter_review/Enterotype_coverage/Graphs'))
qiime_file_v12 = file.path(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt')
qiime_file_v34 = file.path(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt')
qiime_file_v56 = file.path(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt')
locations <- c('all','.IL','.TR','.RE')
stable_enterotype_file <- c(file.path(home.dir,'/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype/Stable_enterotypes.csv'),
                            file.path(home.dir,'/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/Old-the-same/Stable_enterotypes_IL.csv'),
                            file.path(home.dir,'/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/Old-the-same/Stable_enterotypes_TR.csv'),
                            file.path(home.dir,'/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/Old-the-same/Stable_enterotypes_RE.csv'))
# get v1v2, v3v4 and v5v6 samples with poor and good coverage - iteration 1
# v1v2
data_v12_all <- readQiimeSmart(input_file = qiime_file_v12, locat = 'all', percent_threshold = 0, remove.outliers = 0,to.scale = FALSE)
colnames(data_v12_all) <- gsub('.V[0-9][0-9]','',colnames(data_v12_all))
data_v12_all.poor <- colnames(data_v12_all[,colSums(data_v12_all)<23000])
data_v12_all.good <- colnames(data_v12_all[,colSums(data_v12_all)>=23000])
data_v12_IL <- readQiimeSmart(input_file = qiime_file_v12, locat = '.IL', percent_threshold = 0, remove.outliers = 0,to.scale = FALSE)
colnames(data_v12_IL) <- gsub('.V[0-9][0-9]','',colnames(data_v12_IL))
data_v12_IL.poor <- colnames(data_v12_IL[,colSums(data_v12_IL)<23000])
data_v12_IL.good <- colnames(data_v12_IL[,colSums(data_v12_IL)>=23000])
data_v12_TR <- readQiimeSmart(input_file = qiime_file_v12, locat = '.TR', percent_threshold = 0, remove.outliers = 0,to.scale = FALSE)
colnames(data_v12_TR) <- gsub('.V[0-9][0-9]','',colnames(data_v12_TR))
data_v12_TR.poor <- colnames(data_v12_TR[,colSums(data_v12_TR)<23000])
data_v12_TR.good <- colnames(data_v12_TR[,colSums(data_v12_TR)>=23000])
data_v12_RE <- readQiimeSmart(input_file = qiime_file_v12, locat = '.RE', percent_threshold = 0, remove.outliers = 0,to.scale = FALSE)
colnames(data_v12_RE) <- gsub('.V[0-9][0-9]','',colnames(data_v12_RE))
data_v12_RE.poor <- colnames(data_v12_RE[,colSums(data_v12_RE)<22500])
data_v12_RE.good <- colnames(data_v12_RE[,colSums(data_v12_RE)>=22500])
data_v12.poor <- list(data_v12_all.poor,data_v12_IL.poor,data_v12_TR.poor,data_v12_RE.poor)
data_v12.good <- list(data_v12_all.good,data_v12_IL.good,data_v12_TR.good,data_v12_RE.good)
# v3v4
data_v34_all <- readQiimeSmart(input_file = qiime_file_v34, locat = 'all', percent_threshold = 0, remove.outliers = 0,to.scale = FALSE)
colnames(data_v34_all) <- gsub('.V[0-9][0-9]','',colnames(data_v34_all))
data_v34_all.poor <- colnames(data_v34_all[,colSums(data_v34_all)<17000])
data_v34_all.good <- colnames(data_v34_all[,colSums(data_v34_all)>=17000])
data_v34_IL <- readQiimeSmart(input_file = qiime_file_v34, locat = '.IL', percent_threshold = 0, remove.outliers = 0,to.scale = FALSE)
colnames(data_v34_IL) <- gsub('.V[0-9][0-9]','',colnames(data_v34_IL))
data_v34_IL.poor <- colnames(data_v34_IL[,colSums(data_v34_IL)<17000])
data_v34_IL.good <- colnames(data_v34_IL[,colSums(data_v34_IL)>=17000])
data_v34_TR <- readQiimeSmart(input_file = qiime_file_v34, locat = '.TR', percent_threshold = 0, remove.outliers = 0,to.scale = FALSE)
colnames(data_v34_TR) <- gsub('.V[0-9][0-9]','',colnames(data_v34_TR))
data_v34_TR.poor <- colnames(data_v34_TR[,colSums(data_v34_TR)<17000])
data_v34_TR.good <- colnames(data_v34_TR[,colSums(data_v34_TR)>=17000])
data_v34_RE <- readQiimeSmart(input_file = qiime_file_v34, locat = '.RE', percent_threshold = 0, remove.outliers = 0,to.scale = FALSE)
colnames(data_v34_RE) <- gsub('.V[0-9][0-9]','',colnames(data_v34_RE))
data_v34_RE.poor <- colnames(data_v34_RE[,colSums(data_v34_RE)<17500])
data_v34_RE.good <- colnames(data_v34_RE[,colSums(data_v34_RE)>=17500])
data_v34.poor <- list(data_v34_all.poor,data_v34_IL.poor,data_v34_TR.poor,data_v34_RE.poor)
data_v34.good <- list(data_v34_all.good,data_v34_IL.good,data_v34_TR.good,data_v34_RE.good)
# v5v6
data_v56_all <- readQiimeSmart(input_file = qiime_file_v56, locat = 'all', percent_threshold = 0, remove.outliers = 0,to.scale = FALSE)
colnames(data_v56_all) <- gsub('.V[0-9][0-9]','',colnames(data_v56_all))
data_v56_all.poor <- colnames(data_v56_all[,colSums(data_v56_all)<36000])
data_v56_all.good <- colnames(data_v56_all[,colSums(data_v56_all)>=36000])
data_v56_IL <- readQiimeSmart(input_file = qiime_file_v56, locat = '.IL', percent_threshold = 0, remove.outliers = 0,to.scale = FALSE)
colnames(data_v56_IL) <- gsub('.V[0-9][0-9]','',colnames(data_v56_IL))
data_v56_IL.poor <- colnames(data_v56_IL[,colSums(data_v56_IL)<36000])
data_v56_IL.good <- colnames(data_v56_IL[,colSums(data_v56_IL)>=36000])
data_v56_TR <- readQiimeSmart(input_file = qiime_file_v56, locat = '.TR', percent_threshold = 0, remove.outliers = 0,to.scale = FALSE)
colnames(data_v56_TR) <- gsub('.V[0-9][0-9]','',colnames(data_v56_TR))
data_v56_TR.poor <- colnames(data_v56_TR[,colSums(data_v56_TR)<34000])
data_v56_TR.good <- colnames(data_v56_TR[,colSums(data_v56_TR)>=34000])
data_v56_RE <- readQiimeSmart(input_file = qiime_file_v56, locat = '.RE', percent_threshold = 0, remove.outliers = 0,to.scale = FALSE)
colnames(data_v56_RE) <- gsub('.V[0-9][0-9]','',colnames(data_v56_RE))
data_v56_RE.poor <- colnames(data_v56_RE[,colSums(data_v56_RE)<31000])
data_v56_RE.good <- colnames(data_v56_RE[,colSums(data_v56_RE)>=31000])
data_v56.poor <- list(data_v56_all.poor,data_v56_IL.poor,data_v56_TR.poor,data_v56_RE.poor)
data_v56.good <- list(data_v56_all.good,data_v56_IL.good,data_v56_TR.good,data_v56_RE.good)
# V1V2
jpeg('v1v2_enterot_stabil_vs_sample_coverage_window_10k_only_v34_v56_good_cov.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2), mar=c(8,10,7,5), mgp=c(6,2,0))
qiime_file = file.path(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt')
for(i in 1:length(locations)){
  # From these, we need colnames(data) = samples. A bit long time, but do not care
  data <- readQiimeSmart(input_file=qiime_file, locat = locations[i], percent_threshold = 0,remove.outliers = 0,to.scale=F)
  colnames(data) <- gsub('.V[0-9][0-9]','',colnames(data))
  
  # get stably enterotyped samples  
  stable.enterotypes.current <- stable_enterotype_file[i] 
  stable.enterotypes.table <- read.table(stable.enterotypes.current)
  stable.samples.downsamped <- as.vector(stable.enterotypes.table[,1])[!(as.vector(stable.enterotypes.table[,1]) %in% union(data_v34.poor[[i]],data_v56.poor[[i]]))]
  # get unstably enterotyped samples
  unstable.samples <-  colnames(data_v12_all)[!colnames(data_v12_all) %in% as.vector(stable.enterotypes.table[,1])] #temprorarily
  # if locat != all we need to delete other locations
  if(locations[i]!='all'){unstable.samples <- unstable.samples[unstable.samples %like% locations[i]]}
  unstable.samples.downsamped <- unstable.samples[!(unstable.samples %in% union(data_v34.poor[[i]],data_v56.poor[[i]]))]
  # unstable.samples_downsamped <- as.vector(stable.enterotypes.table[,1])[!(as.vector(stable.enterotypes.table[,1]) %in% union(data_v34_all.good,data_v56_all.good))]
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=20,to.append=FALSE,color='black',y_axis_length_factor = 1.5, x_axis_length_factor = 1.03,
                  my.subset=stable.samples.downsamped,adjust.factor=(length(stable.samples.downsamped)+length(unstable.samples.downsamped))/length(stable.samples.downsamped))
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=20,to.append=TRUE,color='red',my.subset=unstable.samples.downsamped,adjust.factor=(length(stable.samples.downsamped)+length(unstable.samples.downsamped))/length(unstable.samples.downsamped))
}
dev.off()
# V3V4
jpeg('v3v4_enterot_stabil_vs_sample_coverage_window_10k_only_v12_v56_good_cov.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2), mar=c(8,10,7,5), mgp=c(6,2,0))
qiime_file = file.path(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt')
for(i in 1:length(locations)){
  # From these, we need colnames(data) = samples. A bit long time, but do not care
  data <- readQiimeSmart(input_file=qiime_file, locat = locations[i], percent_threshold = 0,remove.outliers = 0,to.scale=F)
  colnames(data) <- gsub('.V[0-9][0-9]','',colnames(data))
  
  # get stably enterotyped samples  
  stable.enterotypes.current <- stable_enterotype_file[i] 
  stable.enterotypes.table <- read.table(stable.enterotypes.current)
  stable.samples.downsamped <- as.vector(stable.enterotypes.table[,1])[!(as.vector(stable.enterotypes.table[,1]) %in% union(data_v12.poor[[i]],data_v56.poor[[i]]))]
  # get unstably enterotyped samples
  unstable.samples <-  colnames(data_v34_all)[!colnames(data_v34_all) %in% as.vector(stable.enterotypes.table[,1])] #temprorarily
  # if locat != all we need to delete other locations
  if(locations[i]!='all'){unstable.samples <- unstable.samples[unstable.samples %like% locations[i]]}
  unstable.samples.downsamped <- unstable.samples[!(unstable.samples %in% union(data_v12.poor[[i]],data_v56.poor[[i]]))]
  # unstable.samples_downsamped <- as.vector(stable.enterotypes.table[,1])[!(as.vector(stable.enterotypes.table[,1]) %in% union(data_v34_all.good,data_v56_all.good))]
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=20,to.append=FALSE,color='black',y_axis_length_factor = 1.5, x_axis_length_factor = 1.03,
                  my.subset=stable.samples.downsamped,adjust.factor=(length(stable.samples.downsamped)+length(unstable.samples.downsamped))/length(stable.samples.downsamped))
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=20,to.append=TRUE,color='red',my.subset=unstable.samples.downsamped,adjust.factor=(length(stable.samples.downsamped)+length(unstable.samples.downsamped))/length(unstable.samples.downsamped))
}
dev.off()
# V5V6
jpeg('v5v6_enterot_stabil_vs_sample_coverage_window_10k_only_v12_v34_good_cov.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2), mar=c(8,10,7,5), mgp=c(6,2,0))
qiime_file = file.path(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt')
for(i in 1:length(locations)){
  # From these, we need colnames(data) = samples. A bit long time, but do not care
  data <- readQiimeSmart(input_file=qiime_file, locat = locations[i], percent_threshold = 0,remove.outliers = 0,to.scale=F)
  colnames(data) <- gsub('.V[0-9][0-9]','',colnames(data))
  
  # get stably enterotyped samples  
  stable.enterotypes.current <- stable_enterotype_file[i] 
  stable.enterotypes.table <- read.table(stable.enterotypes.current)
  stable.samples.downsamped <- as.vector(stable.enterotypes.table[,1])[!(as.vector(stable.enterotypes.table[,1]) %in% union(data_v12.poor[[i]],data_v34.poor[[i]]))]
  # get unstably enterotyped samples
  unstable.samples <-  colnames(data_v12_all)[!colnames(data_v12_all) %in% as.vector(stable.enterotypes.table[,1])] #temprorarily
  # if locat != all we need to delete other locations
  if(locations[i]!='all'){unstable.samples <- unstable.samples[unstable.samples %like% locations[i]]}
  unstable.samples.downsamped <- unstable.samples[!(unstable.samples %in% union(data_v12.poor[[i]],data_v34.poor[[i]]))]
  # unstable.samples_downsamped <- as.vector(stable.enterotypes.table[,1])[!(as.vector(stable.enterotypes.table[,1]) %in% union(data_v34_all.good,data_v56_all.good))]
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=20,to.append=FALSE,color='black',y_axis_length_factor = 1.5, x_axis_length_factor = 1.03,
                  my.subset=stable.samples.downsamped,adjust.factor=(length(stable.samples.downsamped)+length(unstable.samples.downsamped))/length(stable.samples.downsamped))
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=20,to.append=TRUE,color='red',my.subset=unstable.samples.downsamped,adjust.factor=(length(stable.samples.downsamped)+length(unstable.samples.downsamped))/length(unstable.samples.downsamped))
}
dev.off()

# 6. Now, to determine the filtration threshold most precisely, we need to compare graphs for different rolling window size
coverages <- c(2,4,10,20)
# V1V2
i=1 #location
jpeg('v1v2_all_enterot_stabil_vs_sample_coverage_window_1k_2k_5k_10k_only_v34_v56_good_cov.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2), mar=c(8,10,7,5), mgp=c(6,2,0))
qiime_file = file.path(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt')
for(j in 1:length(coverages)){
  # From these, we need colnames(data) = samples. A bit long time, but do not care
  data <- readQiimeSmart(input_file=qiime_file, locat = locations[i], percent_threshold = 0,remove.outliers = 0,to.scale=F)
  colnames(data) <- gsub('.V[0-9][0-9]','',colnames(data))
  # get stably enterotyped samples  
  stable.enterotypes.current <- stable_enterotype_file[i] 
  stable.enterotypes.table <- read.table(stable.enterotypes.current)
  stable.samples.downsamped <- as.vector(stable.enterotypes.table[,1])[!(as.vector(stable.enterotypes.table[,1]) %in% union(data_v34.poor[[i]],data_v56.poor[[i]]))]
  # get unstably enterotyped samples
  unstable.samples <-  colnames(data_v12_all)[!colnames(data_v12_all) %in% as.vector(stable.enterotypes.table[,1])] #temprorarily
  # if locat != all we need to delete other locations
  if(locations[i]!='all'){unstable.samples <- unstable.samples[unstable.samples %like% locations[i]]}
  unstable.samples.downsamped <- unstable.samples[!(unstable.samples %in% union(data_v34.poor[[i]],data_v56.poor[[i]]))]
  # unstable.samples_downsamped <- as.vector(stable.enterotypes.table[,1])[!(as.vector(stable.enterotypes.table[,1]) %in% union(data_v34_all.good,data_v56_all.good))]
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=coverages[j],to.append=FALSE,color='black',y_axis_length_factor = 1.5, x_axis_length_factor = 1.03,
                  my.subset=stable.samples.downsamped,adjust.factor=(length(stable.samples.downsamped)+length(unstable.samples.downsamped))/length(stable.samples.downsamped))
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=coverages[j],to.append=TRUE,color='red',my.subset=unstable.samples.downsamped,adjust.factor=(length(stable.samples.downsamped)+length(unstable.samples.downsamped))/length(unstable.samples.downsamped))
}
dev.off()
i=2 #location IL
jpeg('v1v2_IL_enterot_stabil_vs_sample_coverage_window_1k_2k_5k_10k_only_v34_v56_good_cov.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2), mar=c(8,10,7,5), mgp=c(6,2,0))
qiime_file = file.path(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt')
for(j in 1:length(coverages)){
  # From these, we need colnames(data) = samples. A bit long time, but do not care
  data <- readQiimeSmart(input_file=qiime_file, locat = locations[i], percent_threshold = 0,remove.outliers = 0,to.scale=F)
  colnames(data) <- gsub('.V[0-9][0-9]','',colnames(data))
  # get stably enterotyped samples  
  stable.enterotypes.current <- stable_enterotype_file[i] 
  stable.enterotypes.table <- read.table(stable.enterotypes.current)
  stable.samples.downsamped <- as.vector(stable.enterotypes.table[,1])[!(as.vector(stable.enterotypes.table[,1]) %in% union(data_v34.poor[[i]],data_v56.poor[[i]]))]
  # get unstably enterotyped samples
  unstable.samples <-  colnames(data_v12_all)[!colnames(data_v12_all) %in% as.vector(stable.enterotypes.table[,1])] #temprorarily
  # if locat != all we need to delete other locations
  if(locations[i]!='all'){unstable.samples <- unstable.samples[unstable.samples %like% locations[i]]}
  unstable.samples.downsamped <- unstable.samples[!(unstable.samples %in% union(data_v34.poor[[i]],data_v56.poor[[i]]))]
  # unstable.samples_downsamped <- as.vector(stable.enterotypes.table[,1])[!(as.vector(stable.enterotypes.table[,1]) %in% union(data_v34_all.good,data_v56_all.good))]
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=coverages[j],to.append=FALSE,color='black',y_axis_length_factor = 1.5, x_axis_length_factor = 1.03,
                  my.subset=stable.samples.downsamped,adjust.factor=(length(stable.samples.downsamped)+length(unstable.samples.downsamped))/length(stable.samples.downsamped))
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=coverages[j],to.append=TRUE,color='red',my.subset=unstable.samples.downsamped,adjust.factor=(length(stable.samples.downsamped)+length(unstable.samples.downsamped))/length(unstable.samples.downsamped))
}
dev.off()
i=3 #location TR
jpeg('v1v2_TR_enterot_stabil_vs_sample_coverage_window_1k_2k_5k_10k_only_v34_v56_good_cov.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2), mar=c(8,10,7,5), mgp=c(6,2,0))
qiime_file = file.path(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt')
for(j in 1:length(coverages)){
  # From these, we need colnames(data) = samples. A bit long time, but do not care
  data <- readQiimeSmart(input_file=qiime_file, locat = locations[i], percent_threshold = 0,remove.outliers = 0,to.scale=F)
  colnames(data) <- gsub('.V[0-9][0-9]','',colnames(data))
  # get stably enterotyped samples  
  stable.enterotypes.current <- stable_enterotype_file[i] 
  stable.enterotypes.table <- read.table(stable.enterotypes.current)
  stable.samples.downsamped <- as.vector(stable.enterotypes.table[,1])[!(as.vector(stable.enterotypes.table[,1]) %in% union(data_v34.poor[[i]],data_v56.poor[[i]]))]
  # get unstably enterotyped samples
  unstable.samples <-  colnames(data_v12_all)[!colnames(data_v12_all) %in% as.vector(stable.enterotypes.table[,1])] #temprorarily
  # if locat != all we need to delete other locations
  if(locations[i]!='all'){unstable.samples <- unstable.samples[unstable.samples %like% locations[i]]}
  unstable.samples.downsamped <- unstable.samples[!(unstable.samples %in% union(data_v34.poor[[i]],data_v56.poor[[i]]))]
  # unstable.samples_downsamped <- as.vector(stable.enterotypes.table[,1])[!(as.vector(stable.enterotypes.table[,1]) %in% union(data_v34_all.good,data_v56_all.good))]
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=coverages[j],to.append=FALSE,color='black',y_axis_length_factor = 1.5, x_axis_length_factor = 1.03,
                  my.subset=stable.samples.downsamped,adjust.factor=(length(stable.samples.downsamped)+length(unstable.samples.downsamped))/length(stable.samples.downsamped))
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=coverages[j],to.append=TRUE,color='red',my.subset=unstable.samples.downsamped,adjust.factor=(length(stable.samples.downsamped)+length(unstable.samples.downsamped))/length(unstable.samples.downsamped))
}
dev.off()
i=4 #location RE
jpeg('v1v2_RE_enterot_stabil_vs_sample_coverage_window_1k_2k_5k_10k_only_v34_v56_good_cov.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2), mar=c(8,10,7,5), mgp=c(6,2,0))
qiime_file = file.path(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt')
for(j in 1:length(coverages)){
  # From these, we need colnames(data) = samples. A bit long time, but do not care
  data <- readQiimeSmart(input_file=qiime_file, locat = locations[i], percent_threshold = 0,remove.outliers = 0,to.scale=F)
  colnames(data) <- gsub('.V[0-9][0-9]','',colnames(data))
  # get stably enterotyped samples  
  stable.enterotypes.current <- stable_enterotype_file[i] 
  stable.enterotypes.table <- read.table(stable.enterotypes.current)
  stable.samples.downsamped <- as.vector(stable.enterotypes.table[,1])[!(as.vector(stable.enterotypes.table[,1]) %in% union(data_v34.poor[[i]],data_v56.poor[[i]]))]
  # get unstably enterotyped samples
  unstable.samples <-  colnames(data_v12_all)[!colnames(data_v12_all) %in% as.vector(stable.enterotypes.table[,1])] #temprorarily
  # if locat != all we need to delete other locations
  if(locations[i]!='all'){unstable.samples <- unstable.samples[unstable.samples %like% locations[i]]}
  unstable.samples.downsamped <- unstable.samples[!(unstable.samples %in% union(data_v34.poor[[i]],data_v56.poor[[i]]))]
  # unstable.samples_downsamped <- as.vector(stable.enterotypes.table[,1])[!(as.vector(stable.enterotypes.table[,1]) %in% union(data_v34_all.good,data_v56_all.good))]
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=coverages[j],to.append=FALSE,color='black',y_axis_length_factor = 1.5, x_axis_length_factor = 1.03,
                  my.subset=stable.samples.downsamped,adjust.factor=(length(stable.samples.downsamped)+length(unstable.samples.downsamped))/length(stable.samples.downsamped))
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=coverages[j],to.append=TRUE,color='red',my.subset=unstable.samples.downsamped,adjust.factor=(length(stable.samples.downsamped)+length(unstable.samples.downsamped))/length(unstable.samples.downsamped))
}
dev.off()
# V3V4
i=1 #location
jpeg('v3v4_all_enterot_stabil_vs_sample_coverage_window_1k_2k_5k_10k_only_v12_v56_good_cov.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2), mar=c(8,10,7,5), mgp=c(6,2,0))
qiime_file = file.path(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt')
for(j in 1:length(coverages)){
  # From these, we need colnames(data) = samples. A bit long time, but do not care
  data <- readQiimeSmart(input_file=qiime_file, locat = locations[i], percent_threshold = 0,remove.outliers = 0,to.scale=F)
  colnames(data) <- gsub('.V[0-9][0-9]','',colnames(data))
  # get stably enterotyped samples  
  stable.enterotypes.current <- stable_enterotype_file[i] 
  stable.enterotypes.table <- read.table(stable.enterotypes.current)
  stable.samples.downsamped <- as.vector(stable.enterotypes.table[,1])[!(as.vector(stable.enterotypes.table[,1]) %in% union(data_v12.poor[[i]],data_v56.poor[[i]]))]
  # get unstably enterotyped samples
  unstable.samples <-  colnames(data_v34_all)[!colnames(data_v34_all) %in% as.vector(stable.enterotypes.table[,1])] #temprorarily
  # if locat != all we need to delete other locations
  if(locations[i]!='all'){unstable.samples <- unstable.samples[unstable.samples %like% locations[i]]}
  unstable.samples.downsamped <- unstable.samples[!(unstable.samples %in% union(data_v12.poor[[i]],data_v56.poor[[i]]))]
  # unstable.samples_downsamped <- as.vector(stable.enterotypes.table[,1])[!(as.vector(stable.enterotypes.table[,1]) %in% union(data_v34_all.good,data_v56_all.good))]
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=coverages[j],to.append=FALSE,color='black',y_axis_length_factor = 1.5, x_axis_length_factor = 1.03,
                  my.subset=stable.samples.downsamped,adjust.factor=(length(stable.samples.downsamped)+length(unstable.samples.downsamped))/length(stable.samples.downsamped))
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=coverages[j],to.append=TRUE,color='red',my.subset=unstable.samples.downsamped,adjust.factor=(length(stable.samples.downsamped)+length(unstable.samples.downsamped))/length(unstable.samples.downsamped))
}
dev.off()
# V3V4
i=2 #location IL
jpeg('v3v4_IL_enterot_stabil_vs_sample_coverage_window_1k_2k_5k_10k_only_v12_v56_good_cov.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2), mar=c(8,10,7,5), mgp=c(6,2,0))
qiime_file = file.path(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt')
for(j in 1:length(coverages)){
  # From these, we need colnames(data) = samples. A bit long time, but do not care
  data <- readQiimeSmart(input_file=qiime_file, locat = locations[i], percent_threshold = 0,remove.outliers = 0,to.scale=F)
  colnames(data) <- gsub('.V[0-9][0-9]','',colnames(data))
  # get stably enterotyped samples  
  stable.enterotypes.current <- stable_enterotype_file[i] 
  stable.enterotypes.table <- read.table(stable.enterotypes.current)
  stable.samples.downsamped <- as.vector(stable.enterotypes.table[,1])[!(as.vector(stable.enterotypes.table[,1]) %in% union(data_v12.poor[[i]],data_v56.poor[[i]]))]
  # get unstably enterotyped samples
  unstable.samples <-  colnames(data_v34_all)[!colnames(data_v34_all) %in% as.vector(stable.enterotypes.table[,1])] #temprorarily
  # if locat != all we need to delete other locations
  if(locations[i]!='all'){unstable.samples <- unstable.samples[unstable.samples %like% locations[i]]}
  unstable.samples.downsamped <- unstable.samples[!(unstable.samples %in% union(data_v12.poor[[i]],data_v56.poor[[i]]))]
  # unstable.samples_downsamped <- as.vector(stable.enterotypes.table[,1])[!(as.vector(stable.enterotypes.table[,1]) %in% union(data_v34_all.good,data_v56_all.good))]
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=coverages[j],to.append=FALSE,color='black',y_axis_length_factor = 1.5, x_axis_length_factor = 1.03,
                  my.subset=stable.samples.downsamped,adjust.factor=(length(stable.samples.downsamped)+length(unstable.samples.downsamped))/length(stable.samples.downsamped))
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=coverages[j],to.append=TRUE,color='red',my.subset=unstable.samples.downsamped,adjust.factor=(length(stable.samples.downsamped)+length(unstable.samples.downsamped))/length(unstable.samples.downsamped))
}
dev.off()
# V3V4
i=3 #location TR
jpeg('v3v4_TR_enterot_stabil_vs_sample_coverage_window_1k_2k_5k_10k_only_v12_v56_good_cov.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2), mar=c(8,10,7,5), mgp=c(6,2,0))
qiime_file = file.path(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt')
for(j in 1:length(coverages)){
  # From these, we need colnames(data) = samples. A bit long time, but do not care
  data <- readQiimeSmart(input_file=qiime_file, locat = locations[i], percent_threshold = 0,remove.outliers = 0,to.scale=F)
  colnames(data) <- gsub('.V[0-9][0-9]','',colnames(data))
  # get stably enterotyped samples  
  stable.enterotypes.current <- stable_enterotype_file[i] 
  stable.enterotypes.table <- read.table(stable.enterotypes.current)
  stable.samples.downsamped <- as.vector(stable.enterotypes.table[,1])[!(as.vector(stable.enterotypes.table[,1]) %in% union(data_v12.poor[[i]],data_v56.poor[[i]]))]
  # get unstably enterotyped samples
  unstable.samples <-  colnames(data_v34_all)[!colnames(data_v34_all) %in% as.vector(stable.enterotypes.table[,1])] #temprorarily
  # if locat != all we need to delete other locations
  if(locations[i]!='all'){unstable.samples <- unstable.samples[unstable.samples %like% locations[i]]}
  unstable.samples.downsamped <- unstable.samples[!(unstable.samples %in% union(data_v12.poor[[i]],data_v56.poor[[i]]))]
  # unstable.samples_downsamped <- as.vector(stable.enterotypes.table[,1])[!(as.vector(stable.enterotypes.table[,1]) %in% union(data_v34_all.good,data_v56_all.good))]
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=coverages[j],to.append=FALSE,color='black',y_axis_length_factor = 1.5, x_axis_length_factor = 1.03,
                  my.subset=stable.samples.downsamped,adjust.factor=(length(stable.samples.downsamped)+length(unstable.samples.downsamped))/length(stable.samples.downsamped))
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=coverages[j],to.append=TRUE,color='red',my.subset=unstable.samples.downsamped,adjust.factor=(length(stable.samples.downsamped)+length(unstable.samples.downsamped))/length(unstable.samples.downsamped))
}
dev.off()
# V3V4
i=4 #location RE
jpeg('v3v4_RE_enterot_stabil_vs_sample_coverage_window_1k_2k_5k_10k_only_v12_v56_good_cov.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2), mar=c(8,10,7,5), mgp=c(6,2,0))
qiime_file = file.path(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt')
for(j in 1:length(coverages)){
  # From these, we need colnames(data) = samples. A bit long time, but do not care
  data <- readQiimeSmart(input_file=qiime_file, locat = locations[i], percent_threshold = 0,remove.outliers = 0,to.scale=F)
  colnames(data) <- gsub('.V[0-9][0-9]','',colnames(data))
  # get stably enterotyped samples  
  stable.enterotypes.current <- stable_enterotype_file[i] 
  stable.enterotypes.table <- read.table(stable.enterotypes.current)
  stable.samples.downsamped <- as.vector(stable.enterotypes.table[,1])[!(as.vector(stable.enterotypes.table[,1]) %in% union(data_v12.poor[[i]],data_v56.poor[[i]]))]
  # get unstably enterotyped samples
  unstable.samples <-  colnames(data_v34_all)[!colnames(data_v34_all) %in% as.vector(stable.enterotypes.table[,1])] #temprorarily
  # if locat != all we need to delete other locations
  if(locations[i]!='all'){unstable.samples <- unstable.samples[unstable.samples %like% locations[i]]}
  unstable.samples.downsamped <- unstable.samples[!(unstable.samples %in% union(data_v12.poor[[i]],data_v56.poor[[i]]))]
  # unstable.samples_downsamped <- as.vector(stable.enterotypes.table[,1])[!(as.vector(stable.enterotypes.table[,1]) %in% union(data_v34_all.good,data_v56_all.good))]
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=coverages[j],to.append=FALSE,color='black',y_axis_length_factor = 1.5, x_axis_length_factor = 1.03,
                  my.subset=stable.samples.downsamped,adjust.factor=(length(stable.samples.downsamped)+length(unstable.samples.downsamped))/length(stable.samples.downsamped))
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=coverages[j],to.append=TRUE,color='red',my.subset=unstable.samples.downsamped,adjust.factor=(length(stable.samples.downsamped)+length(unstable.samples.downsamped))/length(unstable.samples.downsamped))
}
dev.off()
# V5V6
i=1 #location all
jpeg('v5v6_all_enterot_stabil_vs_sample_coverage_window_1k_2k_5k_10k_only_v12_v34_good_cov.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2), mar=c(8,10,7,5), mgp=c(6,2,0))
qiime_file = file.path(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt')
for(j in 1:length(coverages)){
  # From these, we need colnames(data) = samples. A bit long time, but do not care
  data <- readQiimeSmart(input_file=qiime_file, locat = locations[i], percent_threshold = 0,remove.outliers = 0,to.scale=F)
  colnames(data) <- gsub('.V[0-9][0-9]','',colnames(data))
  # get stably enterotyped samples  
  stable.enterotypes.current <- stable_enterotype_file[i] 
  stable.enterotypes.table <- read.table(stable.enterotypes.current)
  stable.samples.downsamped <- as.vector(stable.enterotypes.table[,1])[!(as.vector(stable.enterotypes.table[,1]) %in% union(data_v12.poor[[i]],data_v34.poor[[i]]))]
  # get unstably enterotyped samples
  unstable.samples <-  colnames(data_v56_all)[!colnames(data_v56_all) %in% as.vector(stable.enterotypes.table[,1])] #temprorarily
  # if locat != all we need to delete other locations
  if(locations[i]!='all'){unstable.samples <- unstable.samples[unstable.samples %like% locations[i]]}
  unstable.samples.downsamped <- unstable.samples[!(unstable.samples %in% union(data_v12.poor[[i]],data_v34.poor[[i]]))]
  # unstable.samples_downsamped <- as.vector(stable.enterotypes.table[,1])[!(as.vector(stable.enterotypes.table[,1]) %in% union(data_v34_all.good,data_v56_all.good))]
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=coverages[j],to.append=FALSE,color='black',y_axis_length_factor = 1.5, x_axis_length_factor = 1.03,
                  my.subset=stable.samples.downsamped,adjust.factor=(length(stable.samples.downsamped)+length(unstable.samples.downsamped))/length(stable.samples.downsamped))
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=coverages[j],to.append=TRUE,color='red',my.subset=unstable.samples.downsamped,adjust.factor=(length(stable.samples.downsamped)+length(unstable.samples.downsamped))/length(unstable.samples.downsamped))
}
dev.off()
# V5V6
i=2 #location IL
jpeg('v5v6_IL_enterot_stabil_vs_sample_coverage_window_1k_2k_5k_10k_only_v12_v34_good_cov.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2), mar=c(8,10,7,5), mgp=c(6,2,0))
qiime_file = file.path(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt')
for(j in 1:length(coverages)){
  # From these, we need colnames(data) = samples. A bit long time, but do not care
  data <- readQiimeSmart(input_file=qiime_file, locat = locations[i], percent_threshold = 0,remove.outliers = 0,to.scale=F)
  colnames(data) <- gsub('.V[0-9][0-9]','',colnames(data))
  # get stably enterotyped samples  
  stable.enterotypes.current <- stable_enterotype_file[i] 
  stable.enterotypes.table <- read.table(stable.enterotypes.current)
  stable.samples.downsamped <- as.vector(stable.enterotypes.table[,1])[!(as.vector(stable.enterotypes.table[,1]) %in% union(data_v12.poor[[i]],data_v34.poor[[i]]))]
  # get unstably enterotyped samples
  unstable.samples <-  colnames(data_v56_all)[!colnames(data_v56_all) %in% as.vector(stable.enterotypes.table[,1])] #temprorarily
  # if locat != all we need to delete other locations
  if(locations[i]!='all'){unstable.samples <- unstable.samples[unstable.samples %like% locations[i]]}
  unstable.samples.downsamped <- unstable.samples[!(unstable.samples %in% union(data_v12.poor[[i]],data_v34.poor[[i]]))]
  # unstable.samples_downsamped <- as.vector(stable.enterotypes.table[,1])[!(as.vector(stable.enterotypes.table[,1]) %in% union(data_v34_all.good,data_v56_all.good))]
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=coverages[j],to.append=FALSE,color='black',y_axis_length_factor = 1.5, x_axis_length_factor = 1.03,
                  my.subset=stable.samples.downsamped,adjust.factor=(length(stable.samples.downsamped)+length(unstable.samples.downsamped))/length(stable.samples.downsamped))
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=coverages[j],to.append=TRUE,color='red',my.subset=unstable.samples.downsamped,adjust.factor=(length(stable.samples.downsamped)+length(unstable.samples.downsamped))/length(unstable.samples.downsamped))
}
dev.off()
# V5V6
i=3 #location TR
jpeg('v5v6_TR_enterot_stabil_vs_sample_coverage_window_1k_2k_5k_10k_only_v12_v34_good_cov.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2), mar=c(8,10,7,5), mgp=c(6,2,0))
qiime_file = file.path(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt')
for(j in 1:length(coverages)){
  # From these, we need colnames(data) = samples. A bit long time, but do not care
  data <- readQiimeSmart(input_file=qiime_file, locat = locations[i], percent_threshold = 0,remove.outliers = 0,to.scale=F)
  colnames(data) <- gsub('.V[0-9][0-9]','',colnames(data))
  # get stably enterotyped samples  
  stable.enterotypes.current <- stable_enterotype_file[i] 
  stable.enterotypes.table <- read.table(stable.enterotypes.current)
  stable.samples.downsamped <- as.vector(stable.enterotypes.table[,1])[!(as.vector(stable.enterotypes.table[,1]) %in% union(data_v12.poor[[i]],data_v34.poor[[i]]))]
  # get unstably enterotyped samples
  unstable.samples <-  colnames(data_v56_all)[!colnames(data_v56_all) %in% as.vector(stable.enterotypes.table[,1])] #temprorarily
  # if locat != all we need to delete other locations
  if(locations[i]!='all'){unstable.samples <- unstable.samples[unstable.samples %like% locations[i]]}
  unstable.samples.downsamped <- unstable.samples[!(unstable.samples %in% union(data_v12.poor[[i]],data_v34.poor[[i]]))]
  # unstable.samples_downsamped <- as.vector(stable.enterotypes.table[,1])[!(as.vector(stable.enterotypes.table[,1]) %in% union(data_v34_all.good,data_v56_all.good))]
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=coverages[j],to.append=FALSE,color='black',y_axis_length_factor = 1.5, x_axis_length_factor = 1.03,
                  my.subset=stable.samples.downsamped,adjust.factor=(length(stable.samples.downsamped)+length(unstable.samples.downsamped))/length(stable.samples.downsamped))
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=coverages[j],to.append=TRUE,color='red',my.subset=unstable.samples.downsamped,adjust.factor=(length(stable.samples.downsamped)+length(unstable.samples.downsamped))/length(unstable.samples.downsamped))
}
dev.off()
# V5V6
i=4 #location RE
jpeg('v5v6_RE_enterot_stabil_vs_sample_coverage_window_1k_2k_5k_10k_only_v12_v34_good_cov.jpg',width = 1600,height = 1600)
par(mfrow=c(2,2), mar=c(8,10,7,5), mgp=c(6,2,0))
qiime_file = file.path(home.dir,'/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt')
for(j in 1:length(coverages)){
  # From these, we need colnames(data) = samples. A bit long time, but do not care
  data <- readQiimeSmart(input_file=qiime_file, locat = locations[i], percent_threshold = 0,remove.outliers = 0,to.scale=F)
  colnames(data) <- gsub('.V[0-9][0-9]','',colnames(data))
  # get stably enterotyped samples  
  stable.enterotypes.current <- stable_enterotype_file[i] 
  stable.enterotypes.table <- read.table(stable.enterotypes.current)
  stable.samples.downsamped <- as.vector(stable.enterotypes.table[,1])[!(as.vector(stable.enterotypes.table[,1]) %in% union(data_v12.poor[[i]],data_v34.poor[[i]]))]
  # get unstably enterotyped samples
  unstable.samples <-  colnames(data_v56_all)[!colnames(data_v56_all) %in% as.vector(stable.enterotypes.table[,1])] #temprorarily
  # if locat != all we need to delete other locations
  if(locations[i]!='all'){unstable.samples <- unstable.samples[unstable.samples %like% locations[i]]}
  unstable.samples.downsamped <- unstable.samples[!(unstable.samples %in% union(data_v12.poor[[i]],data_v34.poor[[i]]))]
  # unstable.samples_downsamped <- as.vector(stable.enterotypes.table[,1])[!(as.vector(stable.enterotypes.table[,1]) %in% union(data_v34_all.good,data_v56_all.good))]
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=coverages[j],to.append=FALSE,color='black',y_axis_length_factor = 1.5, x_axis_length_factor = 1.03,
                  my.subset=stable.samples.downsamped,adjust.factor=(length(stable.samples.downsamped)+length(unstable.samples.downsamped))/length(stable.samples.downsamped))
  exploreCoverage(input_type='QIIME', input_file=qiime_file, locat = locations[i], percent_threshold = 0, remove.outliers = 0, cumulative = FALSE,
                  usredneniye=coverages[j],to.append=TRUE,color='red',my.subset=unstable.samples.downsamped,adjust.factor=(length(stable.samples.downsamped)+length(unstable.samples.downsamped))/length(unstable.samples.downsamped))
}
dev.off()


#6. get v1v2, v3v4 and v5v6 samples with poor and good coverage - iteration 2
# v1v2
data_v12_all <- readQiimeSmart(input_file = qiime_file_v12, locat = 'all', percent_threshold = 0, remove.outliers = 0,to.scale = FALSE)
colnames(data_v12_all) <- gsub('.V[0-9][0-9]','',colnames(data_v12_all))
data_v12_all.poor <- colnames(data_v12_all[,colSums(data_v12_all)<19000])
data_v12_all.good <- colnames(data_v12_all[,colSums(data_v12_all)>=19000])
data_v12_IL <- readQiimeSmart(input_file = qiime_file_v12, locat = '.IL', percent_threshold = 0, remove.outliers = 0,to.scale = FALSE)
colnames(data_v12_IL) <- gsub('.V[0-9][0-9]','',colnames(data_v12_IL))
data_v12_IL.poor <- colnames(data_v12_IL[,colSums(data_v12_IL)<19000])
data_v12_IL.good <- colnames(data_v12_IL[,colSums(data_v12_IL)>=19000])
data_v12_TR <- readQiimeSmart(input_file = qiime_file_v12, locat = '.TR', percent_threshold = 0, remove.outliers = 0,to.scale = FALSE)
colnames(data_v12_TR) <- gsub('.V[0-9][0-9]','',colnames(data_v12_TR))
data_v12_TR.poor <- colnames(data_v12_TR[,colSums(data_v12_TR)<24000])
data_v12_TR.good <- colnames(data_v12_TR[,colSums(data_v12_TR)>=24000])
data_v12_RE <- readQiimeSmart(input_file = qiime_file_v12, locat = '.RE', percent_threshold = 0, remove.outliers = 0,to.scale = FALSE)
colnames(data_v12_RE) <- gsub('.V[0-9][0-9]','',colnames(data_v12_RE))
data_v12_RE.poor <- colnames(data_v12_RE[,colSums(data_v12_RE)<23000])
data_v12_RE.good <- colnames(data_v12_RE[,colSums(data_v12_RE)>=23000])
data_v12.poor <- list(data_v12_all.poor,data_v12_IL.poor,data_v12_TR.poor,data_v12_RE.poor)
data_v12.good <- list(data_v12_all.good,data_v12_IL.good,data_v12_TR.good,data_v12_RE.good)
# v3v4
data_v34_all <- readQiimeSmart(input_file = qiime_file_v34, locat = 'all', percent_threshold = 0, remove.outliers = 0,to.scale = FALSE)
colnames(data_v34_all) <- gsub('.V[0-9][0-9]','',colnames(data_v34_all))
data_v34_all.poor <- colnames(data_v34_all[,colSums(data_v34_all)<10000])
data_v34_all.good <- colnames(data_v34_all[,colSums(data_v34_all)>=10000])
data_v34_IL <- readQiimeSmart(input_file = qiime_file_v34, locat = '.IL', percent_threshold = 0, remove.outliers = 0,to.scale = FALSE)
colnames(data_v34_IL) <- gsub('.V[0-9][0-9]','',colnames(data_v34_IL))
data_v34_IL.poor <- colnames(data_v34_IL[,colSums(data_v34_IL)<12000])
data_v34_IL.good <- colnames(data_v34_IL[,colSums(data_v34_IL)>=12000])
data_v34_TR <- readQiimeSmart(input_file = qiime_file_v34, locat = '.TR', percent_threshold = 0, remove.outliers = 0,to.scale = FALSE)
colnames(data_v34_TR) <- gsub('.V[0-9][0-9]','',colnames(data_v34_TR))
data_v34_TR.poor <- colnames(data_v34_TR[,colSums(data_v34_TR)<15000])
data_v34_TR.good <- colnames(data_v34_TR[,colSums(data_v34_TR)>=15000])
data_v34_RE <- readQiimeSmart(input_file = qiime_file_v34, locat = '.RE', percent_threshold = 0, remove.outliers = 0,to.scale = FALSE)
colnames(data_v34_RE) <- gsub('.V[0-9][0-9]','',colnames(data_v34_RE))
data_v34_RE.poor <- colnames(data_v34_RE[,colSums(data_v34_RE)<14000])
data_v34_RE.good <- colnames(data_v34_RE[,colSums(data_v34_RE)>=14000])
data_v34.poor <- list(data_v34_all.poor,data_v34_IL.poor,data_v34_TR.poor,data_v34_RE.poor)
data_v34.good <- list(data_v34_all.good,data_v34_IL.good,data_v34_TR.good,data_v34_RE.good)
# v5v6
data_v56_all <- readQiimeSmart(input_file = qiime_file_v56, locat = 'all', percent_threshold = 0, remove.outliers = 0,to.scale = FALSE)
colnames(data_v56_all) <- gsub('.V[0-9][0-9]','',colnames(data_v56_all))
data_v56_all.poor <- colnames(data_v56_all[,colSums(data_v56_all)<20000])
data_v56_all.good <- colnames(data_v56_all[,colSums(data_v56_all)>=20000])
data_v56_IL <- readQiimeSmart(input_file = qiime_file_v56, locat = '.IL', percent_threshold = 0, remove.outliers = 0,to.scale = FALSE)
colnames(data_v56_IL) <- gsub('.V[0-9][0-9]','',colnames(data_v56_IL))
data_v56_IL.poor <- colnames(data_v56_IL[,colSums(data_v56_IL)<25000])
data_v56_IL.good <- colnames(data_v56_IL[,colSums(data_v56_IL)>=25000])
data_v56_TR <- readQiimeSmart(input_file = qiime_file_v56, locat = '.TR', percent_threshold = 0, remove.outliers = 0,to.scale = FALSE)
colnames(data_v56_TR) <- gsub('.V[0-9][0-9]','',colnames(data_v56_TR))
data_v56_TR.poor <- colnames(data_v56_TR[,colSums(data_v56_TR)<20000])
data_v56_TR.good <- colnames(data_v56_TR[,colSums(data_v56_TR)>=20000])
data_v56_RE <- readQiimeSmart(input_file = qiime_file_v56, locat = '.RE', percent_threshold = 0, remove.outliers = 0,to.scale = FALSE)
colnames(data_v56_RE) <- gsub('.V[0-9][0-9]','',colnames(data_v56_RE))
data_v56_RE.poor <- colnames(data_v56_RE[,colSums(data_v56_RE)<22000])
data_v56_RE.good <- colnames(data_v56_RE[,colSums(data_v56_RE)>=22000])
data_v56.poor <- list(data_v56_all.poor,data_v56_IL.poor,data_v56_TR.poor,data_v56_RE.poor)
data_v56.good <- list(data_v56_all.good,data_v56_IL.good,data_v56_TR.good,data_v56_RE.good)

# these are samples that have good coverage for all amplicons
good.cov.all <- intersect(intersect(data_v12.good[[1]],data_v34.good[[1]]),data_v56.good[[1]])
good.cov.IL <- intersect(intersect(data_v12.good[[2]],data_v34.good[[2]]),data_v56.good[[2]])
good.cov.TR <- intersect(intersect(data_v12.good[[3]],data_v34.good[[3]]),data_v56.good[[3]])
good.cov.RE <- intersect(intersect(data_v12.good[[4]],data_v34.good[[4]]),data_v56.good[[4]])
good.cov <- list(good.cov.all,good.cov.IL,good.cov.TR,good.cov.RE)
# these are samples that have poor coverage for all amplicons
poor.cov.all <- union(union(data_v12.poor[[1]],data_v34.poor[[1]]),data_v56.poor[[1]])
poor.cov.IL <- union(union(data_v12.poor[[2]],data_v34.poor[[2]]),data_v56.poor[[2]])
poor.cov.TR <- union(union(data_v12.poor[[3]],data_v34.poor[[3]]),data_v56.poor[[3]])
poor.cov.RE <- union(union(data_v12.poor[[4]],data_v34.poor[[4]]),data_v56.poor[[4]])
poor.cov <- list(poor.cov.all,poor.cov.IL,poor.cov.TR,poor.cov.RE)

# now we can understand which percentage of samples that have good coverage for all amplicons are stable across all amplicons
locations <- c('all','.IL','.TR','.RE')
stable_enterotype_file <- c(file.path(home.dir,'/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype/Stable_enterotypes.csv'),
                            file.path(home.dir,'/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/Old-the-same/Stable_enterotypes_IL.csv'),
                            file.path(home.dir,'/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/Old-the-same/Stable_enterotypes_TR.csv'),
                            file.path(home.dir,'/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/Old-the-same/Stable_enterotypes_RE.csv'))
stability.in.good.cov.samples <- numeric(4)
for (i in 1:length(locations)){
stable.enterotypes.current <- stable_enterotype_file[i] 
stable.enterotypes.table <- read.table(stable.enterotypes.current)
# intersection of well covered with stable
stable.samples <- as.vector(stable.enterotypes.table[,1])
stable.samples.intersect.good.cov <- intersect(stable.samples,good.cov[[i]])
# intersection of well covered with unstable
unstable.samples <-  colnames(data_v56_all)[!colnames(data_v56_all) %in% stable.samples] 
if(i!=1){unstable.samples <- unstable.samples[unstable.samples %like% locations[i]]}
unstable.samples.intersect.good.cov <- intersect(unstable.samples,good.cov[[i]])
# get fraction of covered samples that are stable
stability.in.good.cov.samples[i] <- signif(length(stable.samples.intersect.good.cov)/(length(stable.samples.intersect.good.cov)+length(unstable.samples.intersect.good.cov)),digits = 3)
print(paste0('Stability of ',locations[i],' samples with coverage >= threshold: ',100*stability.in.good.cov.samples[i],'%, N=',
             length(good.cov[[i]])))
}

# now - stability for samples with poor cov at least for 1 amplicon
stability.in.poor.cov.samples <- numeric(4)
for (i in 1:length(locations)){
  stable.enterotypes.current <- stable_enterotype_file[i] 
  stable.enterotypes.table <- read.table(stable.enterotypes.current)
  # intersection of well covered with stable
  stable.samples <- as.vector(stable.enterotypes.table[,1])
  stable.samples.intersect.poor.cov <- intersect(stable.samples,poor.cov[[i]])
  # intersection of well covered with unstable
  unstable.samples <-  colnames(data_v56_all)[!colnames(data_v56_all) %in% stable.samples] 
  if(i!=1){unstable.samples <- unstable.samples[unstable.samples %like% locations[i]]}
  unstable.samples.intersect.poor.cov <- intersect(unstable.samples,poor.cov[[i]])
  # get fraction of covered samples that are stable
  stability.in.poor.cov.samples[i] <- signif(length(stable.samples.intersect.poor.cov)/(length(stable.samples.intersect.poor.cov)+length(unstable.samples.intersect.poor.cov)),digits = 3)
  print(paste0('Stability of ',locations[i],' samples with coverage < threshold: ',100*stability.in.poor.cov.samples[i],'%, N=',
               length(poor.cov[[i]])))
}

# Average enterotype stability per location without coverage filtration
stab.location <- numeric(length(locations))
for (i in 1:length(locations)){
  stab.location[i] <- signif((stability.in.good.cov.samples[i]*length(good.cov[[i]])+stability.in.poor.cov.samples[i]*length(poor.cov[[i]]))/(length(poor.cov[[i]])+length(good.cov[[i]])),
                             digits = 3)
print(paste0('Average stability for ',locations[i],': ',stab.location[i]))
}

# 7. Try to change cov thresholds for IL and see the results
# cutpoint set I. 
data_v12_IL <- readQiimeSmart(input_file = qiime_file_v12, locat = '.IL', percent_threshold = 0, remove.outliers = 0,to.scale = FALSE)
colnames(data_v12_IL) <- gsub('.V[0-9][0-9]','',colnames(data_v12_IL))
data_v12_IL.poor <- colnames(data_v12_IL[,colSums(data_v12_IL)<34000])
data_v12_IL.good <- colnames(data_v12_IL[,colSums(data_v12_IL)>=34000])
data_v12.poor <- list(data_v12_all.poor,data_v12_IL.poor,data_v12_TR.poor,data_v12_RE.poor)
data_v12.good <- list(data_v12_all.good,data_v12_IL.good,data_v12_TR.good,data_v12_RE.good)


data_v34_IL <- readQiimeSmart(input_file = qiime_file_v34, locat = '.IL', percent_threshold = 0, remove.outliers = 0,to.scale = FALSE)
colnames(data_v34_IL) <- gsub('.V[0-9][0-9]','',colnames(data_v34_IL))
data_v34_IL.poor <- colnames(data_v34_IL[,colSums(data_v34_IL)<20000])
data_v34_IL.good <- colnames(data_v34_IL[,colSums(data_v34_IL)>=20000])
data_v34_TR <- readQiimeSmart(input_file = qiime_file_v34, locat = '.TR', percent_threshold = 0, remove.outliers = 0,to.scale = FALSE)
data_v34.poor <- list(data_v34_all.poor,data_v34_IL.poor,data_v34_TR.poor,data_v34_RE.poor)
data_v34.good <- list(data_v34_all.good,data_v34_IL.good,data_v34_TR.good,data_v34_RE.good)

data_v56_IL <- readQiimeSmart(input_file = qiime_file_v56, locat = '.IL', percent_threshold = 0, remove.outliers = 0,to.scale = FALSE)
colnames(data_v56_IL) <- gsub('.V[0-9][0-9]','',colnames(data_v56_IL))
data_v56_IL.poor <- colnames(data_v56_IL[,colSums(data_v56_IL)<39000])
data_v56_IL.good <- colnames(data_v56_IL[,colSums(data_v56_IL)>=39000])
data_v56.poor <- list(data_v56_all.poor,data_v56_IL.poor,data_v56_TR.poor,data_v56_RE.poor)
data_v56.good <- list(data_v56_all.good,data_v56_IL.good,data_v56_TR.good,data_v56_RE.good)

# cutpoint set II.
data_v12_IL <- readQiimeSmart(input_file = qiime_file_v12, locat = '.IL', percent_threshold = 0, remove.outliers = 0,to.scale = FALSE)
colnames(data_v12_IL) <- gsub('.V[0-9][0-9]','',colnames(data_v12_IL))
data_v12_IL.poor <- colnames(data_v12_IL[,colSums(data_v12_IL)<19000] | data_v12_IL[,colSums(data_v12_IL)>85000])
data_v12_IL.good <- colnames(data_v12_IL[,colSums(data_v12_IL)>=19000] & data_v12_IL[,colSums(data_v12_IL)<=85000])
data_v12.poor <- list(data_v12_all.poor,data_v12_IL.poor,data_v12_TR.poor,data_v12_RE.poor)
data_v12.good <- list(data_v12_all.good,data_v12_IL.good,data_v12_TR.good,data_v12_RE.good)


data_v34_IL <- readQiimeSmart(input_file = qiime_file_v34, locat = '.IL', percent_threshold = 0, remove.outliers = 0,to.scale = FALSE)
colnames(data_v34_IL) <- gsub('.V[0-9][0-9]','',colnames(data_v34_IL))
data_v34_IL.poor <- colnames(data_v34_IL[,colSums(data_v34_IL)<12000] | colSums(data_v34_IL)>83000)
data_v34_IL.good <- colnames(data_v34_IL[,colSums(data_v34_IL)>=12000] & colSums(data_v34_IL)<=83000)
data_v34_TR <- readQiimeSmart(input_file = qiime_file_v34, locat = '.TR', percent_threshold = 0, remove.outliers = 0,to.scale = FALSE)
data_v34.poor <- list(data_v34_all.poor,data_v34_IL.poor,data_v34_TR.poor,data_v34_RE.poor)
data_v34.good <- list(data_v34_all.good,data_v34_IL.good,data_v34_TR.good,data_v34_RE.good)

data_v56_IL <- readQiimeSmart(input_file = qiime_file_v56, locat = '.IL', percent_threshold = 0, remove.outliers = 0,to.scale = FALSE)
colnames(data_v56_IL) <- gsub('.V[0-9][0-9]','',colnames(data_v56_IL))
data_v56_IL.poor <- colnames(data_v56_IL[,colSums(data_v56_IL)<25000])
data_v56_IL.good <- colnames(data_v56_IL[,colSums(data_v56_IL)>=25000])
data_v56.poor <- list(data_v56_all.poor,data_v56_IL.poor,data_v56_TR.poor,data_v56_RE.poor)
data_v56.good <- list(data_v56_all.good,data_v56_IL.good,data_v56_TR.good,data_v56_RE.good)
