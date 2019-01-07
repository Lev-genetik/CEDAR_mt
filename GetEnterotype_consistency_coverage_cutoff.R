# Aim: measure eneterotype stability when applying certain coverage threshold for each location-enterotype pair
# INPUT:
# locat = location (all, .IL, .TR or .RE)
# stable_enterotype_file = e.g. '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/Old-the-same/Stable_enterotypes_IL.csv'
# thresholds = vector of thresholds for amplicons (v1v2,v3v4,v5v6) - po umolch (10 000,10 000, 10 000)
# qiime_files = QIIME output files
# percent_threshold - for readQiimeSmart noise removal of rare species, we set = 0 as want to see total coverage 
# remove.outliers - for readQiimeSmart, also = 0 for not removing any samples (anyway, that does not affect enterotyping! Keep it = 0)
# OUTPUT:
# a data frame with fraction and total number of stable enterotypes among well covered samples, other samples, all samples
getEnterotypeConsistCoverCutoff <- function(locat, stable_enterotype_file, thresholds = c(10000,10000,10000), qiime_files = c(
  '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt', 
  '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt', 
  '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt'),
                                            percent_threshold = 0, remove.outliers = 0){
home.dir <- '/media/lev-genetik/980E73270E72FD96/Liege'
progr.dir <- '/home/lev-genetik/Desktop/Projects/liege/src/Lev'
source(file.path(progr.dir,'/Functions_general_non_enterotyping.R'))
# V12
data_v12 <- readQiimeSmart(input_file = qiime_files[1], locat = locat, percent_threshold = percent_threshold, 
                           remove.outliers = remove.outliers,to.scale = FALSE)
colnames(data_v12) <- gsub('.V[0-9][0-9]','',colnames(data_v12))
data_v12.poor <- colnames(data_v12[,colSums(data_v12)<thresholds[1]])
data_v12.good <- colnames(data_v12[,colSums(data_v12)>=thresholds[1]])
# V34
data_v34 <- readQiimeSmart(input_file = qiime_files[2], locat = locat, percent_threshold = percent_threshold, 
                           remove.outliers = remove.outliers,to.scale = FALSE)
colnames(data_v34) <- gsub('.V[0-9][0-9]','',colnames(data_v34))
data_v34.poor <- colnames(data_v34[,colSums(data_v34)<thresholds[2]])
data_v34.good <- colnames(data_v34[,colSums(data_v34)>=thresholds[2]])
# V56
data_v56 <- readQiimeSmart(input_file = qiime_files[3], locat = locat, percent_threshold = percent_threshold, 
                           remove.outliers = remove.outliers,to.scale = FALSE)
colnames(data_v56) <- gsub('.V[0-9][0-9]','',colnames(data_v56))
data_v56.poor <- colnames(data_v56[,colSums(data_v56)<thresholds[3]])
data_v56.good <- colnames(data_v56[,colSums(data_v56)>=thresholds[3]])
# these are samples that have good coverage for all amplicons
good.cov <- intersect(intersect(data_v12.good,data_v34.good),data_v56.good)
# these are samples that have poor coverage for at least 1 amplicon
poor.cov <- union(union(data_v12.poor,data_v34.poor),data_v56.poor)

# read the enterotype file
stable.enterotypes.table <- read.table(stable_enterotype_file)

# intersection of well covered with stable
stable.samples <- as.vector(stable.enterotypes.table[,1])
stable.samples.intersect.good.cov <- intersect(stable.samples,good.cov)
# intersection of well covered with unstable
unstable.samples <-  colnames(data_v12)[!colnames(data_v12) %in% stable.samples] 
if(locat!='all'){unstable.samples <- unstable.samples[unstable.samples %like% locat]}
unstable.samples.intersect.good.cov <- intersect(unstable.samples,good.cov)
# intersection of poorly covered with stable and unstable
stable.samples.intersect.poor.cov <- intersect(stable.samples,poor.cov)
unstable.samples.intersect.poor.cov <- intersect(unstable.samples,poor.cov)
# get fraction of covered samples that are stable
stability.in.good.cov.samples <- signif(length(stable.samples.intersect.good.cov)/(length(stable.samples.intersect.good.cov)+length(unstable.samples.intersect.good.cov)),digits = 3)
stability.in.poor.cov.samples <- signif(length(stable.samples.intersect.poor.cov)/(length(stable.samples.intersect.poor.cov)+length(unstable.samples.intersect.poor.cov)),digits = 3)
average.stab <- signif((stability.in.good.cov.samples*length(good.cov)+stability.in.poor.cov.samples*length(poor.cov))/(length(poor.cov)+length(good.cov)),
                       digits = 3)
df_result <- rbind('stability.fraction' = c(stability.in.good.cov.samples, stability.in.poor.cov.samples,average.stab),'N' = c(length(good.cov),length(poor.cov),length(good.cov)+length(poor.cov)))
colnames(df_result) <- c('Good cov','Poor cov','All samples')
return(df_result)
}


