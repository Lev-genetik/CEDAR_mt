# Get stable enterotypes across V12 V34 and V56. All samples enterotypes together + separately by location

# Samples enterotyped separately by location
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Enterotyping_compar_functions_v5.R')
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/GetStablyEnterotypedSamples.R')
setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype')
# read the files with original enterotypes separately calculated for IL, TR and RE
df_1.IL <- read.table ('V1V2_3e_IL_sample-enterotype.csv')
df_2.IL <- read.table ('V3V4_3e_IL_sample-enterotype.csv')
df_3.IL <- read.table ('V5V6_3e_IL_sample-enterotype.csv')
df_1.TR <- read.table ('V1V2_3e_TR_sample-enterotype.csv')
df_2.TR <- read.table ('V3V4_3e_TR_sample-enterotype.csv')
df_3.TR <- read.table ('V5V6_3e_TR_sample-enterotype.csv')
df_1.RE <- read.table ('V1V2_3e_RE_sample-enterotype.csv')
df_2.RE <- read.table ('V3V4_3e_RE_sample-enterotype.csv')
df_3.RE <- read.table ('V5V6_3e_RE_sample-enterotype.csv')

# Here, we check that cluster correspondences of 1-3 are the same as if we combine 1-2 with 2-3

# IL uncomment when needed to design the correspondence matrix
# CompEntPred.pair(df_1.IL,df_3.IL)$Optimal_cluster_correspondences
# CompEntPred.pair(df_1.IL,df_2.IL)$Optimal_cluster_correspondences
# CompEntPred.pair(df_2.IL,df_3.IL)$Optimal_cluster_correspondences

# # TR uncomment when needed to design the correspondence matrix
# CompEntPred.pair(df_1.TR,df_3.TR)$Optimal_cluster_correspondences
# CompEntPred.pair(df_1.TR,df_2.TR)$Optimal_cluster_correspondences
# CompEntPred.pair(df_2.TR,df_3.TR)$Optimal_cluster_correspondences

# # RE uncomment when needed to design the correspondence matrix
# CompEntPred.pair(df_1.RE,df_3.RE)$Optimal_cluster_correspondences
# CompEntPred.pair(df_1.RE,df_2.RE)$Optimal_cluster_correspondences
# CompEntPred.pair(df_2.RE,df_3.RE)$Optimal_cluster_correspondences

# For IL, cluster correspondences are:
correspondence.IL.12.56 <- c('k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Bacteroidaceae.g__Bacteroides','k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Prevotellaceae.g__Prevotella','k__Bacteria.p__Firmicutes.c__Bacilli.o__Bacillales.f__Planococcaceae.g__')
correspondence.IL.12.12 <- c('k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Bacteroidaceae.g__Bacteroides','k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Prevotellaceae.g__Prevotella','k__Bacteria.p__Actinobacteria.c__Actinobacteria.o__Actinomycetales.f__Corynebacteriaceae.g__Corynebacterium')
correspondence.IL.12.34 <- c('k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pasteurellales.f__Pasteurellaceae.g__Aggregatibacter','k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Lachnospiraceae.g__','k__Bacteria.p__Firmicutes.c__Bacilli.o__Gemellales.f__Gemellaceae.g__Gemella')
correspondence.IL <- cbind(correspondence.IL.12.12,correspondence.IL.12.34,correspondence.IL.12.56)

# For TR, cluster correspondences are:
correspondence.TR.12.12 <- c('k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Bacteroidaceae.g__Bacteroides','k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Prevotellaceae.g__Prevotella','k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Alteromonadales.f__Psychromonadaceae.g__Psychromonas')
correspondence.TR.12.34 <- c('k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Bacteroidaceae.g__Bacteroides','k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Prevotellaceae.g__Prevotella','k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Aeromonadales.f__Aeromonadaceae.g__')
correspondence.TR.12.56 <- c('k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Bacteroidaceae.g__Bacteroides','k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Prevotellaceae.g__Prevotella','k__Bacteria.p__Proteobacteria.c__Betaproteobacteria.o__Burkholderiales.f__Comamonadaceae.g__Delftia')
correspondence.TR <- cbind(correspondence.TR.12.12,correspondence.TR.12.34,correspondence.TR.12.56)

# For RE, cluster correspondences are:
correspondence.RE.12.12 <- c('k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Bacteroidaceae.g__Bacteroides','k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Prevotellaceae.g__Prevotella','k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Aeromonadales.f__Succinivibrionaceae.g__Succinivibrio')
correspondence.RE.12.34 <- c('k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Bacteroidaceae.g__Bacteroides','k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Prevotellaceae.g__Prevotella','k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pseudomonadales.f__Moraxellaceae.g__Acinetobacter')
correspondence.RE.12.56 <- c('k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Bacteroidaceae.g__Bacteroides','k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Prevotellaceae.g__Prevotella','k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Pseudomonadales.f__Moraxellaceae.g__Enhydrobacter')
correspondence.RE <- cbind(correspondence.RE.12.12,correspondence.RE.12.34,correspondence.RE.12.56)

result.IL <- getStableSamples(df1=df_1.IL,df2=df_2.IL,df3=df_3.IL,correspondence.IL)
result.TR <- getStableSamples(df1=df_1.TR,df2=df_2.TR,df3=df_3.TR,correspondence.TR)
result.RE <- getStableSamples(df1=df_1.RE,df2=df_2.RE,df3=df_3.RE,correspondence.RE)

# write the results to files
setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype')
write.table(result.IL,'Stable_enterotypes_IL.csv')
write.table(result.TR,'Stable_enterotypes_TR.csv')
write.table(result.RE,'Stable_enterotypes_RE.csv')

# Get stable enterotypes v12-v34-v56 for all locations enterotyped together
# Done without a script

# Get stable enterotypes across V12 and V56 only for all locations enterotyped together
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/GetStablyEnterotypedSamples.R')
correspondence.without.34.12 <- c("k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Bacteroidaceae.g__Bacteroides","k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Prevotellaceae.g__Prevotella","k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Oceanospirillales.f__Halomonadaceae.g__")
correspondence.without.34.56 <- c("k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Bacteroidaceae.g__Bacteroides","k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Prevotellaceae.g__Prevotella","k__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Alteromonadales.f__Shewanellaceae.g__Shewanella")
correspondence.without.34 <- cbind(correspondence.without.34.12,correspondence.without.34.12,correspondence.without.34.56)
df_1 <- read.table ('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype/V1V2_3e_sample-enterotype.csv')
df_3 <- read.table ('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype/V5V6_3e_sample-enterotype.csv')
result <- getStableSamples(df1=df_1,df2=df_1,df3=df_3,correspondence.without.34)
write.table(result,'/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype/Stable_enterotypes_v12_v56.csv')