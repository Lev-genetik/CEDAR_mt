#needs VPN to S01 server activated
#constructs a matrix (data table? data frame?) ent.combined with all enterotypes for each sample
require('data.table')
require('stringr')
require('pracma')
source('src/functions/metagenome_func.R')
source('src/functions/read_abundance_data_v1v2.R')
source('src/Lev/compare_enterotype_predictions.R')
ent12 <- GetEnterotypesCustom(data.multifactor$g$table, plot.indexes = F, k = 2)
# source('src/enterotypes/explore_enterotypes.R')
print('v1v2 enterotyping finished')

source('src/functions/read_abundance_data.R')
ent34 <- GetEnterotypesCustom(data.multifactor$g$table, plot.indexes = F, k = 2)
print('v3v4 enterotyping finished')

source('src/functions/read_abundance_data_v5v6.R')
ent56 <- GetEnterotypesCustom(data.multifactor$g$table, plot.indexes = F, k = 2)
print('v5v6 enterotyping finished')

CompEntPred(ent12,ent34,ent56)
# ent.combined <- cbind(ent12.data.sample.enterotype,ent34.data.sample.enterotype,ent56.data.sample.enterotype)