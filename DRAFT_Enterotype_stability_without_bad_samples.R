# First strange attempt to solve the question of enterotype stability when removing bad samples
prob.v56 <- perm.v56[perm.v56[,3]==1,]
dim(prob.v56)
prob.v56[,1]
stable.samples <- '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype/Stable_enterotypes.csv'
stable.samples <- read.table('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype/Stable_enterotypes.csv')
head(stable.samples[,1])
# perform selection of v5v6 samples with stabil > thershold
good.prob.56 <- prob.v56[prob.v56[,6]>=0.976,]
good.prob.56[,1] <- gsub('\\.V[0-9][0-9]','',good.prob.56[,1])
good.prob.56.stable <- good.prob.56 %in% stable.samples[,1]



