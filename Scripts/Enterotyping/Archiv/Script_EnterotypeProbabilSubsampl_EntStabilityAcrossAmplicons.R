# Script to get stable enterotypes for each amplicon by subsampling. Get enterotypes for half of the samples
# and asssign stable enterotypes only to samples that have the same enterotype in >=90% of cases
# all locations enterotyped together
# for v1v2, v3v4 and v5v6
# Than, draw graphs of the enterotype probability density
# Than, for each amplicon, sort samples according to enterotype probability in descending order and check
# if the probability of having stable enterotype across amplicons also decreases.

setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities')
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Enterotype_Probabilities_functions.R')
# v1v2 probability enterotyping
sample.enterotype.probabil.subsampl.v12 <- getStablEntSubsampl(input_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt',N=1000)
write.table(sample.enterotype.probabil.subsampl.v12,'V1V2_3e_sample-enterotype_probabilities.csv')
# v3v4 probability enterotyping
sample.enterotype.probabil.subsampl.v34 <- getStablEntSubsampl(input_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt',N=1000)
write.table(sample.enterotype.probabil.subsampl.v34,'V3V4_3e_sample-enterotype_probabilities.csv')
# v5v6 probability enterotyping
sample.enterotype.probabil.subsampl.v56 <- getStablEntSubsampl(input_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt',N=1000)
write.table(sample.enterotype.probabil.subsampl.v56,'V5V6_3e_sample-enterotype_probabilities.csv')

# now draw graphs
library(graphics)
sample.enterotype.probabil.subsampl.v12 <- read.table('V1V2_3e_sample-enterotype_probabilities.csv')
sample.enterotype.probabil.subsampl.v34 <- read.table('V3V4_3e_sample-enterotype_probabilities.csv')
sample.enterotype.probabil.subsampl.v56 <- read.table('V5V6_3e_sample-enterotype_probabilities.csv')
setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities/Graphs_enterotype_probability')

jpeg(filename = 'Ent_probabil_density_12_black_34_red_56_blue_095_100.jpg',width=1000,height = 600)
plot(density(sample.enterotype.probabil.subsampl.v12$Probability),xlim = c(0.95,1),xlab = 'Enterotype probability',main = '')
lines(density(sample.enterotype.probabil.subsampl.v34$Probability), col='red')
lines(density(sample.enterotype.probabil.subsampl.v56$Probability), col='blue')
dev.off()

jpeg(filename = 'Ent_probabil_density_12_black_34_red_56_blue_050_100.jpg',width=1000,height = 600)
plot(density(sample.enterotype.probabil.subsampl.v12$Probability),xlim = c(0.5,1),xlab = 'Enterotype probability',main = '')
lines(density(sample.enterotype.probabil.subsampl.v34$Probability), col='red')
lines(density(sample.enterotype.probabil.subsampl.v56$Probability), col='blue')
dev.off()

jpeg(filename = 'Ent_probabil_density_v1v2_0_1.jpg',width=1000,height = 600)
plot(density(sample.enterotype.probabil.subsampl.v12$Probability),xlim = c(0,1),xlab = 'v1v2 enterotype probability',col='black',main = '')
dev.off()
jpeg(filename = 'Ent_probabil_density_v3v4_0_1.jpg',width=1000,height = 600)
plot(density(sample.enterotype.probabil.subsampl.v34$Probability),xlim = c(0,1),xlab = 'v3v4 enterotype probability',col='red',main = '')
dev.off()
jpeg(filename = 'Ent_probabil_density_v5v6_0_1.jpg',width=1000,height = 600)
plot(density(sample.enterotype.probabil.subsampl.v56$Probability),xlim = c(0,1),xlab = 'v5v6 enterotype probability',col='blue',main = '')
dev.off()


# Than, for each amplicon, sort samples according to enterotype probability in descending order and check
# if the probability of having stable enterotype across amplicons also decreases.
library(data.table)
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Enterotype_Probabilities_functions.R')
setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities/Graphs_enterotype_probability')

# For V12
ent.prob.vs.stability.v12 <- subsVsLoc(file1 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities/V1V2_3e_sample-enterotype_probabilities.csv',file2='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype/Stable_enterotypes.csv', folder = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities/Graphs_enterotype_probability',graph.file = 'v1v2_ent_probability_by_sampl_vs_v12v34v56.jpg',window=200)
# get the first number of sample where probability of enterotype is already <1.000   (n)
sample.enterotype.probabil.subsampl.v12 <- read.table('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities/V1V2_3e_sample-enterotype_probabilities.csv')
sample.enterotype.probabil.subsampl.v12.sorted <- sample.enterotype.probabil.subsampl.v12[order(-sample.enterotype.probabil.subsampl.v12$Probability),]
for (i in 1:dim(sample.enterotype.probabil.subsampl.v12.sorted)[1]){
  if(sample.enterotype.probabil.subsampl.v12.sorted$Probability[i]==1){n=i}else{return()}
}
n=n+1
jpeg('Enterotype_probabil_v12_descend_cumulat_vs_ent_stability_v12_v34_v56.jpg',width = 1000, height = 600)
plot(ent.prob.vs.stability.v12$ent.stability.besttoworst.sum.first.n.samples[n:dim(sample.enterotype.probabil.subsampl.v12.sorted)[1]], ylim=c(0.5,0.75), xlab='Sampl sort enterotype probabil in descend, prob 100% excl', ylab='Probability of stable enterotype')
dev.off()

# get enterotype probability vs v1v2 v5v6 stability graph (window for counting average consistency=200)
ent.prob.vs.stability_v12v56.v12 <- subsVsLoc(file1 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities/V1V2_3e_sample-enterotype_probabilities.csv',file2='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype/Stable_enterotypes_v12_v56.csv', folder = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities/Graphs_enterotype_probability',graph.file = 'v1v2_ent_probability_by_sampl_vs_v12v56.jpg',window=200)
jpeg('Enterotype_probabil_v12_descend_cumulat_vs_ent_stability_v12_v56.jpg',width = 1000, height = 600)
plot(ent.prob.vs.stability_v12v56.v12$ent.stability.besttoworst.sum.first.n.samples[n:dim(sample.enterotype.probabil.subsampl.v12.sorted)[1]], ylim=c(0.6,0.8), xlab='Sampl sort enterotype probabil in descend, prob 100% excl', ylab='Probability of stable enterotype')
dev.off()

# to manually see after which sample, the enterotype stability decreses. 
order(-ent.prob.vs.stability_v12v56.v12$ent.stability.besttoworst.sum.first.n.samples[n:dim(sample.enterotype.probabil.subsampl.v12.sorted)[1]])
# Here, n=530+109=639
# For first 639 samples, v12-v56 stable enterotype found in 74.6%.
# For last 148 samples, v12-v56 stable enterotype found in 56.1%.
mean(ent.prob.vs.stability_v12v56.v12$identity.vector[1:639])
mean(ent.prob.vs.stability_v12v56.v12$identity.vector[640:787])