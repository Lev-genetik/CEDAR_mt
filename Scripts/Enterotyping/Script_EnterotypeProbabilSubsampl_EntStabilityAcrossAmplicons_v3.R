# Script to get stable enterotypes for each amplicon by subsampling (= probatypes). Get enterotypes for half of the samples
# and asssign stable enterotypes only to samples that have the same enterotype in >=90% of cases
# all locations enterotyped together
# for v1v2, v3v4 and v5v6
# Than, draw graphs of the enterotype probability density
# Than, for each amplicon, sort samples according to enterotype probability in descending order and check
# if the probability of having stable enterotype across amplicons also decreases.

setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities')
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Enterotype_Probabilities_functions_v2.R')
# v1v2 probability enterotyping
sample.enterotype.probabil.subsampl.v12 <- getStablEntSubsampl(input_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt',enterot_file='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype/V1V2_3e_sample-enterotype.csv',N=1000)
write.table(sample.enterotype.probabil.subsampl.v12,'V1V2_3e_sample-enterotype_probabilities.csv')
# v3v4 probability enterotyping
sample.enterotype.probabil.subsampl.v34 <- getStablEntSubsampl(input_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt',enterot_file='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype/V3V4_3e_sample-enterotype.csv',N=1000)
write.table(sample.enterotype.probabil.subsampl.v34,'V3V4_3e_sample-enterotype_probabilities.csv')
# v5v6 probability enterotyping
sample.enterotype.probabil.subsampl.v56 <- getStablEntSubsampl(input_file = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt',enterot_file='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype/V5V6_3e_sample-enterotype.csv',N=1000)
write.table(sample.enterotype.probabil.subsampl.v56,'V5V6_3e_sample-enterotype_probabilities.csv')

# now draw graphs of enterotype probability density
library(graphics)
setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities')
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
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Enterotype_Probabilities_functions_v2.R')
setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities/Graphs_enterotype_probability')

# FOR V12 vs v12-v34-v56
# get number of samples with permatype=enterotype
perm.v12 <- read.table('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities/V1V2_3e_sample-enterotype_probabilities.csv')
ent.perm.eq.v12 <- sum(perm.v12$Enterotype_Probatype_consist)
print(paste0('Total number of samples where enterotype = permatype: ',ent.perm.eq.v12))
# get the result (stability.window, identity.vector, ent.stability.besttoworst.sum.first.n.samples) and draw 
# graph of enterotype stability depending on enterotype probability in descend order
ent.prob.vs.stability.v12 <- subsVsLoc(file1 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities/V1V2_3e_sample-enterotype_probabilities.csv',file2='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype/Stable_enterotypes.csv', folder = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities/Graphs_enterotype_probability',graph.file1 = 'v1v2_ent_probability_by_sampl_vs_v12v34v56.jpg',graph.file2='Enterotype_probabil_v12_descend_cumulat_vs_ent_stability_v12_v34_v56.jpg',window=200)
# Now according to graph of v1v2_ent_probability_by_sampl_vs_v12v34v56.jpg and ent.prob.vs.stability.v12$stability.window we see
ent.prob.vs.stability.v12$stability.window
# that the huge quality descend starts after 484+200-1=683 that is the last good sample.
# Get the Enterotype (probatype) Probability threshold for rank_cutoff=683
getProbabilCutoff(file='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities/V1V2_3e_sample-enterotype_probabilities.csv',rank_cutoff=683)
# here: 0.865
# Now among the samples with Probatype=Enterotype, see what is the difference between Enterotype stability before and after 
# enterotype probability rank_cutoff
last.samurai <- 683
before <- signif(mean(ent.prob.vs.stability.v12$identity.vector[1:last.samurai]),digits = 3)
print(paste0('For V12 vs v12-v34-v56 the enterotype stability for samples of good-quality-probatypes: ',before))
after <- signif(mean(ent.prob.vs.stability.v12$identity.vector[(last.samurai+1):length(ent.prob.vs.stability.v12$identity.vector)]),digits = 3)
print(paste0('For V12 vs v12-v34-v56 the enterotype stability for samples of bad-quality-probatypes: ',after))
last.samurai <-ent.perm.eq.v12
before <- signif(mean(ent.prob.vs.stability.v12$identity.vector[1:last.samurai]),digits = 3)
print(paste0('For v12 vs v12-v34-v56 the enterotype stability for all samples where enterotype=probatype: ',before))
# Now get maximal enterotype stability cutoff
number.100.stab <- 512  # this is number of first samples with probatype stability 100%
array <- ent.prob.vs.stability.v12$ent.stability.besttoworst.sum.first.n.samples[number.100.stab:length(ent.prob.vs.stability.v12$ent.stability.besttoworst.sum.first.n.samples)]
dobavka <- max(which(array==max(array)))
last.samurai <- number.100.stab + dobavka - 1
max(ent.prob.vs.stability.v12$ent.stability.besttoworst.sum.first.n.samples[number.100.stab:length(ent.prob.vs.stability.v12$ent.stability.besttoworst.sum.first.n.samples)])
getProbabilCutoff(file='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities/V1V2_3e_sample-enterotype_probabilities.csv',rank_cutoff=last.samurai)
print(paste0('last of the first best samples is: ',last.samurai))
before <- signif(mean(ent.prob.vs.stability.v12$identity.vector[1:last.samurai]),digits = 3)
print(paste0('For v12 vs v12-v34-v56 the enterotype stability for samples of best-quality-probatypes: ',before))
after <- signif(mean(ent.prob.vs.stability.v12$identity.vector[(last.samurai+1):length(ent.prob.vs.stability.v12$identity.vector)]),digits = 3)
print(paste0('For v12 vs v12-v34-v56 the enterotype stability for samples of bad-quality-probatypes: ',after))

# The same for V12 vs v12-v56
ent.prob.vs.stability_v12v56.v12 <- subsVsLoc(file1 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities/V1V2_3e_sample-enterotype_probabilities.csv',file2='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype/Stable_enterotypes_v12_v56.csv', folder = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities/Graphs_enterotype_probability',graph.file1 = 'v1v2_ent_probability_by_sampl_vs_v12v56.jpg',graph.file2='Enterotype_probabil_v12_descend_cumulat_vs_ent_stability_v12_v56.jpg',window=200)
# to determine last peak accuartely
ent.prob.vs.stability_v12v56.v12$stability.window
# here, cutoff rank = 494+200-1=693
last.samurai <- 693
getProbabilCutoff(file='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities/V1V2_3e_sample-enterotype_probabilities.csv',rank_cutoff=last.samurai)
# here: 0.82
before <- signif(mean(ent.prob.vs.stability_v12v56.v12$identity.vector[1:last.samurai]),digits = 3)
print(paste0('For V12 vs v12-v56 the enterotype stability for samples of good-quality-probatypes: ',before))
after <- signif(mean(ent.prob.vs.stability_v12v56.v12$identity.vector[(last.samurai+1):length(ent.prob.vs.stability_v12v56.v12$identity.vector)]),digits = 3)
print(paste0('For V12 vs v12-v56 the enterotype stability for samples of bad-quality-probatypes: ',after))
# Now get maximal enterotype stability cutoff
number.100.stab <- 512  # this is number of first samples with probatype stability 100%
# which.max(ent.prob.vs.stability_v12v56.v12$ent.stability.besttoworst.sum.first.n.samples[number.100.stab:length(ent.prob.vs.stability_v12v56.v12$ent.stability.besttoworst.sum.first.n.samples)])
# last.samurai <- number.100.stab + which.max(ent.prob.vs.stability_v12v56.v12$ent.stability.besttoworst.sum.first.n.samples[number.100.stab:length(ent.prob.vs.stability_v12v56.v12$ent.stability.besttoworst.sum.first.n.samples)])
# max(ent.prob.vs.stability_v12v56.v12$ent.stability.besttoworst.sum.first.n.samples[number.100.stab:length(ent.prob.vs.stability_v12v56.v12$ent.stability.besttoworst.sum.first.n.samples)])
# getProbabilCutoff(file='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities/V1V2_3e_sample-enterotype_probabilities.csv',rank_cutoff=last.samurai)
array <- ent.prob.vs.stability_v12v56.v12$ent.stability.besttoworst.sum.first.n.samples[number.100.stab:length(ent.prob.vs.stability_v12v56.v12$ent.stability.besttoworst.sum.first.n.samples)]
dobavka <- max(which(array==max(array)))
last.samurai <- number.100.stab + dobavka - 1
max(ent.prob.vs.stability_v12v56.v12$ent.stability.besttoworst.sum.first.n.samples[number.100.stab:length(ent.prob.vs.stability_v12v56.v12$ent.stability.besttoworst.sum.first.n.samples)])
getProbabilCutoff(file='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities/V1V2_3e_sample-enterotype_probabilities.csv',rank_cutoff=last.samurai)
print(paste0('last of the first best samples is: ',last.samurai))
before <- signif(mean(ent.prob.vs.stability_v12v56.v12$identity.vector[1:last.samurai]),digits = 3)
print(paste0('For v12 vs v12-v56 the enterotype stability for samples of best-quality-probatypes: ',before))
after <- signif(mean(ent.prob.vs.stability_v12v56.v12$identity.vector[(last.samurai+1):length(ent.prob.vs.stability_v12v56.v12$identity.vector)]),digits = 3)
print(paste0('For v12 vs v12-v56 the enterotype stability for samples of bad-quality-probatypes: ',after))



# For V34 vs v12-v34-v56
# get number of samples with permatype=enterotype
perm.v34 <- read.table('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities/V3V4_3e_sample-enterotype_probabilities.csv')
ent.perm.eq.v34 <- sum(perm.v34$Enterotype_Probatype_consist)
print(paste0('Total number of samples where enterotype = permatype: ',ent.perm.eq.v34))
# get the result (stability.window, identity.vector, ent.stability.besttoworst.sum.first.n.samples) and draw 
# graph of enterotype stability depending on enterotype probability in descend order
ent.prob.vs.stability.v34 <- subsVsLoc(file1 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities/V3V4_3e_sample-enterotype_probabilities.csv',file2='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype/Stable_enterotypes.csv', folder = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities/Graphs_enterotype_probability',graph.file1 = 'v3v4_ent_probability_by_sampl_vs_v12v34v56.jpg',graph.file2='Enterotype_probabil_v34_descend_cumulat_vs_ent_stability_v12_v34_v56.jpg',window=200)
# Now according to graph of v3v4_ent_probability_by_sampl_vs_v12v34v56.jpg and ent.prob.vs.stability.v34$stability.window we see
ent.prob.vs.stability.v34$stability.window
# that the huge quality descend starts after 240+200-1 = 439 [337+200-1=536 ]that is the last good sample.
# Get the Enterotype (probatype) Probability threshold for rank_cutoff=439
getProbabilCutoff(file='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities/V3V4_3e_sample-enterotype_probabilities.csv',rank_cutoff=439)
# here: 0.927
# Now among the samples with Probatype=Enterotype, see what is the difference between Enterotype stability before and after 
# enterotype probability rank_cutoff
last.samurai <- 439
before <- signif(mean(ent.prob.vs.stability.v34$identity.vector[1:last.samurai]),digits = 3)
print(paste0('For v34 vs v12-v34-v56 the enterotype stability for samples of good-quality-probatypes: ',before))
after <- signif(mean(ent.prob.vs.stability.v34$identity.vector[(last.samurai+1):length(ent.prob.vs.stability.v34$identity.vector)]),digits = 3)
print(paste0('For v34 vs v12-v34-v56 the enterotype stability for samples of bad-quality-probatypes: ',after))
last.samurai <-ent.perm.eq.v34
before <- signif(mean(ent.prob.vs.stability.v34$identity.vector[1:last.samurai]),digits = 3)
print(paste0('For v34 vs v12-v34-v56 the enterotype stability for all samples where enterotype=probatype: ',before))
# Now get maximal enterotype stability cutoff
number.100.stab <- 300  # this is number of first samples with probatype stability 100%
array <- ent.prob.vs.stability.v34$ent.stability.besttoworst.sum.first.n.samples[number.100.stab:length(ent.prob.vs.stability.v34$ent.stability.besttoworst.sum.first.n.samples)]
dobavka <- max(which(array==max(array)))
last.samurai <- number.100.stab + dobavka - 1
max(ent.prob.vs.stability.v34$ent.stability.besttoworst.sum.first.n.samples[number.100.stab:length(ent.prob.vs.stability.v34$ent.stability.besttoworst.sum.first.n.samples)])
getProbabilCutoff(file='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities/V3V4_3e_sample-enterotype_probabilities.csv',rank_cutoff=last.samurai)
print(paste0('last of the first best samples is: ',last.samurai))
before <- signif(mean(ent.prob.vs.stability.v34$identity.vector[1:last.samurai]),digits = 3)
print(paste0('For v34 vs v12-v34-v56 the enterotype stability for samples of best-quality-probatypes: ',before))
after <- signif(mean(ent.prob.vs.stability.v34$identity.vector[(last.samurai+1):length(ent.prob.vs.stability.v34$identity.vector)]),digits = 3)
print(paste0('For v34 vs v12-v34-v56 the enterotype stability for samples of bad-quality-probatypes: ',after))

# The same for v34 vs v12-v56 Did not do this as it has no much sense v34 vs v12-v56 stable..
ent.prob.vs.stability_v12v56.v34 <- subsVsLoc(file1 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities/v3v4_3e_sample-enterotype_probabilities.csv',file2='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype/Stable_enterotypes_v12_v56.csv', folder = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities/Graphs_enterotype_probability',graph.file1 = 'v3v4_ent_probability_by_sampl_vs_v12v56.jpg',graph.file2='Enterotype_probabil_v34_descend_cumulat_vs_ent_stability_v12_v56.jpg',window=200)
# to determine last peak accuartely
ent.prob.vs.stability_v12v56.v34$stability.window
# here, cutoff rank = 494+200-1=693
last.samurai <- 693
getProbabilCutoff(file='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities/v3v4_3e_sample-enterotype_probabilities.csv',rank_cutoff=last.samurai)
# here: 0.82
before <- signif(mean(ent.prob.vs.stability_v12v56.v34$identity.vector[1:last.samurai]),digits = 3)
print(paste0('For v34 vs v12-v56 the enterotype stability for samples of good-quality-probatypes: ',before))
after <- signif(mean(ent.prob.vs.stability_v12v56.v34$identity.vector[(last.samurai+1):length(ent.prob.vs.stability_v12v56.v34$identity.vector)]),digits = 3)
print(paste0('For v34 vs v12-v56 the enterotype stability for samples of bad-quality-probatypes: ',after))
# Now get maximal enterotype stability cutoff
number.100.stab <- 512  # this is number of first samples with probatype stability 100%
which.max(ent.prob.vs.stability_v12v56.v34$ent.stability.besttoworst.sum.first.n.samples[number.100.stab:length(ent.prob.vs.stability_v12v56.v34$ent.stability.besttoworst.sum.first.n.samples)])
last.samurai <- number.100.stab + which.max(ent.prob.vs.stability_v12v56.v34$ent.stability.besttoworst.sum.first.n.samples[number.100.stab:length(ent.prob.vs.stability_v12v56.v34$ent.stability.besttoworst.sum.first.n.samples)])
max(ent.prob.vs.stability_v12v56.v34$ent.stability.besttoworst.sum.first.n.samples[number.100.stab:length(ent.prob.vs.stability_v12v56.v34$ent.stability.besttoworst.sum.first.n.samples)])
getProbabilCutoff(file='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities/v3v4_3e_sample-enterotype_probabilities.csv',rank_cutoff=last.samurai)



# FOR V56 vs v12-v34-v56
# get number of samples with permatype=enterotype
perm.v56 <- read.table('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities/V5V6_3e_sample-enterotype_probabilities.csv')
ent.perm.eq.v56 <- sum(perm.v56$Enterotype_Probatype_consist)
print(paste0('Total number of samples where enterotype = permatype: ',ent.perm.eq.v56))
# get the result (stability.window, identity.vector, ent.stability.besttoworst.sum.first.n.samples) and draw 
# graph of enterotype stability depending on enterotype probability in descend order
ent.prob.vs.stability.v56 <- subsVsLoc(file1 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities/V5V6_3e_sample-enterotype_probabilities.csv',file2='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype/Stable_enterotypes.csv', folder = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities/Graphs_enterotype_probability',graph.file1 = 'V5V6_ent_probability_by_sampl_vs_v12v34v56.jpg',graph.file2='Enterotype_probabil_v56_descend_cumulat_vs_ent_stability_v12_v34_v56.jpg',window=200)
# Now according to graph of V5V6_ent_probability_by_sampl_vs_v12v34v56.jpg and ent.prob.vs.stability.v56$stability.window we see
ent.prob.vs.stability.v56$stability.window
# that the huge quality descend starts after 406+200-1=605 that is the last good sample.
# Get the Enterotype (probatype) Probability threshold for rank_cutoff=605
getProbabilCutoff(file='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities/V5V6_3e_sample-enterotype_probabilities.csv',rank_cutoff=605)
# here: 0.976
# Now among the samples with Probatype=Enterotype, see what is the difference between Enterotype stability before and after 
# enterotype probability rank_cutoff
last.samurai <- 605
before <- signif(mean(ent.prob.vs.stability.v56$identity.vector[1:last.samurai]),digits = 3)
print(paste0('For v56 vs v12-v34-v56 the enterotype stability for samples of good-quality-probatypes: ',before))
after <- signif(mean(ent.prob.vs.stability.v56$identity.vector[(last.samurai+1):length(ent.prob.vs.stability.v56$identity.vector)]),digits = 3)
print(paste0('For v56 vs v12-v34-v56 the enterotype stability for samples of bad-quality-probatypes: ',after))
last.samurai <-ent.perm.eq.v56
before <- signif(mean(ent.prob.vs.stability.v56$identity.vector[1:last.samurai]),digits = 3)
print(paste0('For v56 vs v12-v34-v56 the enterotype stability for all samples where enterotype=probatype: ',before))
# Now get maximal enterotype stability cutoff
number.100.stab <- 481  # this is number of first samples with probatype stability 100%
array <- ent.prob.vs.stability.v56$ent.stability.besttoworst.sum.first.n.samples[number.100.stab:length(ent.prob.vs.stability.v56$ent.stability.besttoworst.sum.first.n.samples)]
dobavka <- max(which(array==max(array)))
last.samurai <- number.100.stab + dobavka - 1
max(ent.prob.vs.stability.v56$ent.stability.besttoworst.sum.first.n.samples[number.100.stab:length(ent.prob.vs.stability.v56$ent.stability.besttoworst.sum.first.n.samples)])
print(paste0('last of the first best samples is: ',last.samurai))
getProbabilCutoff(file='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities/V5V6_3e_sample-enterotype_probabilities.csv',rank_cutoff=last.samurai)
before <- signif(mean(ent.prob.vs.stability.v56$identity.vector[1:last.samurai]),digits = 3)
print(paste0('For v56 vs v12-v34-v56 the enterotype stability for samples of best-quality-probatypes: ',before))
after <- signif(mean(ent.prob.vs.stability.v56$identity.vector[(last.samurai+1):length(ent.prob.vs.stability.v56$identity.vector)]),digits = 3)
print(paste0('For v56 vs v12-v34-v56 the enterotype stability for samples of bad-quality-probatypes: ',after))

# The same for v56 vs v12-v56
ent.prob.vs.stability_v12v56.v56 <- subsVsLoc(file1 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities/V5V6_3e_sample-enterotype_probabilities.csv',file2='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype/Stable_enterotypes_v12_v56.csv', folder = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities/Graphs_enterotype_probability',graph.file1 = 'V5V6_ent_probability_by_sampl_vs_v12v56.jpg',graph.file2='Enterotype_probabil_v56_descend_cumulat_vs_ent_stability_v12_v56.jpg',window=200)
# to determine last peak accuartely
ent.prob.vs.stability_v12v56.v56$stability.window
# here, cutoff rank = 418+200-1=617
last.samurai <- 617
getProbabilCutoff(file='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities/V5V6_3e_sample-enterotype_probabilities.csv',rank_cutoff=last.samurai)
# here: 0.966
before <- signif(mean(ent.prob.vs.stability_v12v56.v56$identity.vector[1:last.samurai]),digits = 3)
print(paste0('For v56 vs v12-v56 the enterotype stability for samples of good-quality-probatypes: ',before))
after <- signif(mean(ent.prob.vs.stability_v12v56.v56$identity.vector[(last.samurai+1):length(ent.prob.vs.stability_v12v56.v56$identity.vector)]),digits = 3)
print(paste0('For v56 vs v12-v56 the enterotype stability for samples of bad-quality-probatypes: ',after))
# Now get maximal enterotype stability cutoff
number.100.stab <- 481  # this is number of first samples with probatype stability 100%
array <- ent.prob.vs.stability_v12v56.v56$ent.stability.besttoworst.sum.first.n.samples[number.100.stab:length(ent.prob.vs.stability_v12v56.v56$ent.stability.besttoworst.sum.first.n.samples)]
dobavka <- max(which(array==max(array)))
last.samurai <- number.100.stab + dobavka - 1
max(ent.prob.vs.stability_v12v56.v56$ent.stability.besttoworst.sum.first.n.samples[number.100.stab:length(ent.prob.vs.stability_v12v56.v56$ent.stability.besttoworst.sum.first.n.samples)])
print(paste0('last of the first best samples is: ',last.samurai))
getProbabilCutoff(file='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Enterotype_Probabilities/V5V6_3e_sample-enterotype_probabilities.csv',rank_cutoff=last.samurai)
before <- signif(mean(ent.prob.vs.stability_v12v56.v56$identity.vector[1:last.samurai]),digits = 3)
print(paste0('For v56 vs v12-v56 the enterotype stability for samples of best-quality-probatypes: ',before))
after <- signif(mean(ent.prob.vs.stability_v12v56.v56$identity.vector[(last.samurai+1):length(ent.prob.vs.stability_v12v56.v56$identity.vector)]),digits = 3)
print(paste0('For v56 vs v12-v56 the enterotype stability for samples of bad-quality-probatypes: ',after))


last.samurai <-ent.perm.eq.v56
before <- signif(mean(ent.prob.vs.stability_v12v56.v56$identity.vector[1:last.samurai]),digits = 3)
print(paste0('For v56 vs v12-v34-v56 the enterotype stability for all samples where enterotype=probatype: ',before))

