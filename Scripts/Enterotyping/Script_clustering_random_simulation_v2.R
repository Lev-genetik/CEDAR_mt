#check if the enterotyping with v1v2 v3v4 and v5v6 is not noise: do samples tend to cluster similarly?
#ent.Raes12.3e etc are results of /home/lev-genetik/Desktop/Projects/liege/src/Lev/Enterotyping_Raes_noise_removed_v6.R
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Enterotyping_compar_functions_v3.R')
freq_12 <- getEntFreq(ent.Raes12.2e[[2]],as.vector(ent.Raes12.2e[[1]][,1]))
freq_34 <- getEntFreq(ent.Raes34.2e[[2]],as.vector(ent.Raes34.2e[[1]][,1]))
freq_56 <- getEntFreq(ent.Raes56.2e[[2]],as.vector(ent.Raes56.2e[[1]][,1]))

source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/clustering_random_simulation_v2.R')
simulation_V12_V34 <- simulatn(freq_12,freq_34,permut.n = 10000)
simulation_V12_V56 <- simulatn(freq_12,freq_56,permut.n = 10000)
simulation_V34_V56 <- simulatn(freq_34,freq_56,permut.n = 10000)
