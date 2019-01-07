# This program draws a scatter plot with points representing enterotypes. Set1 - stable across amplicons
# set 2 - stable across couple of amplicons: 12-34, 34-56, 12-56
# set 3 - not stable across amplicons at all
plotEntStabil <- function(file1 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/V1V2_3e_IL_sample-enterotype.csv', 
                          file2 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/V3V4_3e_IL_sample-enterotype.csv', 
                          file3 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/V5V6_3e_IL_sample-enterotype.csv',
                          qiime_out = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt',
                          outdir='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs/BCA',
                          location='.IL',
                          podrovniat_files = FALSE, 
                          outliers.removed = FALSE) {
  
library(ade4)
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/GetStablyEnterotypedSamples_across_amplicons_smart_v3.R')
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Enterotyping_Raes_noise_removed_v10.R')
if (outliers.removed){podrovniat_files=T}
  #for naming the resultant plots, we need to preserve file names
  file1.name <- file1
  file2.name <- file2
  file3.name <- file3
  
  
No_Bacteroides_in_correspondence = F
No_Prevotella_in_correspondence = F
  
# this used in case we have files with non-100%-identical samples
if(podrovniat_files){
  df1 <- read.table(file1)
  df2 <- read.table(file2)
  df3 <- read.table(file3)
  # delete .V12 .V34 .V56
  df1[,1] <- gsub('\\.V[0-9][0-9]','',df1[,1])
  df2[,1] <- gsub('\\.V[0-9][0-9]','',df2[,1])
  df3[,1] <- gsub('\\.V[0-9][0-9]','',df3[,1])
  common.samples <- intersect(intersect(df1[,1],df2[,1]),df3[,1])
  df1 <- df1[df1[,1]%in%common.samples,]
  df2 <- df2[df2[,1]%in%common.samples,]
  df3 <- df3[df3[,1]%in%common.samples,]
  file1 = df1
  file2 = df2
  file3 = df3
}
# do the main job of enterotype classification. We get 5 data frames
# if podrovniat_files==TRUE than we inut here the data frames
classify.locat <- getStabAcrAmpl(file1,file2,file3,files=!podrovniat_files)

#we still need to do the job of enterotyping for the qiime output: as we need bca file 
ent.Raes.3e.curr.locat <- Enterotyping_Raes(qiime_out, '', clust.number=3, noiserem=T, locat=location, plots=F)
sample.enterotype <- ent.Raes.3e.curr.locat$sample.enterotype
sample.enterotype[,1] <- gsub('\\.V[0-9][0-9]','',sample.enterotype[,1])

# replace enterotype names by B, P and 3
if (length(levels(sample.enterotype[,2])[levels(sample.enterotype[,2])%like%'g__Bacteroides'])!=0){
  # usually it is this, simplest case
  levels(sample.enterotype[,2])[levels(sample.enterotype[,2])%like%'g__Bacteroides'] <- 'B'
} else{ #shit no Bacteroides, we need to get the enterotype corresponding to Bacteroides and name is e1
  if (!is.null(classify.locat$driver.correspondence[classify.locat$driver.correspondence[,3]%like%'g__Bacteroides',])){
  ent_instead_Bact.set <- classify.locat$driver.correspondence[classify.locat$driver.correspondence[,3]%like%'g__Bacteroides',]
  } else if (!is.null(classify.locat$driver.correspondence[classify.locat$driver.correspondence[,2]%like%'g__Bacteroides',])){
  ent_instead_Bact.set <- classify.locat$driver.correspondence[classify.locat$driver.correspondence[,2]%like%'g__Bacteroides',]
  } else if (!is.null(classify.locat$driver.correspondence[classify.locat$driver.correspondence[,1]%like%'g__Bacteroides',])){
    ent_instead_Bact.set <- classify.locat$driver.correspondence[classify.locat$driver.correspondence[,1]%like%'g__Bacteroides',]
  } else {
    # and this is real shit - we have no Bacteroides at all
    No_Bacteroides_in_correspondence = TRUE
  }
  if(!No_Bacteroides_in_correspondence){
  ent_instead_Bact <- ent_instead_Bact.set[ent_instead_Bact.set %in% levels(sample.enterotype[,2])]
  levels(sample.enterotype[,2])[levels(sample.enterotype[,2])==ent_instead_Bact] <- 'e1'
  }
}
if (length(levels(sample.enterotype[,2])[levels(sample.enterotype[,2])%like%'g__Prevotella'])!=0){
  # usually it is this, simplest case
levels(sample.enterotype[,2])[levels(sample.enterotype[,2])%like%'g__Prevotella'] <- 'P'
} else{ #shit
  if (!is.null(classify.locat$driver.correspondence[classify.locat$driver.correspondence[,3]%like%'g__Prevotella',])){
  ent_instead_Prev.set <- classify.locat$driver.correspondence[classify.locat$driver.correspondence[,3]%like%'g__Prevotella',]
  } else if (!is.null(classify.locat$driver.correspondence[classify.locat$driver.correspondence[,2]%like%'g__Prevotella',])){
    ent_instead_Prev.set <- classify.locat$driver.correspondence[classify.locat$driver.correspondence[,2]%like%'g__Prevotella',]
  } else if (!is.null(classify.locat$driver.correspondence[classify.locat$driver.correspondence[,1]%like%'g__Prevotella',])){
    ent_instead_Prev.set <- classify.locat$driver.correspondence[classify.locat$driver.correspondence[,1]%like%'g__Prevotella',]
  } else{
    No_Prevotella_in_correspondence = TRUE
  }
  if(!No_Prevotella_in_correspondence){
  ent_instead_Prev <- ent_instead_Prev.set[ent_instead_Prev.set %in% levels(sample.enterotype[,2])]
  levels(sample.enterotype[,2])[levels(sample.enterotype[,2])==ent_instead_Prev] <- 'e2'
  }
}
#  This if - very rare, if no Bacteroides/Prevotella in correspondence
if ((No_Prevotella_in_correspondence+No_Bacteroides_in_correspondence)==2){
  levels(sample.enterotype[,2])[1] <- 'e1'
  levels(sample.enterotype[,2])[2] <- 'e2'
} else if ((No_Prevotella_in_correspondence+No_Bacteroides_in_correspondence)==1){
  if (No_Prevotella_in_correspondence == 1){
    levels(sample.enterotype[,2])[!levels(sample.enterotype[,2])%in% c('e1','B')][1] <- 'e2'
  } else if (No_Bacteroides_in_correspondence == 1){
    levels(sample.enterotype[,2])[!levels(sample.enterotype[,2])%in% c('e2','P')][1] <- 'e1'
  }
}
# and finally e3
levels(sample.enterotype[,2])[!levels(sample.enterotype[,2])%in%c('P','B','e1','e2')] <- 'e3'

# rows of samples that are stable for all amplicons in IL
stable.locat.12_34_56_numbers <- sample.enterotype[,1] %in% classify.locat$result.v12v34v56[,1]
stable.locat.12_34_numbers <- sample.enterotype[,1] %in% classify.locat$result.v12v34[,1]
stable.locat.12_56_numbers <- sample.enterotype[,1] %in% classify.locat$result.v12v56[,1]
stable.locat.34_56_numbers <- sample.enterotype[,1] %in% classify.locat$result.v34v56[,1]
unstable.locat_numbers <- sample.enterotype[,1] %in% classify.locat$result.nonconserved[,1]

# drawing scatter plot
# first, compose the name of the file: amplicon taken froma qiime_out, number of ent, location - from file1.
vector_tmp <- word(word(file1.name,start = -1,sep = '/'),start=c(1,2,3),sep = '_')
vector_tmp[1] <- word(qiime_out,start = -2,sep = '/')
vector_tmp1 <- paste(vector_tmp[1], vector_tmp[2], vector_tmp[3], sep ='_')
if(outliers.removed){vector_tmp1 <- paste0(vector_tmp1,'_outliers_remov')}
file_out1 <- paste0(vector_tmp1,'_bca_first_2_eignenvectors_thin.jpg')

setwd(outdir)
jpeg(filename = file_out1,width = 5000,height = 4000)
legend <- 'Red - stable all ampl., green (+,X,*) - stable v12-v34, v12-v56, v34-56 resp., blue - unstable'
xlim = c(min(ent.Raes.3e.curr.locat$obs_bet_result$ls[,1])-0.5,max(ent.Raes.3e.curr.locat$obs_bet_result$ls[,1])+0.5)
ylim = c(min(ent.Raes.3e.curr.locat$obs_bet_result$ls[,2])-1,max(ent.Raes.3e.curr.locat$obs_bet_result$ls[,2])+1)

# samples stable across 3 amplicons = red
s.class(ent.Raes.3e.curr.locat$obs_bet_result$ls[stable.locat.12_34_56_numbers,], 
        fac=as.factor(sample.enterotype$enterotype[stable.locat.12_34_56_numbers]), 
        col=c('red3','red3','red3','red3'), 
        xlim=xlim, 
        ylim=ylim, 
        cpoint = 4, 
        grid=F, 
        sub=legend, 
        cstar = 1, 
        csub = 10, 
        cellipse = 1, 
        pch = 12)
# samples stable v12-v34
s.class(ent.Raes.3e.curr.locat$obs_bet_result$ls[stable.locat.12_34_numbers,], fac=as.factor(sample.enterotype$enterotype[stable.locat.12_34_numbers]), col=c('springgreen4','springgreen4','springgreen4'), cpoint = 4, grid=F, cellipse = 0, cstar = 0, label = NULL, pch = 3, add.plot = T)
# samples stable v12-v56
s.class(ent.Raes.3e.curr.locat$obs_bet_result$ls[stable.locat.12_56_numbers,], fac=as.factor(sample.enterotype$enterotype[stable.locat.12_56_numbers]), col=c('springgreen4','springgreen4','springgreen4'), cpoint = 4, grid=F, cellipse = 0, cstar = 0, label = NULL, pch = 4, add.plot = T)
# samples stable v34-v56
s.class(ent.Raes.3e.curr.locat$obs_bet_result$ls[stable.locat.34_56_numbers,], fac=as.factor(sample.enterotype$enterotype[stable.locat.34_56_numbers]), col=c('springgreen4','springgreen4','springgreen4'), cpoint = 4, grid=F, cellipse = 0, cstar = 0, label = NULL, pch = 8, add.plot = T)
# unstable samples = blue
s.class(ent.Raes.3e.curr.locat$obs_bet_result$ls[unstable.locat_numbers,], fac=as.factor(sample.enterotype$enterotype[unstable.locat_numbers]), col=c('blue','blue','blue'), cpoint = 4, grid=F, cellipse = 0, cstar = 0, label = NULL, pch = 13, add.plot = T)
dev.off()

file_out2 <- paste0(vector_tmp1,'_bca_first_2_eignenvectors_thick.jpg')
# and this old style - good for seing in low resolution
jpeg(filename = file_out2,width = 2000,height = 1600)
legend <- 'Red - stable all ampl., yellow - stable v12-v34, orange v12-v56, green v34-56, blue - unstable'
xlim = c(min(ent.Raes.3e.curr.locat$obs_bet_result$ls[,1])-0.5,max(ent.Raes.3e.curr.locat$obs_bet_result$ls[,1])+0.5)
ylim = c(min(ent.Raes.3e.curr.locat$obs_bet_result$ls[,2])-1,max(ent.Raes.3e.curr.locat$obs_bet_result$ls[,2])+1)
lwd = 10
s.class(ent.Raes.3e.curr.locat$obs_bet_result$ls[stable.locat.12_34_56_numbers,], fac=as.factor(sample.enterotype$enterotype[stable.locat.12_34_56_numbers]), col=c('red','red','red','red'), xlim=xlim, ylim=ylim, cpoint = 2, grid=F, sub=legend, cstar = 1, csub = 4, cellipse = 1)
s.class(ent.Raes.3e.curr.locat$obs_bet_result$ls[stable.locat.12_34_numbers,], fac=as.factor(sample.enterotype$enterotype[stable.locat.12_34_numbers]), col=c('orange','orange','orange','orange'), cpoint = 2, grid=F, cellipse = 0, cstar = 0, label = NULL, add.plot = T)
s.class(ent.Raes.3e.curr.locat$obs_bet_result$ls[stable.locat.12_56_numbers,], fac=as.factor(sample.enterotype$enterotype[stable.locat.12_56_numbers]), col=c('yellow','yellow','yellow','yellow'), cpoint = 2, grid=F, cellipse = 0, cstar = 0, label = NULL, add.plot = T)
s.class(ent.Raes.3e.curr.locat$obs_bet_result$ls[stable.locat.34_56_numbers,], fac=as.factor(sample.enterotype$enterotype[stable.locat.34_56_numbers]), col=c('green','green','green','green'), cpoint = 2, grid=F, cellipse = 0, cstar = 0, label = NULL, add.plot = T)
s.class(ent.Raes.3e.curr.locat$obs_bet_result$ls[unstable.locat_numbers,], fac=as.factor(sample.enterotype$enterotype[unstable.locat_numbers]), col=c('blue','blue','blue','blue'), cpoint = 2, grid=F, cellipse = 0, cstar = 0, label = NULL, add.plot = T)
dev.off()

# s.label(ent.Raes.3e.curr.locat$obs_bet_result$ls[stable.locat.12_34_56_numbers,],sub=legend, csub = 1.5, possub = "bottmoright")

}