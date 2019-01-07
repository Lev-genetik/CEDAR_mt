

# see how different are locations in terms of PCoAs and PCA (PCs)
# UNPUT:
# PCoA.pair=T: do draw the plot of the 2 PCAs from draw.PCoA.pairwise for different locations
# PCA.pair=T: do draw the plot of the 2 PCAs from draw.PCA.pairwise for different locations
# PCoA.individual=T: do draw the plots of individual PCoAs for different locations
# PCA.individual=T: do draw the plots of individual PCs for different locations
# draw.PCoA.sep=seq(1,9): make plots for location difference across PCoAs 1 to 9
# draw.PCA.sep=seq(1,9): make plots for location difference across PCAs 1 to 9
# draw.PCoA.pairwise=c(1,2): make plots for location difference in the coordinates of PCoAs 1 and 2
# draw.PCA.pairwise = c(1,2):  make plots for location difference in the coordinates of PCAs 1 and 2
# sample coverage filter and noise removal can be specified
# qiime_out = qiiime table used as start data
# outdir = folder used for drawing plots
# locations=c('.IL','.TR','.RE')
# outliers.removed = TRUE (remove the 1% most distant samples while enterotyping)
# sample.coverage.filter = determine sample coverage filter. e.g. v1v2 23000, v3v4 13000, v5v6 21000
plotIntestineLocations <- function(qiime_out = file.path(home.dir,'Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt'),
                                   outdir=file.path(home.dir,'Results/Locations_compare/Graphs'),
                                   locations=c('.IL','.TR','.RE'), outliers.removed = TRUE, sample.coverage.filter=10000,
                                   PCoA.pair=T,PCA.pair=T,PCoA.individual=T,PCA.individual=T,
                                   draw.PCoA.sep=seq(1,9), draw.PCA.sep=seq(1,9), draw.PCoA.pairwise = c(1,2),
                                   draw.PCA.pairwise=c(1,2)){
  if(length(draw.PCA.sep)>9){
    stop('length of draw.PCA.sep must be <=9')
  }
  if(length(draw.PCoA.sep)>9){
    stop('length of draw.PCoA.sep must be <=9')
  }
  library(ade4)
  library(stringr)
  source(file.path(progr.dir,'GetStablyEnterotypedSamples_across_amplicons_smart_v5.R'))
  source(file.path(progr.dir,'Enterotyping_Raes_noise_removed_v12.R'))
  tmp.dir <- file.path(home.dir,'tmp')
  dir.create(tmp.dir,showWarnings = FALSE)
  ent.Raes.2e.curr.locat <- Enterotyping_Raes(qiime_out, output_dir = tmp.dir, clust.number=2, noiserem=T, plots=T,
                                              sample.coverage.filter = sample.coverage.filter,
                                              remove.outliers = 0.01*outliers.removed)
  sample.enterotype <- ent.Raes.2e.curr.locat$sample.enterotype
  # sample.enterotype[,1] <- gsub('\\.V[0-9][0-9]','',sample.enterotype[,1])
  
  
  # rows of samples that are stable for all amplicons in IL
  IL.samples <- sample.enterotype[,1] %like% '\\.IL'
  TR.samples <- sample.enterotype[,1] %like% '\\.TR'
  RE.samples <- sample.enterotype[,1] %like% '\\.RE'
  samples.location <- as.character(sample.enterotype[,1])
  samples.location[samples.location %like% '\\.IL'] <- 'ileum'
  samples.location[samples.location %like% '\\.TR'] <- 'transverse colon'
  samples.location[samples.location %like% '\\.RE'] <- 'rectum'

  # first, compose the name of the files
  file_name_start <- word(qiime_out,start = -2,sep = '/')
  file_name_start <- paste0(file_name_start,'_cov_thresh_',round(sample.coverage.filter/1000),'k')
  if(outliers.removed){file_name_start <- paste0(file_name_start,'_outliers_remov')}
  file_out.PCoA.pair <- paste0(file_name_start,'_PCoA_eignenvectors_',draw.PCoA.pairwise[1],'_',draw.PCoA.pairwise[2],'.jpg')
  file_out.PCA.pair <- paste0(file_name_start,'_PCA_eignenvectors_',draw.PCA.pairwise[1],'_',draw.PCA.pairwise[2],'.jpg')
  file_out.PCoA.individual <- paste0(file_name_start,'_PCoAs',str_c(draw.PCoA.sep,collapse = '_'),'.jpg')
  file_out.PCA.individual <- paste0(file_name_start,'_PCAs',str_c(draw.PCA.sep,collapse = '_'),'.jpg')
  # set wd
  setwd(outdir)

  # now, draw the plots
  # PCoA pair
  if(PCoA.pair){
    jpeg(filename = file_out.PCoA.pair,width = 2000,height = 1600)
    legend <- 'red - ileum., green - transverse colon, purple - rectum'
    xlim = c(1.01*min(ent.Raes.2e.curr.locat$obs_pcoa_result$li[,draw.PCoA.pairwise[1]]),1.01*max(ent.Raes.2e.curr.locat$obs_pcoa_result$li[,draw.PCoA.pairwise[1]]))
    ylim = c(1.01*min(ent.Raes.2e.curr.locat$obs_pcoa_result$li[,draw.PCoA.pairwise[2]]),1.01*max(ent.Raes.2e.curr.locat$obs_pcoa_result$li[,draw.PCoA.pairwise[2]]))
    var.1 <- paste0(100*signif(ent.Raes.2e.curr.locat$obs_pcoa_result$eig[draw.PCoA.pairwise[1]]/sum(ent.Raes.2e.curr.locat$obs_pcoa_result$eig),digits = 3),'%')
    var.2 <- paste0(100*signif(ent.Raes.2e.curr.locat$obs_pcoa_result$eig[draw.PCoA.pairwise[2]]/sum(ent.Raes.2e.curr.locat$obs_pcoa_result$eig),digits = 3),'%')
    subtitle.PCoA <- paste0('PC ',draw.PCoA.pairwise[1],' (',var.1,'), PC ',draw.PCoA.pairwise[2], ' (',var.2,')')
    s.class(ent.Raes.2e.curr.locat$obs_pcoa_result$li[,draw.PCoA.pairwise], 
            fac=as.factor(samples.location), #wierd 
            col=c('red','purple','darkgreen'), 
            xlim=xlim, 
            ylim=ylim, 
            cpoint = 4, 
            clabel=4,
            grid=F, 
            sub=paste0(subtitle.PCoA,'; ',legend), 
            cstar = 0, 
            csub = 4, 
            cellipse = 1,
            pch = 20)
    dev.off()
  }

  # PCA.pair
  if(PCA.pair){
  jpeg(filename = file_out.PCA.pair,width = 2000,height = 1600)
  legend <- 'red - ileum., green - transverse colon, purple - rectum'
  xlim = c(1.01*min(ent.Raes.2e.curr.locat$obs_pca_result$li[,draw.PCA.pairwise[1]]),1.01*max(ent.Raes.2e.curr.locat$obs_pca_result$li[,draw.PCA.pairwise[1]]))
  ylim = c(1.01*min(ent.Raes.2e.curr.locat$obs_pca_result$li[,draw.PCA.pairwise[2]]),1.01*max(ent.Raes.2e.curr.locat$obs_pca_result$li[,draw.PCA.pairwise[2]]))
  var.1 <- paste0(100*signif(ent.Raes.2e.curr.locat$obs_pca_result$eig[draw.PCA.pairwise[1]]/sum(ent.Raes.2e.curr.locat$obs_pca_result$eig),digits = 3),'%')
  var.2 <- paste0(100*signif(ent.Raes.2e.curr.locat$obs_pca_result$eig[draw.PCA.pairwise[2]]/sum(ent.Raes.2e.curr.locat$obs_pca_result$eig),digits = 3),'%')
  subtitle.PCA <- paste0('PC ',draw.PCA.pairwise[1],' (',var.1,'), PC ',draw.PCA.pairwise[2], ' (',var.2,')')
  s.class(ent.Raes.2e.curr.locat$obs_pca_result$li[,draw.PCA.pairwise], 
          fac=as.factor(samples.location), #wierd 
          col=c('red','purple','darkgreen'), 
          xlim=xlim, 
          ylim=ylim,
          cpoint = 2, 
          grid=F, 
          sub=paste0(subtitle.PCA,'; ',legend), 
          cstar = 0, 
          csub = 4, 
          cellipse = 1,
          pch = 20)
  dev.off()
  }
  
  # PCoA individual PCs
  if(PCoA.individual){
    # draw samples' density vs PCos
    jpeg(file_out.PCoA.individual,width=4000,height=3200)
    par(mfrow=c(3,3),adj=0.5,mar=c(15,17,10,2),mgp=c(11,5,0))
    for (i in draw.PCoA.sep){
      TR.IL <- (mean(ent.Raes.2e.curr.locat$obs_pcoa_result$li[TR.samples,i])-
                  mean(ent.Raes.2e.curr.locat$obs_pcoa_result$li[IL.samples,i]))/
        (sqrt(sd(ent.Raes.2e.curr.locat$obs_pcoa_result$li[IL.samples,i])*
                sd(ent.Raes.2e.curr.locat$obs_pcoa_result$li[TR.samples,i])))
      RE.IL <- (mean(ent.Raes.2e.curr.locat$obs_pcoa_result$li[RE.samples,i])-
                  mean(ent.Raes.2e.curr.locat$obs_pcoa_result$li[IL.samples,i]))/
        (sqrt(sd(ent.Raes.2e.curr.locat$obs_pcoa_result$li[IL.samples,i])*
                sd(ent.Raes.2e.curr.locat$obs_pcoa_result$li[RE.samples,i])))
      RE.TR <- (mean(ent.Raes.2e.curr.locat$obs_pcoa_result$li[RE.samples,i])-
                  mean(ent.Raes.2e.curr.locat$obs_pcoa_result$li[TR.samples,i]))/
        (sqrt(sd(ent.Raes.2e.curr.locat$obs_pcoa_result$li[TR.samples,i])*
                sd(ent.Raes.2e.curr.locat$obs_pcoa_result$li[RE.samples,i])))
      subtitle <- paste0('TR-IL: ',signif(TR.IL,digits = 2),
                         ', RE-IL: ',signif(RE.IL,digits = 2),
                         ', RE-TR: ',signif(RE.TR,digits = 2))
      subtitle.color = ifelse(max(abs(TR.IL),abs(RE.IL),abs(RE.TR))>=0.5,'red','black')
      plot(density(ent.Raes.2e.curr.locat$obs_pcoa_result$li[IL.samples,i]),col='darkred',
           xlim=c(min(ent.Raes.2e.curr.locat$obs_pcoa_result$li[,i]),max(ent.Raes.2e.curr.locat$obs_pcoa_result$li[,i])),
           ylim=c(0,1.05*max(density(ent.Raes.2e.curr.locat$obs_pcoa_result$li[IL.samples,i])$y,
                             density(ent.Raes.2e.curr.locat$obs_pcoa_result$li[TR.samples,i])$y,
                             density(ent.Raes.2e.curr.locat$obs_pcoa_result$li[RE.samples,i])$y)),
           main = paste0('PCo: ',i, ' (',
                         100*signif(ent.Raes.2e.curr.locat$obs_pcoa_result$eig[i]/sum(ent.Raes.2e.curr.locat$obs_pcoa_result$eig),digits=3),
                         '% var), '),
           xlab='PCo value', ylab = 'Density',cex.main=7,cex.lab=7,cex.axis=7,lwd=5
      )
      title(sub=subtitle, cex.sub=5,line=-83,col.sub=subtitle.color)
      lines(density(ent.Raes.2e.curr.locat$obs_pcoa_result$li[TR.samples,i]),col='darkgreen',lwd=5)
      lines(density(ent.Raes.2e.curr.locat$obs_pcoa_result$li[RE.samples,i]),col='blue3',lwd=5)
    }
    dev.off()
  }
  
  # PCA individual PCs
  if(PCA.individual){
    # draw samples' density vs PCos
    jpeg(file_out.PCA.individual,width=4000,height=3200)
    par(mfrow=c(3,3),adj=0.5,mar=c(15,17,10,2),mgp=c(11,5,0))
    for (i in draw.PCA.sep){
      TR.IL <- (mean(ent.Raes.2e.curr.locat$obs_pca_result$li[TR.samples,i])-
                  mean(ent.Raes.2e.curr.locat$obs_pca_result$li[IL.samples,i]))/
        (sqrt(sd(ent.Raes.2e.curr.locat$obs_pca_result$li[IL.samples,i])*
                sd(ent.Raes.2e.curr.locat$obs_pca_result$li[TR.samples,i])))
      RE.IL <- (mean(ent.Raes.2e.curr.locat$obs_pca_result$li[RE.samples,i])-
                  mean(ent.Raes.2e.curr.locat$obs_pca_result$li[IL.samples,i]))/
        (sqrt(sd(ent.Raes.2e.curr.locat$obs_pca_result$li[IL.samples,i])*
                sd(ent.Raes.2e.curr.locat$obs_pca_result$li[RE.samples,i])))
      RE.TR <- (mean(ent.Raes.2e.curr.locat$obs_pca_result$li[RE.samples,i])-
                  mean(ent.Raes.2e.curr.locat$obs_pca_result$li[TR.samples,i]))/
        (sqrt(sd(ent.Raes.2e.curr.locat$obs_pca_result$li[TR.samples,i])*
                sd(ent.Raes.2e.curr.locat$obs_pca_result$li[RE.samples,i])))
      subtitle <- paste0('TR-IL: ',signif(TR.IL,digits = 2),
                         ', RE-IL: ',signif(RE.IL,digits = 2),
                         ', RE-TR: ',signif(RE.TR,digits = 2))
      subtitle.color = ifelse(max(abs(TR.IL),abs(RE.IL),abs(RE.TR))>=0.5,'red','black')
      plot(density(ent.Raes.2e.curr.locat$obs_pca_result$li[IL.samples,i]),col='darkred',
           xlim=c(min(ent.Raes.2e.curr.locat$obs_pca_result$li[,i]),max(ent.Raes.2e.curr.locat$obs_pca_result$li[,i])),
           ylim=c(0,1.05*max(density(ent.Raes.2e.curr.locat$obs_pca_result$li[IL.samples,i])$y,
                             density(ent.Raes.2e.curr.locat$obs_pca_result$li[TR.samples,i])$y,
                             density(ent.Raes.2e.curr.locat$obs_pca_result$li[RE.samples,i])$y)),
           main = paste0('PC: ',i, ' (',
                         100*signif(ent.Raes.2e.curr.locat$obs_pca_result$eig[i]/sum(ent.Raes.2e.curr.locat$obs_pca_result$eig),digits=3),
                         '% var), '),
           xlab='PC value', ylab = 'Density',cex.main=7,cex.lab=7,cex.axis=7,lwd=5
      )
      title(sub=subtitle, cex.sub=5,line=-83,col.sub=subtitle.color)
      lines(density(ent.Raes.2e.curr.locat$obs_pca_result$li[TR.samples,i]),col='darkgreen',lwd=5)
      lines(density(ent.Raes.2e.curr.locat$obs_pca_result$li[RE.samples,i]),col='blue3',lwd=5)
    }
    dev.off()
  }
  
  
  # jpeg(filename = file_out2,width = 2000,height = 1600)
  # par(mfrow=c(3,3),adj=0.5)
  # legend <- 'red - ileum., green - transverse colon, purple - rectum'
  # j=2 #we see that 2nd eigen vector is better than first in separating locations
  # for (i in c(1,3:10)){
  #   a <- min(i,j)
  #   b<- max(i,j)
  #   xlim = c(min(ent.Raes.2e.curr.locat$obs_pcoa_result$li[,a])-0.1,max(ent.Raes.2e.curr.locat$obs_pcoa_result$li[,a])+0.1)
  #   ylim = c(min(ent.Raes.2e.curr.locat$obs_pcoa_result$li[,b])-0.1,max(ent.Raes.2e.curr.locat$obs_pcoa_result$li[,b])+0.1)
  #   # now make the plot symmetrical over x=0 and y=0
  #   xlim <- c(-max(abs(xlim)),max(abs(xlim)))
  #   ylim <- c(-max(abs(ylim)),max(abs(ylim)))
  #   # what PCo's we have
  #   PCos <- paste0('PCo: ',a, ' (',
  #                  100*signif(ent.Raes.2e.curr.locat$obs_pcoa_result$eig[a]/sum(ent.Raes.2e.curr.locat$obs_pcoa_result$eig),digits=3),
  #                  '% var), ',b, ' (',
  #                  100*signif(ent.Raes.2e.curr.locat$obs_pcoa_result$eig[b]/sum(ent.Raes.2e.curr.locat$obs_pcoa_result$eig),digits=3),
  #                  '% var)')
  #   
  # s.class(ent.Raes.2e.curr.locat$obs_pcoa_result$li[,c(a,b)], 
  #         fac=as.factor(revalue(samples.location,c('ileum'='IL','transverse colon'='TR','rectum'='RE'))),  
  #         col=c('red','purple','darkgreen'), 
  #         xlim=xlim, 
  #         ylim=ylim, 
  #         cpoint = 2, 
  #         grid=F, 
  #         sub=PCos, 
  #         cstar = 0, 
  #         csub = 5, 
  #         cellipse = 1,
  #         pch = 20)
  # }
  # dev.off()
  # 
  
  # this is for calculating mean and SD for PCs for each location
  # IL.samples.eig1.mean <- mean(ent.Raes.2e.curr.locat$obs_pca_result$li[IL.samples,1])
  # TR.samples.eig1.mean <- mean(ent.Raes.2e.curr.locat$obs_pca_result$li[TR.samples,1])
  # RE.samples.eig1.mean <- mean(ent.Raes.2e.curr.locat$obs_pca_result$li[RE.samples,1])
  # IL.samples.eig1.sd <- sd(ent.Raes.2e.curr.locat$obs_pca_result$li[IL.samples,1])
  # TR.samples.eig1.sd <- sd(ent.Raes.2e.curr.locat$obs_pca_result$li[TR.samples,1])
  # RE.samples.eig1.sd <- sd(ent.Raes.2e.curr.locat$obs_pca_result$li[RE.samples,1])
  # IL.samples.eig2.mean <- mean(ent.Raes.2e.curr.locat$obs_pca_result$li[IL.samples,2])
  # TR.samples.eig2.mean <- mean(ent.Raes.2e.curr.locat$obs_pca_result$li[TR.samples,2])
  # RE.samples.eig2.mean <- mean(ent.Raes.2e.curr.locat$obs_pca_result$li[RE.samples,2])
  # IL.samples.eig2.sd <- sd(ent.Raes.2e.curr.locat$obs_pca_result$li[IL.samples,2])
  # TR.samples.eig2.sd <- sd(ent.Raes.2e.curr.locat$obs_pca_result$li[TR.samples,2])
  # RE.samples.eig2.sd <- sd(ent.Raes.2e.curr.locat$obs_pca_result$li[RE.samples,2])
  # mean.pca.1 <- c(IL.samples.eig1.mean,TR.samples.eig1.mean,RE.samples.eig1.mean)
  # sd.pca.1 <- c(IL.samples.eig1.sd,TR.samples.eig1.sd,RE.samples.eig1.sd)
  # mean.pca.2 <- c(IL.samples.eig2.mean,TR.samples.eig2.mean,RE.samples.eig2.mean)
  # sd.pca.2 <- c(IL.samples.eig2.sd,TR.samples.eig2.sd,RE.samples.eig2.sd)
  # rownames(location.info.pca) <- c('IL','TR','RE')
  # location.info.pca <- signif(cbind(mean.pca.1,sd.pca.1,mean.pca.2,sd.pca.2),digits = 3)
  # print(location.info.pca)
  
}



# Compare taxa in three locations (or amplicons)
# INPUT:
# (need to open connection before running the program and close after!)
# 3 files, 3 locations (e.g. file=file2=file3 (v1v2), loc1 = '.IL', loc2 = '.TR',loc3 = '.RE')
# percent_threshold - to remove noise
# remove.outliers - p-value threshold to remove outliers farest from the other samples
# sample.cov.filter = sample.cov.filter for readQiimeSmart
# OUTPUT:
# stacked barplot with lenend 
compTaxa <- function(file1,loc1,file2,loc2,file3,loc3,percent_threshold = 0.01,
                     remove.outliers = 0.01,sample.cov.filter=0,main_ = ''){
  library(stringr)
  library(ggplot2)
  source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Functions_general_non_enterotyping.R')
  df1 <- readQiimeSmart(input_file = file1, locat = loc1, percent_threshold = percent_threshold,
                        remove.outliers = remove.outliers,sample.cov.filter=sample.cov.filter)
  df2 <- readQiimeSmart(input_file = file2, locat = loc2, percent_threshold = percent_threshold,
                        remove.outliers = remove.outliers,sample.cov.filter=sample.cov.filter)
  df3 <- readQiimeSmart(input_file = file3, locat = loc3, percent_threshold = percent_threshold, 
                        remove.outliers = remove.outliers,sample.cov.filter=sample.cov.filter)
  # appendix - this is for labels only
  locat_trim <- function(locat){
    if(locat=='.IL'){result = 'IL'}
    else if(locat=='.TR'){result = 'TR'}
    else if(locat=='.RE'){result = 'RE'}
    else if(locat=='all'){result = 'all'}
    return(result)
  }
  
  # in each column, need to obtain fraction of each bacteria
  # df1
  df1_col_sum <- colSums(df1)
  df1_norm <- df1
  for (i in 1:dim(df1)[2]){
    df1_norm[,i] <- df1[,i]/df1_col_sum[i]
  }
  # df2
  df2_col_sum <- colSums(df2)
  df2_norm <- df2
  for (i in 1:dim(df2)[2]){
    df2_norm[,i] <- df2[,i]/df2_col_sum[i]
  }
  # df3
  df3_col_sum <- colSums(df3)
  df3_norm <- df3
  for (i in 1:dim(df3)[2]){
    df3_norm[,i] <- df3[,i]/df3_col_sum[i]
  }
  
  # calculate average abundance and SD per taxon
  # get matrix taxon-average-SD
  # v1v2
  df1_taxon_aver_abund <- matrix(0,dim(df1)[1],4)
  tag_df1 <- paste0(word(file1,start=-2,sep='/'),'_',locat_trim(loc1))
  for (i in 1:dim(df1)[1]){
    df1_taxon_aver_abund[i,1] <- tag_df1
    df1_taxon_aver_abund[i,2] <- rownames(df1_norm)[i] 
    df1_taxon_aver_abund[i,3] <- mean(as.numeric(df1_norm[i,]))
    df1_taxon_aver_abund[i,4] <- sd(as.numeric(df1_norm[i,]))
  }
  # v3v4
  df2_taxon_aver_abund <- matrix(0,dim(df2)[1],4)
  tag_df2 <- paste0(word(file2,start=-2,sep='/'),'_',locat_trim(loc2))
  for (i in 1:dim(df2)[1]){
    df2_taxon_aver_abund[i,1] <- tag_df2
    df2_taxon_aver_abund[i,2] <- rownames(df2_norm)[i] 
    df2_taxon_aver_abund[i,3] <- mean(as.numeric(df2_norm[i,]))
    df2_taxon_aver_abund[i,4] <- sd(as.numeric(df2_norm[i,]))
  }
  # v5v6
  df3_taxon_aver_abund <- matrix(0,dim(df3)[1],4)
  tag_df3 <- paste0(word(file3,start=-2,sep='/'),'_',locat_trim(loc3))
  for (i in 1:dim(df3)[1]){
    df3_taxon_aver_abund[i,1] <- tag_df3
    df3_taxon_aver_abund[i,2] <- rownames(df3_norm)[i] 
    df3_taxon_aver_abund[i,3] <- mean(as.numeric(df3_norm[i,]))
    df3_taxon_aver_abund[i,4] <- sd(as.numeric(df3_norm[i,]))
  }
  # sort the matrices
  df1_taxon_aver_abund <- df1_taxon_aver_abund[order(as.numeric(df1_taxon_aver_abund[,3]),decreasing = TRUE),]
  df2_taxon_aver_abund <- df2_taxon_aver_abund[order(as.numeric(df2_taxon_aver_abund[,3]),decreasing = TRUE),]
  df3_taxon_aver_abund <- df3_taxon_aver_abund[order(as.numeric(df3_taxon_aver_abund[,3]),decreasing = TRUE),]
  # save first N elements
  N=10
  df1_taxon_aver_abund.top <- df1_taxon_aver_abund[1:N,]
  df2_taxon_aver_abund.top <- df2_taxon_aver_abund[1:N,]
  df3_taxon_aver_abund.top <- df3_taxon_aver_abund[1:N,]
  
  top.taxa.united  <- union(df3_taxon_aver_abund.top[,2],union(df1_taxon_aver_abund.top[,2],df2_taxon_aver_abund.top[,2]))
  top_taxa_aver_abund_across_df1_df2_df3 <- cbind(top.taxa.united,numeric(length(top.taxa.united)),
                                                  numeric(length(top.taxa.united)),numeric(length(top.taxa.united)))
  colnames(top_taxa_aver_abund_across_df1_df2_df3) <- c('Taxon',locat_trim(loc1),
                                                        locat_trim(loc2),locat_trim(loc3))
  top_taxa_aver_abund_across_df1_df2_df3[,2] <- df1_taxon_aver_abund[match(top.taxa.united,df1_taxon_aver_abund[,2]),3]
  top_taxa_aver_abund_across_df1_df2_df3[,3] <- df2_taxon_aver_abund[match(top.taxa.united,df2_taxon_aver_abund[,2]),3]
  top_taxa_aver_abund_across_df1_df2_df3[,4] <- df3_taxon_aver_abund[match(top.taxa.united,df3_taxon_aver_abund[,2]),3]
  
  our.legend <- character(dim(top_taxa_aver_abund_across_df1_df2_df3)[1])
  for (i in 1:length(our.legend)){
    if(word(top_taxa_aver_abund_across_df1_df2_df3[i,1],start = -1,sep=';') =='g__'){
      our.legend[i] <- word(top_taxa_aver_abund_across_df1_df2_df3[i,1],start = -2,end=-2,sep=';')
    } else{
      our.legend[i] <- word(top_taxa_aver_abund_across_df1_df2_df3[i,1],start = -1,end=-1,sep=';')
    }
  }
  
  par(cex=4,mar=c(5,5,2,15))
  barplot(top_taxa_aver_abund_across_df1_df2_df3[,2:4], col=colors()[c(23,89,12,7,52,227,145,441,373,424,432,370)], 
          border="white", space=0.75, font.axis=2, xlab="location",ylab='taxon fraction',
          cex.lab=2, cex.names = 2, legend.text = our.legend,xlim=c(0,5),
          main=main_,
          args.legend = c(cex=1,xjust=-0.2))
}


# Compare taxa in three locations - draw BOXPLOTS
# INPUT:
# file1,2,3 and loc1,2,3 make file-locations pairs to compare between (file1+loc1, file2+loc2 etc.)
# file_out = output file, must be .jpg
# 3 files, 3 locations (e.g. file=file2=file3 (v1v2), loc1 = '.IL', loc2 = '.TR',loc3 = '.RE')
# percent_threshold - to remove noise - taxa with average abundance across 
#                     all samples < x% removed (beware: 0.01 means 0.01%!)
# remove.outliers - p-value threshold to remove outliers farest from the other samples
# sample.cov.filter - for reading QIIME output, min number of reads to consider sample
# mfrow_ = c(3,3) - 3x3=9 plots per file
# mar_=c(6,6,6,2) - margins for each plot. main and axis font size correlate with them
# mgp_ = c(6,3,0) - margins for plots (main, axis title, axis text)
# taxon.abundance.threshold.plot=0.01 - we select only samples with average abundance in at least
#                                       1 location >= it
# width_=2000,height_=1600 - parameters for output files
# border_ = c('red','darkgreen','purple') - colors for locations
# min.abundance.any.sample = 0.001 - we do not plot if the upper whisker in maximal abundance of this taxon 
#                                    in any sample and location < this value. Upper whisker is 1.5*IQR+Q3, IQR is difference between quantiles Q75-Q25, Q3 = Q75
# OUTPUT:
# files with boxplots for taxa across locations
compTaxaBoxplots <- function(file1,loc1,file2,loc2,file3,loc3,file_out,
                             out_dir = file.path(home.dir,'Results/Locations_compare/Abundance_boxplots'),
                             percent_threshold = 0.01,remove.outliers = 0.01,sample.cov.filter=0,
                             mfrow_ = c(3,3),mar_=c(6,8,6,2),mgp_ = c(6,3,0), taxon.abundance.threshold.plot=0.01,
                             width_=2000,height_=1600, border_ = c('red','darkgreen','purple'),
                             min.abundance.any.sample = 0.001){
  library(stringr)
  library(ggplot2)
  source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Functions_general_non_enterotyping.R')
  df1 <- readQiimeSmart(input_file = file1, locat = loc1, percent_threshold = percent_threshold,
                        remove.outliers = remove.outliers,sample.cov.filter=sample.cov.filter)
  df2 <- readQiimeSmart(input_file = file2, locat = loc2, percent_threshold = percent_threshold,
                        remove.outliers = remove.outliers,sample.cov.filter=sample.cov.filter)
  df3 <- readQiimeSmart(input_file = file3, locat = loc3, percent_threshold = percent_threshold, 
                        remove.outliers = remove.outliers,sample.cov.filter=sample.cov.filter)
  # appendix - this is for labels only
  locat_trim <- function(locat){
    if(locat=='.IL'){result = 'IL'}
    else if(locat=='.TR'){result = 'TR'}
    else if(locat=='.RE'){result = 'RE'}
    else if(locat=='all'){result = 'all'}
    return(result)
  }
  
  # in each column, need to obtain fraction of each bacteria
  # df1
  df1_col_sum <- colSums(df1)
  df1_norm <- df1
  for (i in 1:dim(df1)[2]){
    df1_norm[,i] <- df1[,i]/df1_col_sum[i]
  }
  # df2
  df2_col_sum <- colSums(df2)
  df2_norm <- df2
  for (i in 1:dim(df2)[2]){
    df2_norm[,i] <- df2[,i]/df2_col_sum[i]
  }
  # df3
  df3_col_sum <- colSums(df3)
  df3_norm <- df3
  for (i in 1:dim(df3)[2]){
    df3_norm[,i] <- df3[,i]/df3_col_sum[i]
  }
  
  # calculate average abundance and SD per taxon
  # get matrix taxon-average-SD
  # v1v2
  df1_taxon_aver_abund <- matrix(0,dim(df1)[1],4)
  tag_df1 <- paste0(word(file1,start=-2,sep='/'),'_',locat_trim(loc1))
  for (i in 1:dim(df1)[1]){
    df1_taxon_aver_abund[i,1] <- tag_df1
    df1_taxon_aver_abund[i,2] <- rownames(df1_norm)[i] 
    df1_taxon_aver_abund[i,3] <- mean(as.numeric(df1_norm[i,]))
    df1_taxon_aver_abund[i,4] <- sd(as.numeric(df1_norm[i,]))
  }
  # v3v4
  df2_taxon_aver_abund <- matrix(0,dim(df2)[1],4)
  tag_df2 <- paste0(word(file2,start=-2,sep='/'),'_',locat_trim(loc2))
  for (i in 1:dim(df2)[1]){
    df2_taxon_aver_abund[i,1] <- tag_df2
    df2_taxon_aver_abund[i,2] <- rownames(df2_norm)[i] 
    df2_taxon_aver_abund[i,3] <- mean(as.numeric(df2_norm[i,]))
    df2_taxon_aver_abund[i,4] <- sd(as.numeric(df2_norm[i,]))
  }
  # v5v6
  df3_taxon_aver_abund <- matrix(0,dim(df3)[1],4)
  tag_df3 <- paste0(word(file3,start=-2,sep='/'),'_',locat_trim(loc3))
  for (i in 1:dim(df3)[1]){
    df3_taxon_aver_abund[i,1] <- tag_df3
    df3_taxon_aver_abund[i,2] <- rownames(df3_norm)[i] 
    df3_taxon_aver_abund[i,3] <- mean(as.numeric(df3_norm[i,]))
    df3_taxon_aver_abund[i,4] <- sd(as.numeric(df3_norm[i,]))
  }
  # sort the matrices
  df1_taxon_aver_abund <- df1_taxon_aver_abund[order(as.numeric(df1_taxon_aver_abund[,3]),decreasing = TRUE),]
  df2_taxon_aver_abund <- df2_taxon_aver_abund[order(as.numeric(df2_taxon_aver_abund[,3]),decreasing = TRUE),]
  df3_taxon_aver_abund <- df3_taxon_aver_abund[order(as.numeric(df3_taxon_aver_abund[,3]),decreasing = TRUE),]
  # now, draw the boxplots
  setwd(out_dir)
  common.taxa <- intersect(df1_taxon_aver_abund[,2],intersect(df2_taxon_aver_abund[,2],
                                                              df3_taxon_aver_abund[,2]))
  all.taxa <- union(df1_taxon_aver_abund[,2],union(df2_taxon_aver_abund[,2],
                                                   df3_taxon_aver_abund[,2]))
  first.taxa.df1 <- df1_taxon_aver_abund[df1_taxon_aver_abund[,3]>=taxon.abundance.threshold.plot,2]
  first.taxa.df2 <- df2_taxon_aver_abund[df2_taxon_aver_abund[,3]>=taxon.abundance.threshold.plot,2]
  first.taxa.df3 <- df3_taxon_aver_abund[df3_taxon_aver_abund[,3]>=taxon.abundance.threshold.plot,2]
  first.taxa.any.df <- union(first.taxa.df1,union(first.taxa.df2,first.taxa.df3))
  # Now need to rarify the first.taxa.any.df by deleting the ones where max 
  #       whisker in maximal abundance of this taxon in any sample and location < min.abundance.any.sample
  for(i in 1:length(first.taxa.any.df)){
    if(max(boxplot(as.vector(as.numeric(df1_norm[rownames(df1_norm) %in% first.taxa.any.df[i],])),
                   as.vector(as.numeric(df2_norm[rownames(df2_norm) %in% first.taxa.any.df[i],])),
                   as.vector(as.numeric(df3_norm[rownames(df3_norm) %in% first.taxa.any.df[i],])),
                   outline = F,plot = F)$stats, na.rm = T) < min.abundance.any.sample){
      first.taxa.any.df[i] <- NA
    }
  }
  first.taxa.any.df <- first.taxa.any.df[!is.na(first.taxa.any.df)]
  if(mfrow_[1]<1 | mfrow_[2]<1){stop('mfrow_ must be an integer vector')}
  if(mar_[1]<1 | mar_[2]<1 | mar_[3]<1 | mar_[4]<1){stop('mar_ must be an integer vector')}
  number.of.files = ceiling(length(first.taxa.any.df)/(mfrow_[1]*mfrow_[2]))
  file.name.root <- word(file_out,sep='\\.',start = 1,end=-2)
  for (j in 1:number.of.files){
    jpeg(paste0(file.name.root,'_',j,'.jpg'),width = width_,height = height_)
    par(mfrow = mfrow_,mar=mar_,mgp = mgp_)
    for (i in 1:(mfrow_[1]*mfrow_[2])){
      if((j-1)*(mfrow_[1]*mfrow_[2])+i<= length(first.taxa.any.df)){
        main_=word(first.taxa.any.df[(j-1)*(mfrow_[1]*mfrow_[2])+i],start = -1,sep=';')
        # in case unknown family
        if(main_=='f__'){main_<-word(first.taxa.any.df[(j-1)*(mfrow_[1]*mfrow_[2])+i],start = -2,end=-1,sep=';')}
        ylim_max = max(boxplot(as.vector(as.numeric(df1_norm[rownames(df1_norm) %in% first.taxa.any.df[(j-1)*(mfrow_[1]*mfrow_[2])+i],])),
                               as.vector(as.numeric(df2_norm[rownames(df2_norm) %in% first.taxa.any.df[(j-1)*(mfrow_[1]*mfrow_[2])+i],])),
                               as.vector(as.numeric(df3_norm[rownames(df3_norm) %in% first.taxa.any.df[(j-1)*(mfrow_[1]*mfrow_[2])+i],])),
                               outline = F,plot = F)$stats, na.rm = T)
        if(length(as.vector(as.numeric(df1_norm[rownames(df1_norm) %in% first.taxa.any.df[(j-1)*(mfrow_[1]*mfrow_[2])+i],]))[!is.na(as.vector(as.numeric(df1_norm[rownames(df1_norm) %in% first.taxa.any.df[(j-1)*(mfrow_[1]*mfrow_[2])+i],])))]!=0)&
           length(as.vector(as.numeric(df2_norm[rownames(df2_norm) %in% first.taxa.any.df[(j-1)*(mfrow_[1]*mfrow_[2])+i],]))[!is.na(as.vector(as.numeric(df2_norm[rownames(df2_norm) %in% first.taxa.any.df[(j-1)*(mfrow_[1]*mfrow_[2])+i],])))]!=0)&
           length(as.vector(as.numeric(df3_norm[rownames(df3_norm) %in% first.taxa.any.df[(j-1)*(mfrow_[1]*mfrow_[2])+i],]))[!is.na(as.vector(as.numeric(df3_norm[rownames(df3_norm) %in% first.taxa.any.df[(j-1)*(mfrow_[1]*mfrow_[2])+i],])))]!=0)){
        p_val = kruskal.test(list(as.vector(as.numeric(df1_norm[rownames(df1_norm) %in% first.taxa.any.df[(j-1)*(mfrow_[1]*mfrow_[2])+i],]))[!is.na(as.vector(as.numeric(df1_norm[rownames(df1_norm) %in% first.taxa.any.df[(j-1)*(mfrow_[1]*mfrow_[2])+i],])))],
                                  as.vector(as.numeric(df2_norm[rownames(df2_norm) %in% first.taxa.any.df[(j-1)*(mfrow_[1]*mfrow_[2])+i],]))[!is.na(as.vector(as.numeric(df2_norm[rownames(df2_norm) %in% first.taxa.any.df[(j-1)*(mfrow_[1]*mfrow_[2])+i],])))],
                                  as.vector(as.numeric(df3_norm[rownames(df3_norm) %in% first.taxa.any.df[(j-1)*(mfrow_[1]*mfrow_[2])+i],]))[!is.na(as.vector(as.numeric(df3_norm[rownames(df3_norm) %in% first.taxa.any.df[(j-1)*(mfrow_[1]*mfrow_[2])+i],])))])
                             )$p.value
        } else {p_val=NA}
        boxplot(as.vector(as.numeric(df1_norm[rownames(df1_norm) %in% first.taxa.any.df[(j-1)*(mfrow_[1]*mfrow_[2])+i],])),
                as.vector(as.numeric(df2_norm[rownames(df2_norm) %in% first.taxa.any.df[(j-1)*(mfrow_[1]*mfrow_[2])+i],])),
                as.vector(as.numeric(df3_norm[rownames(df3_norm) %in% first.taxa.any.df[(j-1)*(mfrow_[1]*mfrow_[2])+i],])),
                main=main_,cex.main=mar_[3],
                cex.axis = mar_[1]-1,
                outline = T,notch=F,names=c(locat_trim(loc1),locat_trim(loc2),locat_trim(loc3)),
                border=border_,ylim = c(0,1.01*ylim_max))
        title(sub = paste0('p-value: ',signif(p_val,digits = 2)),line=-37,cex.sub = max(mar_[1]-1,3),col.sub='grey')
      }
    }
    dev.off()
  }
}


