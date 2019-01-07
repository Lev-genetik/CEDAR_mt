

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