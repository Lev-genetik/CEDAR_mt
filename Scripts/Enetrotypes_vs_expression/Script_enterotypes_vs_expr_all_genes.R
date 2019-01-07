# script for running on the server

setwd('/home/mass/ifilesets/URT/UGU_CROHN/Lev/Assoc_gene_expr_enterotype_pvalues')
# genotypes
CEDAR.genotypes.imputed <- fread('~/Research/GEN/UAG/CROHN/CEDAR_DATA/CEDAR_GENO/CEDAR_Imputed.dos',header=T)
# find out the SNP
library('RCurl') #on a server can not be downloaded
MahMoud <- getURL('https://edc.ulg.ac.be/merci/eQTL_Statistic_all_0.05_A_9d40418dca8dc55bcba1_.txt')
cyseQTLs <-  read.table(text=MahMoud, header=T, sep=' ')
write.table(cyseQTLs,file='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotype_expr_association/cys_eQTL/Momo2018_cysQTL_FDR_0_05.txt')


# Expression data not corrected for cys-eQTL
# for IL: 
# This file ~/Research/GEN/UAG/CROHN/CEDAR_DATA/CEDAR_GE/GE_Corr_4PCs_Cis/IL_GE_Corrected4_Covars_PCs.txt
# If I work on the server
library(data.table)
# below, if use fread, have problems with expr.IL[,i]
expr.IL <- read.table('~/Research/GEN/UAG/CROHN/CEDAR_DATA/CEDAR_GE/GE_Corr_4PCs_Cis/IL_GE_Corrected4_Covars_PCs.txt',header=T,sep=' ')
expr.TR <- read.table('~/Research/GEN/UAG/CROHN/CEDAR_DATA/CEDAR_GE/GE_Corr_4PCs_Cis/TR_GE_Corrected4_Covars_PCs.txt',header=T,sep=' ')
expr.RE <- read.table('~/Research/GEN/UAG/CROHN/CEDAR_DATA/CEDAR_GE/GE_Corr_4PCs_Cis/RE_GE_Corrected4_Covars_PCs.txt',header=T,sep=' ')
expr.CD4 <- read.table('~/Research/GEN/UAG/CROHN/CEDAR_DATA/CEDAR_GE/GE_Corr_4PCs_Cis/CD4_GE_Corrected4_Covars_PCs.txt',header=T,sep=' ')
expr.CD8 <- read.table('~/Research/GEN/UAG/CROHN/CEDAR_DATA/CEDAR_GE/GE_Corr_4PCs_Cis/CD8_GE_Corrected4_Covars_PCs.txt',header=T,sep=' ')
expr.CD14 <- read.table('~/Research/GEN/UAG/CROHN/CEDAR_DATA/CEDAR_GE/GE_Corr_4PCs_Cis/CD14_GE_Corrected4_Covars_PCs.txt',header=T,sep=' ')
expr.CD15 <- read.table('~/Research/GEN/UAG/CROHN/CEDAR_DATA/CEDAR_GE/GE_Corr_4PCs_Cis/CD15_GE_Corrected4_Covars_PCs.txt',header=T,sep=' ')
expr.CD19 <- read.table('~/Research/GEN/UAG/CROHN/CEDAR_DATA/CEDAR_GE/GE_Corr_4PCs_Cis/CD19_GE_Corrected4_Covars_PCs.txt',header=T,sep=' ')
expr.PLA <- read.table('~/Research/GEN/UAG/CROHN/CEDAR_DATA/CEDAR_GE/GE_Corr_4PCs_Cis/PLA_GE_Corrected4_Covars_PCs.txt',header=T,sep=' ')

enterotypes <- read.table('/home/mass/ifilesets/URT/UGU_CROHN/Lev/CEDAR_consensus_enterotypes_across_ampl_and_loc.csv',as.is = TRUE)


# for both laptop and server:
enterotypes.Bact <- enterotypes[enterotypes[,3]=='Bacteroides',1]
enterotypes.Prev <- enterotypes[enterotypes[,3]=='Prevotella',1]
# #make subset of Bact and prev - 2/3 of elements
# enterotypes.Bact.subset <- enterotypes.Bact[1:length(enterotypes.Bact)]
# enterotypes.Prev.subset <- enterotypes.Prev[1:length(enterotypes.Prev)]

# now make the t-test fo IL
association.IL <- numeric(dim(expr.IL)[2]-2)
for (i in 3:dim(expr.IL)[2]){
  association.IL[i-2] <- t.test(expr.IL[expr.IL[,2] %in% enterotypes.Bact,i],expr.IL[expr.IL[,2] %in% enterotypes.Prev,i])$p.value
}
names(association.IL) <- colnames(expr.IL)[3:dim(expr.IL)[2]]
assoc.IL.Hotchberg <- p.adjust(association.IL,method='hochberg')
association.IL.sorted <- sort(association.IL)
assoc.IL.sorted.Hotchberg <- p.adjust(association.IL.sorted,method='hochberg')
# now make the t-test fo TR
association.TR <- numeric(dim(expr.TR)[2]-2)
for (i in 3:dim(expr.TR)[2]){
  association.TR[i-2] <- t.test(expr.TR[expr.TR[,2] %in% enterotypes.Bact,i],expr.TR[expr.TR[,2] %in% enterotypes.Prev,i])$p.value
}
names(association.TR) <- colnames(expr.TR)[3:dim(expr.TR)[2]]
assoc.TR.Hotchberg <- p.adjust(association.TR,method='hochberg')
association.TR.sorted <- sort(association.TR)
assoc.TR.sorted.Hotchberg <- p.adjust(association.TR.sorted,method='hochberg')
# now make the t-test fo RE
association.RE <- numeric(dim(expr.RE)[2]-2)
for (i in 3:dim(expr.RE)[2]){
  association.RE[i-2] <- t.test(expr.RE[expr.RE[,2] %in% enterotypes.Bact,i],expr.RE[expr.RE[,2] %in% enterotypes.Prev,i])$p.value
}
names(association.RE) <- colnames(expr.RE)[3:dim(expr.RE)[2]]
assoc.RE.Hotchberg <- p.adjust(association.RE,method='hochberg')
association.RE.sorted <- sort(association.RE)
assoc.RE.sorted.Hotchberg <- p.adjust(association.RE.sorted,method='hochberg')
# now make the t-test fo CD4
association.CD4 <- numeric(dim(expr.CD4)[2]-2)
for (i in 3:dim(expr.CD4)[2]){
  association.CD4[i-2] <- t.test(expr.CD4[expr.CD4[,2] %in% enterotypes.Bact,i],expr.CD4[expr.CD4[,2] %in% enterotypes.Prev,i])$p.value
}
names(association.CD4) <- colnames(expr.CD4)[3:dim(expr.CD4)[2]]
assoc.CD4.Hotchberg <- p.adjust(association.CD4,method='hochberg')
association.CD4.sorted <- sort(association.CD4)
assoc.CD4.sorted.Hotchberg <- p.adjust(association.CD4.sorted,method='hochberg')
# now make the t-test fo CD8
association.CD8 <- numeric(dim(expr.CD8)[2]-2)
for (i in 3:dim(expr.CD8)[2]){
  association.CD8[i-2] <- t.test(expr.CD8[expr.CD8[,2] %in% enterotypes.Bact,i],expr.CD8[expr.CD8[,2] %in% enterotypes.Prev,i])$p.value
}
names(association.CD8) <- colnames(expr.CD8)[3:dim(expr.CD8)[2]]
assoc.CD8.Hotchberg <- p.adjust(association.CD8,method='hochberg')
association.CD8.sorted <- sort(association.CD8)
assoc.CD8.sorted.Hotchberg <- p.adjust(association.CD8.sorted,method='hochberg')
# now make the t-test fo CD14
association.CD14 <- numeric(dim(expr.CD14)[2]-2)
for (i in 3:dim(expr.CD14)[2]){
  association.CD14[i-2] <- t.test(expr.CD14[expr.CD14[,2] %in% enterotypes.Bact,i],expr.CD14[expr.CD14[,2] %in% enterotypes.Prev,i])$p.value
}
names(association.CD14) <- colnames(expr.CD14)[3:dim(expr.CD14)[2]]
assoc.CD14.Hotchberg <- p.adjust(association.CD14,method='hochberg')
association.CD14.sorted <- sort(association.CD14)
assoc.CD14.sorted.Hotchberg <- p.adjust(association.CD14.sorted,method='hochberg')
# now make the t-test fo CD15
association.CD15 <- numeric(dim(expr.CD15)[2]-2)
for (i in 3:dim(expr.CD15)[2]){
  association.CD15[i-2] <- t.test(expr.CD15[expr.CD15[,2] %in% enterotypes.Bact,i],expr.CD15[expr.CD15[,2] %in% enterotypes.Prev,i])$p.value
}
names(association.CD15) <- colnames(expr.CD15)[3:dim(expr.CD15)[2]]
assoc.CD15.Hotchberg <- p.adjust(association.CD15,method='hochberg')
association.CD15.sorted <- sort(association.CD15)
assoc.CD15.sorted.Hotchberg <- p.adjust(association.CD15.sorted,method='hochberg')
# now make the t-test fo CD19
association.CD19 <- numeric(dim(expr.CD19)[2]-2)
for (i in 3:dim(expr.CD19)[2]){
  association.CD19[i-2] <- t.test(expr.CD19[expr.CD19[,2] %in% enterotypes.Bact,i],expr.CD19[expr.CD19[,2] %in% enterotypes.Prev,i])$p.value
}
names(association.CD19) <- colnames(expr.CD19)[3:dim(expr.CD19)[2]]
assoc.CD19.Hotchberg <- p.adjust(association.CD19,method='hochberg')
association.CD19.sorted <- sort(association.CD19)
assoc.CD19.sorted.Hotchberg <- p.adjust(association.CD19.sorted,method='hochberg')
# now make the t-test fo PLA
association.PLA <- numeric(dim(expr.PLA)[2]-2)
for (i in 3:dim(expr.PLA)[2]){
  association.PLA[i-2] <- t.test(expr.PLA[expr.PLA[,2] %in% enterotypes.Bact,i],expr.PLA[expr.PLA[,2] %in% enterotypes.Prev,i])$p.value
}
names(association.PLA) <- colnames(expr.PLA)[3:dim(expr.PLA)[2]]
assoc.PLA.Hotchberg <- p.adjust(association.PLA,method='hochberg')
association.PLA.sorted <- sort(association.PLA)
assoc.PLA.sorted.Hotchberg <- p.adjust(association.PLA.sorted,method='hochberg')

writeLines(as.character(association.IL),con='IL.csv')
writeLines(as.character(association.TR),con='TR.csv')
writeLines(as.character(association.RE),con='RE.csv')
writeLines(as.character(association.CD4),con='CD4.csv')
writeLines(as.character(association.CD8),con='CD8.csv')
writeLines(as.character(association.CD14),con='CD14.csv')
writeLines(as.character(association.CD15),con='CD15.csv')
writeLines(as.character(association.CD19),con='CD19.csv')
writeLines(as.character(association.PLA),con='PLA.csv')

# get number of samples and genes
paste0('IL samples: ',sum(expr.IL[,2] %in% enterotypes.Bact)+sum(expr.IL[,2] %in% enterotypes.Prev))
paste0('IL genes: ',dim(expr.IL)[2])
paste0('TR samples: ',sum(expr.TR[,2] %in% enterotypes.Bact)+sum(expr.TR[,2] %in% enterotypes.Prev))
paste0('TR genes: ',dim(expr.TR)[2])
paste0('RE samples: ',sum(expr.RE[,2] %in% enterotypes.Bact)+sum(expr.RE[,2] %in% enterotypes.Prev))
paste0('RE genes: ',dim(expr.RE)[2])

paste0('CD4 samples: ',sum(expr.CD4[,2] %in% enterotypes.Bact)+sum(expr.CD4[,2] %in% enterotypes.Prev))
paste0('CD4 genes: ',dim(expr.CD4)[2])
paste0('CD8 samples: ',sum(expr.CD8[,2] %in% enterotypes.Bact)+sum(expr.CD8[,2] %in% enterotypes.Prev))
paste0('CD8 genes: ',dim(expr.CD8)[2])
paste0('CD14 samples: ',sum(expr.CD14[,2] %in% enterotypes.Bact)+sum(expr.CD14[,2] %in% enterotypes.Prev))
paste0('CD14 genes: ',dim(expr.CD14)[2])
paste0('CD15 samples: ',sum(expr.CD15[,2] %in% enterotypes.Bact)+sum(expr.CD15[,2] %in% enterotypes.Prev))
paste0('CD15 genes: ',dim(expr.CD15)[2])
paste0('CD19 samples: ',sum(expr.CD19[,2] %in% enterotypes.Bact)+sum(expr.CD19[,2] %in% enterotypes.Prev))
paste0('CD19 genes: ',dim(expr.CD19)[2])
paste0('PLA samples: ',sum(expr.PLA[,2] %in% enterotypes.Bact)+sum(expr.PLA[,2] %in% enterotypes.Prev))
paste0('PLA genes: ',dim(expr.PLA)[2])

# Do on laptop after copying files
setwd(file.path(home.dir,'Results/Enterotype_expr_association/Ent-Gene-Expr-9celltypes-P_values'))
IL.pval.gene.expr <- as.numeric(readLines('IL.csv'))
TR.pval.gene.expr <- as.numeric(readLines('TR.csv'))
RE.pval.gene.expr <- as.numeric(readLines('RE.csv'))
CD4.pval.gene.expr <- as.numeric(readLines('CD4.csv'))
CD8.pval.gene.expr <- as.numeric(readLines('CD8.csv'))
CD14.pval.gene.expr <- as.numeric(readLines('CD14.csv'))
CD15.pval.gene.expr <- as.numeric(readLines('CD15.csv'))
CD19.pval.gene.expr <- as.numeric(readLines('CD19.csv'))
PLA.pval.gene.expr <- as.numeric(readLines('PLA.csv'))

source(file.path(progr.dir,'From_web/QQunif.R'))
qqunif.plot(IL.pval.gene.expr,conf.col="lightgray", conf.alpha=.01,main='IL')
qqunif.plot(TR.pval.gene.expr,conf.col="lightgray", conf.alpha=.01,main='TR')
qqunif.plot(RE.pval.gene.expr,conf.col="lightgray", conf.alpha=.01,main='RE')
qqunif.plot(CD4.pval.gene.expr,conf.col="lightgray", conf.alpha=.01,main='CD4')
qqunif.plot(CD8.pval.gene.expr,conf.col="lightgray", conf.alpha=.01,main='CD8')
qqunif.plot(CD14.pval.gene.expr,conf.col="lightgray", conf.alpha=.01,main='CD14')
qqunif.plot(CD15.pval.gene.expr,conf.col="lightgray", conf.alpha=.01,main='CD15')
qqunif.plot(CD19.pval.gene.expr,conf.col="lightgray", conf.alpha=.01,main='CD19')
qqunif.plot(PLA.pval.gene.expr,conf.col="lightgray", conf.alpha=.01,main='PLA')
dev.off()



signif(min(association.IL),digits = 3)
signif(min(association.TR),digits = 3)
signif(min(association.RE),digits = 3)
signif(min(association.CD4),digits = 3)
signif(min(association.CD8),digits = 3)

                            signif(min(association.CD14),digits = 3)
                                   signif(min(association.CD15),digits = 3)
                                          signif(min(association.CD19),digits = 3)
                                                 signif(min(association.PLA),digits = 3)

                                                 signif(min(assoc.IL.sorted.Hotchberg),digits = 3)
                                                 signif(min(assoc.TR.sorted.Hotchberg),digits = 3)
                                                 signif(min(assoc.RE.sorted.Hotchberg),digits = 3)
                                                 signif(min(assoc.CD4.sorted.Hotchberg),digits = 3)
                                                 signif(min(assoc.CD8.sorted.Hotchberg),digits = 3)
                                                 signif(min(assoc.CD14.sorted.Hotchberg),digits = 3)
                                                 signif(min(assoc.CD15.sorted.Hotchberg),digits = 3)
                                                 signif(min(assoc.CD19.sorted.Hotchberg),digits = 3)
                                                 signif(min(assoc.PLA.sorted.Hotchberg),digits = 3)

                                                 # qqmath(~-log10(association.PLA),xlab = '-log10(qunif(1-x))', ylab='-log10(observed p-values)',main='IL',
                                                 #        distribution=function(x){-log10(qunif(1-x))}
                                                 # );
                                                 # 
                                                 
                                                 
                                                 
                                                 # # if I work on my laptop:
# setwd(file.path(home.dir,'/Work_dir/Expression_enterotype/'))
# expr.IL <- read.table(file.path(home.dir,'/Work_dir/Expression_enterotype/IL_GE_Corrected4_Covars_PCs.txt'),header=T,sep=' ')
# enterotypes <- read.table(file.path(home.dir,'Results/Enterotyping/2018_11_30_enterotyping/Sample-enterotype/CEDAR_consensus_enterotypes_across_ampl_and_loc.csv'),as.is = TRUE)

                                                 
                                                 
                                                 

                                                 
                                                                                                  
                                                 