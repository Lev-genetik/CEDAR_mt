# Input: 
# a file1 with enterotypes: col1 - number (automatically read as rowname), col2 - sample ID, col3 - enterotype_Driver_full, col4 -enterotype number - short
# a file2 with Principal component decomposition of the gene expression data. Column 1 and 2 = sample ID
# (matching to file1), columns 3+ = PCs. 1-200 or 300 or 350 or... respectively

source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/CorrelateEnterotypesExpressionPCs_2018_12.R')

# Part 1: NOT NECESSARY AS WE HAVE >30 SAMPLES IN EACH GROUP (Leuik)
# to calculate if we can use t-test for analyzing if enterotypes influence expression
# do not need it if we use non-parametric test: Leuik says that for >30 samples, we can use t-test anyway

file_1=file.path(home.dir,'Results/Enterotyping/2018_11_30_enterotyping/Sample-enterotype/CEDAR_consensus_enterotypes_across_ampl_and_loc.csv')
checkup.IL <- corrEnterExpr.initiation(file1=file_1, file2=file.path(home.dir,'Work_dir/DATA/Momo-Dmitrieva2018/IL_PCs_after_correction_v2.txt'),locat = 'all')
checkup.TR <- corrEnterExpr.initiation(file1=file_1, file2=file.path(home.dir,'Work_dir/DATA/Momo-Dmitrieva2018/TR_PCs_after_correction_v2.txt'),locat = 'all')
checkup.RE <- corrEnterExpr.initiation(file1=file_1, file2=file.path(home.dir,'Work_dir/DATA/Momo-Dmitrieva2018/RE_PCs_after_correction_v2.txt'),locat = 'all')
checkup.CD4 <- corrEnterExpr.initiation(file1=file_1, file2=file.path(home.dir,'Work_dir/DATA/Momo-Dmitrieva2018/CD4_PCs_after_correction_v2.txt'),locat = 'all')
checkup.CD8 <- corrEnterExpr.initiation(file1=file_1, file2=file.path(home.dir,'Work_dir/DATA/Momo-Dmitrieva2018/CD8_PCs_after_correction_v2.txt'),locat = 'all')
checkup.CD14 <- corrEnterExpr.initiation(file1=file_1, file2=file.path(home.dir,'Work_dir/DATA/Momo-Dmitrieva2018/CD14_PCs_after_correction_v2.txt'),locat = 'all')
checkup.CD15 <- corrEnterExpr.initiation(file1=file_1, file2=file.path(home.dir,'Work_dir/DATA/Momo-Dmitrieva2018/CD15_PCs_after_correction_v2.txt'),locat = 'all')
checkup.CD19 <- corrEnterExpr.initiation(file1=file_1, file2=file.path(home.dir,'Work_dir/DATA/Momo-Dmitrieva2018/CD19_PCs_after_correction_v2.txt'),locat = 'all')
checkup.PLA <- corrEnterExpr.initiation(file1=file_1, file2=file.path(home.dir,'Work_dir/DATA/Momo-Dmitrieva2018/PLA_PCs_after_correction_v2.txt'),locat = 'all')
checkup.IL$can_we_use_t_test[1:3,]
checkup.TR$can_we_use_t_test[1:3,]
checkup.RE$can_we_use_t_test[1:3,]
checkup.CD4$can_we_use_t_test[1:3,]
checkup.CD8$can_we_use_t_test[1:3,]
checkup.CD14$can_we_use_t_test[1:3,]
checkup.CD15$can_we_use_t_test[1:3,]
checkup.CD19$can_we_use_t_test[1:3,]
checkup.PLA$can_we_use_t_test[1:3,]


# Part 2. The statistical test per se

# PCs
file1=file.path(home.dir,'Results/Enterotyping/2018_11_30_enterotyping/Sample-enterotype/CEDAR_consensus_enterotypes_across_ampl_and_loc.csv')
file2.IL <- file.path(home.dir,'Work_dir/DATA/Momo-Dmitrieva2018/IL_PCs_after_correction_v2.txt')
file2.TR <- file.path(home.dir,'Work_dir/DATA/Momo-Dmitrieva2018/TR_PCs_after_correction_v2.txt')
file2.RE <- file.path(home.dir,'Work_dir/DATA/Momo-Dmitrieva2018/RE_PCs_after_correction_v2.txt')
file2.CD4 <- file.path(home.dir,'Work_dir/DATA/Momo-Dmitrieva2018/CD4_PCs_after_correction_v2.txt')
file2.CD8 <- file.path(home.dir,'Work_dir/DATA/Momo-Dmitrieva2018/CD8_PCs_after_correction_v2.txt')
file2.CD14 <- file.path(home.dir,'Work_dir/DATA/Momo-Dmitrieva2018/CD14_PCs_after_correction_v2.txt')
file2.CD15 <- file.path(home.dir,'Work_dir/DATA/Momo-Dmitrieva2018/CD15_PCs_after_correction_v2.txt')
file2.CD19 <- file.path(home.dir,'Work_dir/DATA/Momo-Dmitrieva2018/CD19_PCs_after_correction_v2.txt')
file2.PLA <- file.path(home.dir,'Work_dir/DATA/Momo-Dmitrieva2018/PLA_PCs_after_correction_v2.txt')

# which test do we perform: 'wilcox' or 't.test'
ttype='t.test'
# ttype='wilcox'

# Some samples do not have enterotype - we omit it
options(na.action = 'na.omit')
# do the test
assoc.IL <- corrEnterExpr.test.n(file_1=file1,file_2=file2.IL,locat_='all',n=59,test.type=ttype)
assoc.TR <- corrEnterExpr.test.n(file_1=file1,file_2=file2.TR,locat_='all',n=50,test.type=ttype)
assoc.RE <- corrEnterExpr.test.n(file_1=file1,file_2=file2.RE,locat_='all',n=53,test.type=ttype)
assoc.CD4 <- corrEnterExpr.test.n(file_1=file1,file_2=file2.CD4,locat_='all',n=38,test.type=ttype)
assoc.CD8 <- corrEnterExpr.test.n(file_1=file1,file_2=file2.CD8,locat_='all',n=35,test.type=ttype)
assoc.CD14 <- corrEnterExpr.test.n(file_1=file1,file_2=file2.CD14,locat_='all',n=36,test.type=ttype)
assoc.CD15 <- corrEnterExpr.test.n(file_1=file1,file_2=file2.CD15,locat_='all',n=27,test.type=ttype)
assoc.CD19 <- corrEnterExpr.test.n(file_1=file1,file_2=file2.CD19,locat_='all',n=40,test.type=ttype)
assoc.PLA <- corrEnterExpr.test.n(file_1=file1,file_2=file2.PLA,locat_='all',n=23,test.type=ttype)
# Now get the minimal of the p-values for each tissue
min.pval.per.tissue <- c('IL'=min(assoc.IL$p.values),'TR'=min(assoc.TR$p.values),'RE'=min(assoc.RE$p.values),
                         'CD4'=min(assoc.CD4$p.values),'CD8'=min(assoc.CD8$p.values),'CD14'=min(assoc.CD14$p.values),
                         'CD15'=min(assoc.CD15$p.values),'CD19'=min(assoc.CD19$p.values),'PLA'=min(assoc.PLA$p.values))
min.pval.per.tissue.hochberg <- c('IL'=min(p.adjust(assoc.IL$p.values,method = 'hochberg')),
                                  'TR'=min(p.adjust(assoc.TR$p.values,method = 'hochberg')),
                                  'RE'=min(p.adjust(assoc.RE$p.values,method = 'hochberg')),
                                  'CD4'=min(p.adjust(assoc.CD4$p.values,method = 'hochberg')),
                                  'CD8'=min(p.adjust(assoc.CD8$p.values,method = 'hochberg')),
                                  'CD14'=min(p.adjust(assoc.CD14$p.values,method = 'hochberg')),
                                  'CD15'=min(p.adjust(assoc.CD15$p.values,method = 'hochberg')),
                                  'CD19'=min(p.adjust(assoc.CD19$p.values,method = 'hochberg')),
                                  'PLA'=min(p.adjust(assoc.PLA$p.values,method = 'hochberg'))
                                  )
assoc.enterotype.gene.expr.PCs <- cbind(min.pval.per.tissue,min.pval.per.tissue.hochberg)
signif(assoc.enterotype.gene.expr.PCs,digits = 3)


# Draw QQ-plots
library(lattice);
setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotype_expr_association/Ent-PCs/QQ_plots')
# jpeg('QQ-plot_for_gene_expr_PCs_vs_enterotype.jpg',width=400,height=400)
# par(nfrow=c(3,3))
qqmath(~-log10(assoc.IL$p.values),xlab = '-log10(qunif(1-x))', ylab='-log10(observed p-values)',main='IL',prepanel = prepanel.qqmathline,
       distribution=function(x){-log10(qunif(1-x))}
);
# panel.qqmathline(~-log10(assoc.IL$p.values))
qqmath(~-log10(assoc.TR$p.values),xlab = '-log10(qunif(1-x))', ylab='-log10(observed p-values)',main='TR',
       distribution=function(x){-log10(qunif(1-x))}
);
qqmath(~-log10(assoc.RE$p.values),xlab = '-log10(qunif(1-x))', ylab='-log10(observed p-values)',main='RE',
       distribution=function(x){-log10(qunif(1-x))}
);
qqmath(~-log10(assoc.CD4$p.values),xlab = '-log10(qunif(1-x))', ylab='-log10(observed p-values)',main='CD4',
       distribution=function(x){-log10(qunif(1-x))}
);
qqmath(~-log10(assoc.CD8$p.values),xlab = '-log10(qunif(1-x))', ylab='-log10(observed p-values)',main='CD8',
       distribution=function(x){-log10(qunif(1-x))}
);
qqmath(~-log10(assoc.CD14$p.values),xlab = '-log10(qunif(1-x))', ylab='-log10(observed p-values)',main='CD14',
       distribution=function(x){-log10(qunif(1-x))}
);
qqmath(~-log10(assoc.CD15$p.values),xlab = '-log10(qunif(1-x))', ylab='-log10(observed p-values)',main='CD15',
       distribution=function(x){-log10(qunif(1-x))}
);
qqmath(~-log10(assoc.CD19$p.values),xlab = '-log10(qunif(1-x))', ylab='-log10(observed p-values)',main='CD19',
       distribution=function(x){-log10(qunif(1-x))}
);
qqmath(~-log10(assoc.PLA$p.values),xlab = '-log10(qunif(1-x))', ylab='-log10(observed p-values)',main='PLA',
       distribution=function(x){-log10(qunif(1-x))}
);
dev.off()
# lines(c(0,0),c(1,1))
all.p.values <- as.numeric(c(assoc.IL$p.values,assoc.TR$p.values,assoc.RE$p.values,assoc.CD4$p.values,
                  assoc.CD8$p.values,assoc.CD14$p.values,assoc.CD15$p.values,assoc.CD19$p.values,
                  assoc.PLA$p.values))
# Now, plot all p-values for the PCs on the same plot
jpeg('tmp.jpg',width=600, height = 600)
qqmath(~-log10(all.p.values),xlab = '-log10(qunif(1-x))', ylab='-log10(observed p-values)',
       main='Principal components of gene expression', cex.lab=5,
       distribution=function(x){-log10(qunif(1-x))}
);
lines(c(0,0),c(1,1))
dev.off()
