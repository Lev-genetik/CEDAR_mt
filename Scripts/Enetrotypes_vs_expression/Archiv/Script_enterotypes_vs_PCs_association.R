# Input: 
# a file1 with enterotypes: col1 - number (automatically read as rowname), col2 - sample ID, col3 - enterotype_Driver_full, col4 -enterotype number - short
# a file2 with Principal component decomposition of the gene expression data. Column 1 and 2 = sample ID
# (matching to file1), columns 3+ = PCs. 1-200 or 300 or 350 or... respectively

source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/CorrelateEnterotypesExpressionPCs.R')

# Part 1 to calculate if we can use t-test for analyzing if enterotypes influence expression
# do not need it if we use non-parametric test

# PART1 for v12-v56 stable enterotypes in IL|TR|RE vs corresponding part of intestine expression
file_1='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype/Stable_enterotypes_v12_v56.csv'
checkup.IL.IL <- corrEnterExpr.initiation(file1=file_1, file2='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/Momo-Dmitrieva2018/IL_PCs_after_correction_v2.txt',locat = 'IL')
checkup.RE.RE <- corrEnterExpr.initiation(file1=file_1, file2='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/Momo-Dmitrieva2018/RE_PCs_after_correction_v2.txt',locat = 'RE')
checkup.TR.TR <- corrEnterExpr.initiation(file1=file_1, file2='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/Momo-Dmitrieva2018/TR_PCs_after_correction_v2.txt',locat = 'TR')
write.table(checkup.IL.IL$can_we_use_t_test, '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotype_expr_association/Ent-PCs/2018_09_06_enterotyping/IL_ent_IL_expr_PCs_test_choice.csv',sep='\t',quote = F,append = F)
write.table(checkup.RE.RE$can_we_use_t_test, '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotype_expr_association/Ent-PCs/2018_09_06_enterotyping/RE_ent_RE_expr_PCs_test_choice.csv',sep='\t',quote = F,append = F)
write.table(checkup.TR.TR$can_we_use_t_test, '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotype_expr_association/Ent-PCs/2018_09_06_enterotyping/TR_ent_TR_expr_PCs_test_choice.csv',sep='\t',quote = F,append = F)
t.test.table <- cbind(checkup.IL.IL$can_we_use_t_test[,3],checkup.TR.TR$can_we_use_t_test[,3],checkup.RE.RE$can_we_use_t_test[,3])
colnames(t.test.table) <- c('IL-IL','TR-TR','RE-RE')
t.test.table
write.table(t.test.table, '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotype_expr_association/Ent-PCs/2018_09_06_enterotyping/Ent_vs_expr_where_we can_use_t-test.csv',sep='\t',quote = F,append = F)

# PART1 for v12-v56 absolutely stable enterotypes across all amplicons
file_1='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype/Absolutely_stable_enterotypes_v12_v56_IL_TR_RE.csv'
checkup.all.IL <- corrEnterExpr.initiation(file1=file_1, file2='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/Momo-Dmitrieva2018/IL_PCs_after_correction_v2.txt',locat = 'all')
checkup.all.TR <- corrEnterExpr.initiation(file1=file_1, file2='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/Momo-Dmitrieva2018/TR_PCs_after_correction_v2.txt',locat = 'all')
checkup.all.RE <- corrEnterExpr.initiation(file1=file_1, file2='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/Momo-Dmitrieva2018/RE_PCs_after_correction_v2.txt',locat = 'all')
checkup.all.CD4 <- corrEnterExpr.initiation(file1=file_1, file2='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/Momo-Dmitrieva2018/CD4_PCs_after_correction_v2.txt',locat = 'all')
checkup.all.CD8 <- corrEnterExpr.initiation(file1=file_1, file2='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/Momo-Dmitrieva2018/CD8_PCs_after_correction_v2.txt',locat = 'all')
checkup.all.CD14 <- corrEnterExpr.initiation(file1=file_1, file2='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/Momo-Dmitrieva2018/CD14_PCs_after_correction_v2.txt',locat = 'all')
checkup.all.CD15 <- corrEnterExpr.initiation(file1=file_1, file2='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/Momo-Dmitrieva2018/CD15_PCs_after_correction_v2.txt',locat = 'all')
checkup.all.CD19 <- corrEnterExpr.initiation(file1=file_1, file2='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/Momo-Dmitrieva2018/CD19_PCs_after_correction_v2.txt',locat = 'all')
checkup.all.PLA <- corrEnterExpr.initiation(file1=file_1, file2='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/Momo-Dmitrieva2018/PLA_PCs_after_correction_v2.txt',locat = 'all')
checkup.all.IL$can_we_use_t_test[1:3,]
checkup.all.TR$can_we_use_t_test[1:3,]
checkup.all.RE$can_we_use_t_test[1:3,]
checkup.all.CD4$can_we_use_t_test[1:3,]
checkup.all.CD8$can_we_use_t_test[1:3,]
checkup.all.CD14$can_we_use_t_test[1:3,]
checkup.all.CD15$can_we_use_t_test[1:3,]
checkup.all.CD19$can_we_use_t_test[1:3,]
checkup.all.PLA$can_we_use_t_test[1:3,]


# Part 2. Mann-Whitney U test

# PART2 For v12-v56 stable enterotypes in IL|TR|RE vs corresponding part of intestine expression
# as we can see from the result, it's better to use Mann-Whitney U test
# now get the result! assoc.IL$p.values - are the p-values of the Mann-Whitney U test, we search
# for p-values < 0.01 
# assoc.IL etc. contain the results. We are interested in assoc.IL$p.values etc.
file1 <- '/media/lev-genetik/980E73270E72FD96/Liege/Results/Clustering/Sample_cluster/2018_09_06_enterotyping/Stable_enterotypes_v12_v56.csv'
# these are files by Julia
file2.IL <- '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/Momo-Dmitrieva2018/IL_PCs_after_correction_v2.txt'
file2.TR <- '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/Momo-Dmitrieva2018/TR_PCs_after_correction_v2.txt'
file2.RE <- '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/Momo-Dmitrieva2018/RE_PCs_after_correction_v2.txt'
assoc.IL <- corrEnterExpr.test.n(file_1=file1,file_2=file2.IL,locat_='IL',n=10)
assoc.TR <- corrEnterExpr.test.n(file_1=file1,file_2=file2.TR,locat_='TR',n=10)
assoc.RE <- corrEnterExpr.test.n(file_1=file1,file_2=file2.RE,locat_='RE',n=10)
# assoc.IL$p.values
# assoc.TR$p.values
# assoc.RE$p.values
write.table(round(assoc.IL$p.values,4), file = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotype_expr_association/Ent-PCs/2018_09_06_enterotyping/IL_enter_IL_expr_PCs_assoc.csv')
write.table(round(assoc.TR$p.values,4), file = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotype_expr_association/Ent-PCs/2018_09_06_enterotyping/TR_enter_TR_expr_PCs_assoc.csv')
write.table(round(assoc.RE$p.values,4), file = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotype_expr_association/Ent-PCs/2018_09_06_enterotyping/RE_enter_RE_expr_PCs_assoc.csv')


# PART2 For v12-v56 absolutely stable enterotypes across locations vs all cells' RNA expression
# PCs
file1='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype/Absolutely_stable_enterotypes_v12_v56_IL_TR_RE.csv'
file2.IL <- '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/Momo-Dmitrieva2018/IL_PCs_after_correction_v2.txt'
file2.TR <- '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/Momo-Dmitrieva2018/TR_PCs_after_correction_v2.txt'
file2.RE <- '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/Momo-Dmitrieva2018/RE_PCs_after_correction_v2.txt'
file2.CD4 <- '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/Momo-Dmitrieva2018/CD4_PCs_after_correction_v2.txt'
file2.CD8 <- '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/Momo-Dmitrieva2018/CD8_PCs_after_correction_v2.txt'
file2.CD14 <- '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/Momo-Dmitrieva2018/CD14_PCs_after_correction_v2.txt'
file2.CD15 <- '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/Momo-Dmitrieva2018/CD15_PCs_after_correction_v2.txt'
file2.CD19 <- '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/Momo-Dmitrieva2018/CD19_PCs_after_correction_v2.txt'
file2.PLA <- '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/Momo-Dmitrieva2018/PLA_PCs_after_correction_v2.txt'

# which test do we perform: 'wilcox' or 't.test'
ttype='wilcox'
assoc.IL <- corrEnterExpr.test.n(file_1=file1,file_2=file2.IL,locat_='all',n=2,test.type=ttype)
assoc.TR <- corrEnterExpr.test.n(file_1=file1,file_2=file2.TR,locat_='all',n=2,test.type=ttype)
assoc.RE <- corrEnterExpr.test.n(file_1=file1,file_2=file2.RE,locat_='all',n=2,test.type=ttype)
assoc.CD4 <- corrEnterExpr.test.n(file_1=file1,file_2=file2.CD4,locat_='all',n=2,test.type=ttype)
assoc.CD8 <- corrEnterExpr.test.n(file_1=file1,file_2=file2.CD8,locat_='all',n=2,test.type=ttype)
assoc.CD14 <- corrEnterExpr.test.n(file_1=file1,file_2=file2.CD14,locat_='all',n=2,test.type=ttype)
assoc.CD15 <- corrEnterExpr.test.n(file_1=file1,file_2=file2.CD15,locat_='all',n=2,test.type=ttype)
assoc.CD19 <- corrEnterExpr.test.n(file_1=file1,file_2=file2.CD19,locat_='all',n=2,test.type=ttype)
assoc.PLA <- corrEnterExpr.test.n(file_1=file1,file_2=file2.PLA,locat_='all',n=2,test.type=ttype)



# now get p-values for association of the enterotype with 1st PC for Bact vs Prevotella 
# (we do not consider eterotype_3: too few samples. min n=2, but we do not see the result for 
# PC 2 yet)
assoc.IL$p.values[1,1]
assoc.TR$p.values[1,1]
assoc.RE$p.values[1,1]
assoc.CD4$p.values[1,1]
assoc.CD8$p.values[1,1]
assoc.CD14$p.values[1,1]
assoc.CD15$p.values[1,1]
assoc.CD19$p.values[1,1]
assoc.PLA$p.values[1,1]


