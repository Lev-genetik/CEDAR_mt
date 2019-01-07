# Check how well each amplicon to predicts enterotypes obtained by other 2 amplicons
# 2 to 5 enterotypes
# A. all locations enterotyped together
# 1. percentage of the prediction consistency
# 2. p-value of the prediction consistency
# B. IL, TR and RE separately
# 1. percentage of the prediction consistency 
# 2. p-value of the prediction consistency


source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/WhichAmpliconPredictsOthersBetter.R')
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Enterotyping_compar_functions_v5.R')

# A. all locations enterotyped together 
setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype')

# 1. Get percentages of identical samples
# first, make the calcuation of the identity - we get both numbers and percentages
pred.val.all_loc_together.2e <- getPredNumber('V1V2_2e_sample-enterotype.csv','V3V4_2e_sample-enterotype.csv','V5V6_2e_sample-enterotype.csv')
pred.val.all_loc_together.3e <- getPredNumber('V1V2_3e_sample-enterotype.csv','V3V4_3e_sample-enterotype.csv','V5V6_3e_sample-enterotype.csv')
pred.val.all_loc_together.4e <- getPredNumber('V1V2_4e_sample-enterotype.csv','V3V4_4e_sample-enterotype.csv','V5V6_4e_sample-enterotype.csv')
pred.val.all_loc_together.5e <- getPredNumber('V1V2_5e_sample-enterotype.csv','V3V4_5e_sample-enterotype.csv','V5V6_5e_sample-enterotype.csv')
# for the final table
pred.val.percent.all_loc_together.2e <- pred.val.all_loc_together.2e$ident.percent
pred.val.percent.all_loc_together.3e <- pred.val.all_loc_together.3e$ident.percent
pred.val.percent.all_loc_together.4e <- pred.val.all_loc_together.4e$ident.percent
pred.val.percent.all_loc_together.5e <- pred.val.all_loc_together.5e$ident.percent
# these are technical for all locations taken together
pred.val.number.all_loc_together.2e <- pred.val.all_loc_together.2e$ident.number
pred.val.number.all_loc_together.3e <- pred.val.all_loc_together.3e$ident.number
pred.val.number.all_loc_together.4e <- pred.val.all_loc_together.4e$ident.number
pred.val.number.all_loc_together.5e <- pred.val.all_loc_together.5e$ident.number
pred.val.dolya.all_loc_together.2e <- pred.val.number.all_loc_together.2e/pred.val.all_loc_together.2e$total.sample.number
pred.val.dolya.all_loc_together.3e <- pred.val.number.all_loc_together.3e/pred.val.all_loc_together.3e$total.sample.number
pred.val.dolya.all_loc_together.4e <- pred.val.number.all_loc_together.4e/pred.val.all_loc_together.4e$total.sample.number
pred.val.dolya.all_loc_together.5e <- pred.val.number.all_loc_together.5e/pred.val.all_loc_together.5e$total.sample.number
# calculate all loc taken together
pred.val.percent.all_loc_together.mean.v12 <- paste0(round(100*mean(c(pred.val.dolya.all_loc_together.2e[1],pred.val.dolya.all_loc_together.3e[1],pred.val.dolya.all_loc_together.4e[1],pred.val.dolya.all_loc_together.5e[1])),2),'%')
pred.val.percent.all_loc_together.mean.v34 <- paste0(round(100*mean(c(pred.val.dolya.all_loc_together.2e[2],pred.val.dolya.all_loc_together.3e[2],pred.val.dolya.all_loc_together.4e[3],pred.val.dolya.all_loc_together.5e[1])),2),'%')
pred.val.percent.all_loc_together.mean.v56 <- paste0(round(100*mean(c(pred.val.dolya.all_loc_together.2e[3],pred.val.dolya.all_loc_together.3e[3],pred.val.dolya.all_loc_together.4e[3],pred.val.dolya.all_loc_together.5e[1])),2),'%')
# to calculate the mean for all locations and all enterotype numbers, we need to make one vector:
pred.val.percent.all_loc_together.v12v34v56 <- c(pred.val.dolya.all_loc_together.2e[4],pred.val.dolya.all_loc_together.3e[4],pred.val.dolya.all_loc_together.4e[4],pred.val.dolya.all_loc_together.5e[4])
# and now calculate it
pred.val.percent.all_loc_together.mean.v12v34v56 <- paste0(round(100*mean(pred.val.percent.all_loc_together.v12v34v56),2),'%')
# the last row of the table here:
pred.val.percent.all_loc_together_all_ent <- c(pred.val.percent.all_loc_together.mean.v12,pred.val.percent.all_loc_together.mean.v34,pred.val.percent.all_loc_together.mean.v56,pred.val.percent.all_loc_together.mean.v12v34v56)
# Now, make the table
result.all.loc.together.percent <- rbind(pred.val.percent.all_loc_together.2e,pred.val.percent.all_loc_together.3e,pred.val.percent.all_loc_together.4e,pred.val.percent.all_loc_together.5e,pred.val.percent.all_loc_together_all_ent)
colnames(result.all.loc.together.percent) <- c('v1v2','v3v4','v5v6','average')
rownames(result.all.loc.together.percent) <- c('2 enterotypes','3 enterotypes','4 enterotypes','5 enterotypes','average')

# 2. Get p-values
setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype')
pred.val.pval.all_loc_together.2e <- getPredValue('V1V2_2e_sample-enterotype.csv','V3V4_2e_sample-enterotype.csv','V5V6_2e_sample-enterotype.csv')
pred.val.pval.all_loc_together.3e <- getPredValue('V1V2_3e_sample-enterotype.csv','V3V4_3e_sample-enterotype.csv','V5V6_3e_sample-enterotype.csv')
pred.val.pval.all_loc_together.4e <- getPredValue('V1V2_4e_sample-enterotype.csv','V3V4_4e_sample-enterotype.csv','V5V6_4e_sample-enterotype.csv')
pred.val.pval.all_loc_together.5e <- getPredValue('V1V2_5e_sample-enterotype.csv','V3V4_5e_sample-enterotype.csv','V5V6_5e_sample-enterotype.csv')
# calculate geom means by enterotype
gm.all_loc_together.2e <- exp(mean(log(pred.val.pval.all_loc_together.2e)))
gm.all_loc_together.3e <- exp(mean(log(pred.val.pval.all_loc_together.3e)))
gm.all_loc_together.4e <- exp(mean(log(pred.val.pval.all_loc_together.4e)))
gm.all_loc_together.5e <- exp(mean(log(pred.val.pval.all_loc_together.5e)))
# calculate geom means by location
gm.all_loc_together.v12 <- exp(mean(log(c(pred.val.pval.all_loc_together.2e[1],pred.val.pval.all_loc_together.3e[1],pred.val.pval.all_loc_together.4e[1],pred.val.pval.all_loc_together.5e[1]))))
gm.all_loc_together.v34 <- exp(mean(log(c(pred.val.pval.all_loc_together.2e[2],pred.val.pval.all_loc_together.3e[2],pred.val.pval.all_loc_together.4e[2],pred.val.pval.all_loc_together.5e[2]))))
gm.all_loc_together.v56 <- exp(mean(log(c(pred.val.pval.all_loc_together.2e[3],pred.val.pval.all_loc_together.3e[3],pred.val.pval.all_loc_together.4e[3],pred.val.pval.all_loc_together.5e[3]))))
# now, add geom means by enterotype to the rows of the future table
result.2e.all_loc_together <- c(pred.val.pval.all_loc_together.2e,gm.all_loc_together.2e)
result.3e.all_loc_together <- c(pred.val.pval.all_loc_together.3e,gm.all_loc_together.3e)
result.4e.all_loc_together <- c(pred.val.pval.all_loc_together.4e,gm.all_loc_together.4e)
result.5e.all_loc_together <- c(pred.val.pval.all_loc_together.5e,gm.all_loc_together.5e)
# this is the geom mean of geom means - with it, we can estimate work of the model
gm.all_loc_together.gm <- exp(mean(log(c(gm.all_loc_together.2e,gm.all_loc_together.3e,gm.all_loc_together.4e,gm.all_loc_together.5e))))
# combine the vectors to a result table
result.all_loc_together.pval <- rbind(result.2e.all_loc_together,result.3e.all_loc_together,result.4e.all_loc_together,result.5e.all_loc_together,c(gm.all_loc_together.v12,gm.all_loc_together.v34,gm.all_loc_together.v56,gm.all_loc_together.gm))
result.all_loc_together.pval <- signif(result.all_loc_together.pval,digits = 3)
colnames(result.all_loc_together.pval) <- c('v1v2','v3v4','v5v6','geom.mean.locat')
rownames(result.all_loc_together.pval) <- c('2 enterotypes','3 enterotypes','4 enterotypes','5 enterotypes','geom.mean.enterot')

# 3. write the results of 1. and 2. to files
setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Compare_enterotypings/Comp_amplic_and_locat_pred_value')
write.table(result.all.loc.together.percent, 'all_enterotypes_v12_v34_v56_reproducibility_percent.csv',sep='\t')
write.table(result.all_loc_together.pval, 'all_enterotypes_v12_v34_v56_reproducibility_pval.csv',sep='\t')

# B. all locations enterotyped separately
# 1. Compare average identity percentage between a particular amplicon enterotyping and enterotyping 
# by the other 2 amplicons
setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype')

# 2 enterotypes
pred.percent.2e.IL <- getPredNumber('V1V2_2e_IL_sample-enterotype.csv','V3V4_2e_IL_sample-enterotype.csv','V5V6_2e_IL_sample-enterotype.csv')
pred.percent.2e.TR <- getPredNumber('V1V2_2e_TR_sample-enterotype.csv','V3V4_2e_TR_sample-enterotype.csv','V5V6_2e_TR_sample-enterotype.csv')
pred.percent.2e.RE <- getPredNumber('V1V2_2e_RE_sample-enterotype.csv','V3V4_2e_RE_sample-enterotype.csv','V5V6_2e_RE_sample-enterotype.csv')

# 3 enterotypes
pred.percent.3e.IL <- getPredNumber('V1V2_3e_IL_sample-enterotype.csv','V3V4_3e_IL_sample-enterotype.csv','V5V6_3e_IL_sample-enterotype.csv')
pred.percent.3e.TR <- getPredNumber('V1V2_3e_TR_sample-enterotype.csv','V3V4_3e_TR_sample-enterotype.csv','V5V6_3e_TR_sample-enterotype.csv')
pred.percent.3e.RE <- getPredNumber('V1V2_3e_RE_sample-enterotype.csv','V3V4_3e_RE_sample-enterotype.csv','V5V6_3e_RE_sample-enterotype.csv')

# 4 enterotypes
pred.percent.4e.IL <- getPredNumber('V1V2_4e_IL_sample-enterotype.csv','V3V4_4e_IL_sample-enterotype.csv','V5V6_4e_IL_sample-enterotype.csv')
pred.percent.4e.TR <- getPredNumber('V1V2_4e_TR_sample-enterotype.csv','V3V4_4e_TR_sample-enterotype.csv','V5V6_4e_TR_sample-enterotype.csv')
pred.percent.4e.RE <- getPredNumber('V1V2_4e_RE_sample-enterotype.csv','V3V4_4e_RE_sample-enterotype.csv','V5V6_4e_RE_sample-enterotype.csv')

# 5 enterotypes
pred.percent.5e.IL <- getPredNumber('V1V2_5e_IL_sample-enterotype.csv','V3V4_5e_IL_sample-enterotype.csv','V5V6_5e_IL_sample-enterotype.csv')
pred.percent.5e.TR <- getPredNumber('V1V2_5e_TR_sample-enterotype.csv','V3V4_5e_TR_sample-enterotype.csv','V5V6_5e_TR_sample-enterotype.csv')
pred.percent.5e.RE <- getPredNumber('V1V2_5e_RE_sample-enterotype.csv','V3V4_5e_RE_sample-enterotype.csv','V5V6_5e_RE_sample-enterotype.csv')

# Now calculate average identity percentage per amplicon
aver.pred.percent.2e.v12 <- paste0(round(100*sum(pred.percent.2e.IL$ident.number[1],pred.percent.2e.TR$ident.number[1],pred.percent.2e.RE$ident.number[1])/sum(pred.percent.2e.IL$total.sample.number,pred.percent.2e.TR$total.sample.number,pred.percent.2e.RE$total.sample.number),2),'%')
aver.pred.percent.2e.v34 <- paste0(round(100*sum(pred.percent.2e.IL$ident.number[2],pred.percent.2e.TR$ident.number[2],pred.percent.2e.RE$ident.number[2])/sum(pred.percent.2e.IL$total.sample.number,pred.percent.2e.TR$total.sample.number,pred.percent.2e.RE$total.sample.number),2),'%')
aver.pred.percent.2e.v56 <- paste0(round(100*sum(pred.percent.2e.IL$ident.number[3],pred.percent.2e.TR$ident.number[3],pred.percent.2e.RE$ident.number[3])/sum(pred.percent.2e.IL$total.sample.number,pred.percent.2e.TR$total.sample.number,pred.percent.2e.RE$total.sample.number),2),'%')
aver.pred.percent.2e.all_amp <- paste0(round(100*sum(pred.percent.2e.IL$ident.number[4],pred.percent.2e.TR$ident.number[4],pred.percent.2e.RE$ident.number[4])/sum(pred.percent.2e.IL$total.sample.number,pred.percent.2e.TR$total.sample.number,pred.percent.2e.RE$total.sample.number),2),'%')

aver.pred.percent.3e.v12 <- paste0(round(100*sum(pred.percent.3e.IL$ident.number[1],pred.percent.3e.TR$ident.number[1],pred.percent.3e.RE$ident.number[1])/sum(pred.percent.3e.IL$total.sample.number,pred.percent.3e.TR$total.sample.number,pred.percent.3e.RE$total.sample.number),2),'%')
aver.pred.percent.3e.v34 <- paste0(round(100*sum(pred.percent.3e.IL$ident.number[2],pred.percent.3e.TR$ident.number[2],pred.percent.3e.RE$ident.number[2])/sum(pred.percent.3e.IL$total.sample.number,pred.percent.3e.TR$total.sample.number,pred.percent.3e.RE$total.sample.number),2),'%')
aver.pred.percent.3e.v56 <- paste0(round(100*sum(pred.percent.3e.IL$ident.number[3],pred.percent.3e.TR$ident.number[3],pred.percent.3e.RE$ident.number[3])/sum(pred.percent.3e.IL$total.sample.number,pred.percent.3e.TR$total.sample.number,pred.percent.3e.RE$total.sample.number),2),'%')
aver.pred.percent.3e.all_amp <- paste0(round(100*sum(pred.percent.3e.IL$ident.number[4],pred.percent.3e.TR$ident.number[4],pred.percent.3e.RE$ident.number[4])/sum(pred.percent.3e.IL$total.sample.number,pred.percent.3e.TR$total.sample.number,pred.percent.3e.RE$total.sample.number),2),'%')

aver.pred.percent.4e.v12 <- paste0(round(100*sum(pred.percent.4e.IL$ident.number[1],pred.percent.4e.TR$ident.number[1],pred.percent.4e.RE$ident.number[1])/sum(pred.percent.4e.IL$total.sample.number,pred.percent.4e.TR$total.sample.number,pred.percent.4e.RE$total.sample.number),2),'%')
aver.pred.percent.4e.v34 <- paste0(round(100*sum(pred.percent.4e.IL$ident.number[2],pred.percent.4e.TR$ident.number[2],pred.percent.4e.RE$ident.number[2])/sum(pred.percent.4e.IL$total.sample.number,pred.percent.4e.TR$total.sample.number,pred.percent.4e.RE$total.sample.number),2),'%')
aver.pred.percent.4e.v56 <- paste0(round(100*sum(pred.percent.4e.IL$ident.number[3],pred.percent.4e.TR$ident.number[3],pred.percent.4e.RE$ident.number[3])/sum(pred.percent.4e.IL$total.sample.number,pred.percent.4e.TR$total.sample.number,pred.percent.4e.RE$total.sample.number),2),'%')
aver.pred.percent.4e.all_amp <- paste0(round(100*sum(pred.percent.4e.IL$ident.number[4],pred.percent.4e.TR$ident.number[4],pred.percent.4e.RE$ident.number[4])/sum(pred.percent.4e.IL$total.sample.number,pred.percent.4e.TR$total.sample.number,pred.percent.4e.RE$total.sample.number),2),'%')

aver.pred.percent.5e.v12 <- paste0(round(100*sum(pred.percent.5e.IL$ident.number[1],pred.percent.5e.TR$ident.number[1],pred.percent.5e.RE$ident.number[1])/sum(pred.percent.5e.IL$total.sample.number,pred.percent.5e.TR$total.sample.number,pred.percent.5e.RE$total.sample.number),2),'%')
aver.pred.percent.5e.v34 <- paste0(round(100*sum(pred.percent.5e.IL$ident.number[2],pred.percent.5e.TR$ident.number[2],pred.percent.5e.RE$ident.number[2])/sum(pred.percent.5e.IL$total.sample.number,pred.percent.5e.TR$total.sample.number,pred.percent.5e.RE$total.sample.number),2),'%')
aver.pred.percent.5e.v56 <- paste0(round(100*sum(pred.percent.5e.IL$ident.number[3],pred.percent.5e.TR$ident.number[3],pred.percent.5e.RE$ident.number[3])/sum(pred.percent.5e.IL$total.sample.number,pred.percent.5e.TR$total.sample.number,pred.percent.5e.RE$total.sample.number),2),'%')
aver.pred.percent.5e.all_amp <- paste0(round(100*sum(pred.percent.5e.IL$ident.number[4],pred.percent.5e.TR$ident.number[4],pred.percent.5e.RE$ident.number[4])/sum(pred.percent.5e.IL$total.sample.number,pred.percent.5e.TR$total.sample.number,pred.percent.5e.RE$total.sample.number),2),'%')

# prepare a results' tables: rows = IL,TR,RE and average,
# columns = v1v2, v3v4, v5v6 and average
result.sep.perc.2e <- rbind(pred.percent.2e.IL$ident.percent,pred.percent.2e.TR$ident.percent,pred.percent.2e.RE$ident.percent,c(aver.pred.percent.2e.v12,aver.pred.percent.2e.v34,aver.pred.percent.2e.v56,aver.pred.percent.2e.all_amp))
colnames(result.sep.perc.2e) <- c('v1v2','v3v4','v5v6','average.locat')
rownames(result.sep.perc.2e) <- c('IL','TR','RE','average.locat')

result.sep.perc.3e <- rbind(pred.percent.3e.IL$ident.percent,pred.percent.3e.TR$ident.percent,pred.percent.3e.RE$ident.percent,c(aver.pred.percent.3e.v12,aver.pred.percent.3e.v34,aver.pred.percent.3e.v56,aver.pred.percent.3e.all_amp))
colnames(result.sep.perc.3e) <- c('v1v2','v3v4','v5v6','average.locat')
rownames(result.sep.perc.3e) <- c('IL','TR','RE','average.locat')

result.sep.perc.4e <- rbind(pred.percent.4e.IL$ident.percent,pred.percent.4e.TR$ident.percent,pred.percent.4e.RE$ident.percent,c(aver.pred.percent.4e.v12,aver.pred.percent.4e.v34,aver.pred.percent.4e.v56,aver.pred.percent.4e.all_amp))
colnames(result.sep.perc.4e) <- c('v1v2','v3v4','v5v6','average.locat')
rownames(result.sep.perc.4e) <- c('IL','TR','RE','average.locat')

result.sep.perc.5e <- rbind(pred.percent.5e.IL$ident.percent,pred.percent.5e.TR$ident.percent,pred.percent.5e.RE$ident.percent,c(aver.pred.percent.5e.v12,aver.pred.percent.5e.v34,aver.pred.percent.5e.v56,aver.pred.percent.5e.all_amp))
colnames(result.sep.perc.5e) <- c('v1v2','v3v4','v5v6','average.locat')
rownames(result.sep.perc.5e) <- c('IL','TR','RE','average.locat')

# 2. Compare average identity p-value between a particular amplicon enterotyping and enterotyping 
# by the other 2 amplicons
setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype')

# 2 enterotypes
pred.val.2e.IL <- getPredValue('V1V2_2e_IL_sample-enterotype.csv','V3V4_2e_IL_sample-enterotype.csv','V5V6_2e_IL_sample-enterotype.csv')
pred.val.2e.TR <- getPredValue('V1V2_2e_TR_sample-enterotype.csv','V3V4_2e_TR_sample-enterotype.csv','V5V6_2e_TR_sample-enterotype.csv')
pred.val.2e.RE <- getPredValue('V1V2_2e_RE_sample-enterotype.csv','V3V4_2e_RE_sample-enterotype.csv','V5V6_2e_RE_sample-enterotype.csv')

# 3 enterotypes
pred.val.3e.IL <- getPredValue('V1V2_3e_IL_sample-enterotype.csv','V3V4_3e_IL_sample-enterotype.csv','V5V6_3e_IL_sample-enterotype.csv')
pred.val.3e.TR <- getPredValue('V1V2_3e_TR_sample-enterotype.csv','V3V4_3e_TR_sample-enterotype.csv','V5V6_3e_TR_sample-enterotype.csv')
pred.val.3e.RE <- getPredValue('V1V2_3e_RE_sample-enterotype.csv','V3V4_3e_RE_sample-enterotype.csv','V5V6_3e_RE_sample-enterotype.csv')

# 4 enterotypes
pred.val.4e.IL <- getPredValue('V1V2_4e_IL_sample-enterotype.csv','V3V4_4e_IL_sample-enterotype.csv','V5V6_4e_IL_sample-enterotype.csv')
pred.val.4e.TR <- getPredValue('V1V2_4e_TR_sample-enterotype.csv','V3V4_4e_TR_sample-enterotype.csv','V5V6_4e_TR_sample-enterotype.csv')
pred.val.4e.RE <- getPredValue('V1V2_4e_RE_sample-enterotype.csv','V3V4_4e_RE_sample-enterotype.csv','V5V6_4e_RE_sample-enterotype.csv')

# 5 enterotypes
pred.val.5e.IL <- getPredValue('V1V2_5e_IL_sample-enterotype.csv','V3V4_5e_IL_sample-enterotype.csv','V5V6_5e_IL_sample-enterotype.csv')
pred.val.5e.TR <- getPredValue('V1V2_5e_TR_sample-enterotype.csv','V3V4_5e_TR_sample-enterotype.csv','V5V6_5e_TR_sample-enterotype.csv')
pred.val.5e.RE <- getPredValue('V1V2_5e_RE_sample-enterotype.csv','V3V4_5e_RE_sample-enterotype.csv','V5V6_5e_RE_sample-enterotype.csv')

# Now calculate geometric mean of the prediction values (p-values) by IL, TR and RE in v12, v34 and v56
gm.2e.IL <- exp(mean(log(pred.val.2e.IL)))
gm.2e.TR <- exp(mean(log(pred.val.2e.TR)))
gm.2e.RE <- exp(mean(log(pred.val.2e.RE)))

gm.3e.IL <- exp(mean(log(pred.val.3e.IL)))
gm.3e.TR <- exp(mean(log(pred.val.3e.TR)))
gm.3e.RE <- exp(mean(log(pred.val.3e.RE)))

gm.4e.IL <- exp(mean(log(pred.val.4e.IL)))
gm.4e.TR <- exp(mean(log(pred.val.4e.TR)))
gm.4e.RE <- exp(mean(log(pred.val.4e.RE)))

gm.5e.IL <- exp(mean(log(pred.val.5e.IL)))
gm.5e.TR <- exp(mean(log(pred.val.5e.TR)))
gm.5e.RE <- exp(mean(log(pred.val.5e.RE)))

# Now calculate geometric mean of the prediction values (p-values) by v1v2, v3v4 and v5v6 in different locations

gm.2e.v12 <- exp(mean(log(c(pred.val.2e.IL[1],pred.val.2e.TR[1],pred.val.2e.RE[1]))))
gm.2e.v34 <- exp(mean(log(c(pred.val.2e.IL[2],pred.val.2e.TR[2],pred.val.2e.RE[2]))))
gm.2e.v56 <- exp(mean(log(c(pred.val.2e.IL[3],pred.val.2e.TR[3],pred.val.2e.RE[3]))))

gm.3e.v12 <- exp(mean(log(c(pred.val.3e.IL[1],pred.val.3e.TR[1],pred.val.3e.RE[1]))))
gm.3e.v34 <- exp(mean(log(c(pred.val.3e.IL[2],pred.val.3e.TR[2],pred.val.3e.RE[2]))))
gm.3e.v56 <- exp(mean(log(c(pred.val.3e.IL[3],pred.val.3e.TR[3],pred.val.3e.RE[3]))))

gm.4e.v12 <- exp(mean(log(c(pred.val.4e.IL[1],pred.val.4e.TR[1],pred.val.4e.RE[1]))))
gm.4e.v34 <- exp(mean(log(c(pred.val.4e.IL[2],pred.val.4e.TR[2],pred.val.4e.RE[2]))))
gm.4e.v56 <- exp(mean(log(c(pred.val.4e.IL[3],pred.val.4e.TR[3],pred.val.4e.RE[3]))))

gm.5e.v12 <- exp(mean(log(c(pred.val.5e.IL[1],pred.val.5e.TR[1],pred.val.5e.RE[1]))))
gm.5e.v34 <- exp(mean(log(c(pred.val.5e.IL[2],pred.val.5e.TR[2],pred.val.5e.RE[2]))))
gm.5e.v56 <- exp(mean(log(c(pred.val.5e.IL[3],pred.val.5e.TR[3],pred.val.5e.RE[3]))))

# Now calculate the geometric mean for the prediction values of each amplicon all over enterotypes
gm.v12 <- exp(mean(log(c(gm.2e.v12,gm.3e.v12,gm.4e.v12,gm.5e.v12))))
gm.v34 <- exp(mean(log(c(gm.2e.v34,gm.3e.v34,gm.4e.v34,gm.5e.v34))))
gm.v56 <- exp(mean(log(c(gm.2e.v56,gm.3e.v56,gm.4e.v56,gm.5e.v56))))
gm.loc <- c(gm.v12,gm.v34,gm.v56)
gm.loc <- signif(gm.loc, digits = 3)
names(gm.loc) <- c('v1v2','v3v4','v5v6')


# prepare a results' table: rows = IL,TR,RE and geom. mean,
# columns = v1v2, v3v4, v5v6 and geom. mean
result.2e.IL <- c(pred.val.2e.IL,gm.2e.IL)
result.2e.TR <- c(pred.val.2e.TR,gm.2e.TR)
result.2e.RE <- c(pred.val.2e.RE,gm.2e.RE)
result.2e <- rbind(result.2e.IL,result.2e.TR,result.2e.RE,c(gm.2e.v12,gm.2e.v34,gm.2e.v56,NaN))
colnames(result.2e) <- c('v1v2','v3v4','v5v6','geom.mean.locat')
rownames(result.2e) <- c('IL','TR','RE','geom.mean.locat')

result.3e.IL <- c(pred.val.3e.IL,gm.3e.IL)
result.3e.TR <- c(pred.val.3e.TR,gm.3e.TR)
result.3e.RE <- c(pred.val.3e.RE,gm.3e.RE)
result.3e <- rbind(result.3e.IL,result.3e.TR,result.3e.RE,c(gm.3e.v12,gm.3e.v34,gm.3e.v56,NaN))
colnames(result.3e) <- c('v1v2','v3v4','v5v6','geom.mean.locat')
rownames(result.3e) <- c('IL','TR','RE','geom.mean.locat')

result.4e.IL <- c(pred.val.4e.IL,gm.4e.IL)
result.4e.TR <- c(pred.val.4e.TR,gm.4e.TR)
result.4e.RE <- c(pred.val.4e.RE,gm.4e.RE)
result.4e <- rbind(result.4e.IL,result.4e.TR,result.4e.RE,c(gm.4e.v12,gm.4e.v34,gm.4e.v56,NaN))
colnames(result.4e) <- c('v1v2','v3v4','v5v6','geom.mean.locat')
rownames(result.4e) <- c('IL','TR','RE','geom.mean.locat')

result.5e.IL <- c(pred.val.5e.IL,gm.5e.IL)
result.5e.TR <- c(pred.val.5e.TR,gm.5e.TR)
result.5e.RE <- c(pred.val.5e.RE,gm.5e.RE)
result.5e <- rbind(result.5e.IL,result.5e.TR,result.5e.RE,c(gm.5e.v12,gm.5e.v34,gm.5e.v56,NaN))
colnames(result.5e) <- c('v1v2','v3v4','v5v6','geom.mean.locat')
rownames(result.5e) <- c('IL','TR','RE','geom.mean.locat')

# round the numbers
result.2e <- signif(result.2e,digits = 3)
result.3e <- signif(result.3e,digits = 3)
result.4e <- signif(result.4e,digits = 3)
result.5e <- signif(result.5e,digits = 3)

# 3. write the results of 1. and 2. to files
setwd('/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Compare_enterotypings/comp_amplic_and_locat_pred_value')
# of step 2.
write.table(result.2e,'2_enterotypes_pval.csv',sep='\t')
write.table(result.3e,'3_enterotypes_pval.csv',sep='\t')
write.table(result.4e,'4_enterotypes_pval.csv',sep='\t')
write.table(result.5e,'5_enterotypes_pval.csv',sep='\t')
write.table(as.data.frame(gm.loc),'all_enterotypes_pval.csv',sep='\t')
# of step 1.
write.table(result.sep.perc.2e,'2_enterotypes_percent.csv',sep='\t')
write.table(result.sep.perc.3e,'3_enterotypes_percent.csv',sep='\t')
write.table(result.sep.perc.4e,'4_enterotypes_percent.csv',sep='\t')
write.table(result.sep.perc.5e,'5_enterotypes_percent.csv',sep='\t')

# Finally, compare enterotype predictions by loc together vs loc separately
pred.val.dolya.all_loc_together <- c(pred.val.dolya.all_loc_together.2e[4],pred.val.dolya.all_loc_together.3e[4],pred.val.dolya.all_loc_together.4e[4],pred.val.dolya.all_loc_together.5e[4])
names(pred.val.dolya.all_loc_together) <- c('2_ent','3_ent','4_ent','5_ent')
pred.val.dolya.all_loc_sep.2e <- sum(pred.percent.2e.IL$ident.number[4],pred.percent.2e.TR$ident.number[4],pred.percent.2e.RE$ident.number[4])/sum(pred.percent.2e.IL$total.sample.number,pred.percent.2e.TR$total.sample.number,pred.percent.2e.RE$total.sample.number)
pred.val.dolya.all_loc_sep.3e <- sum(pred.percent.3e.IL$ident.number[4],pred.percent.3e.TR$ident.number[4],pred.percent.3e.RE$ident.number[4])/sum(pred.percent.3e.IL$total.sample.number,pred.percent.3e.TR$total.sample.number,pred.percent.3e.RE$total.sample.number)
pred.val.dolya.all_loc_sep.4e <- sum(pred.percent.4e.IL$ident.number[4],pred.percent.4e.TR$ident.number[4],pred.percent.4e.RE$ident.number[4])/sum(pred.percent.4e.IL$total.sample.number,pred.percent.4e.TR$total.sample.number,pred.percent.4e.RE$total.sample.number)
pred.val.dolya.all_loc_sep.5e <- sum(pred.percent.5e.IL$ident.number[4],pred.percent.5e.TR$ident.number[4],pred.percent.5e.RE$ident.number[4])/sum(pred.percent.5e.IL$total.sample.number,pred.percent.5e.TR$total.sample.number,pred.percent.5e.RE$total.sample.number)
pred.val.dolya.all_loc_sep <- c(pred.val.dolya.all_loc_sep.2e,pred.val.dolya.all_loc_sep.3e,pred.val.dolya.all_loc_sep.4e,pred.val.dolya.all_loc_sep.5e)
names(pred.val.dolya.all_loc_sep) <- c('2_ent','3_ent','4_ent','5_ent')
print('Now, comparison of prediction values for each amplicon: all locations 
      enterotyped together vs all locations enterotyped separately')
signif(pred.val.dolya.all_loc_together,digits = 4)
signif(pred.val.dolya.all_loc_sep,digits = 4)
chisq.test(pred.val.dolya.all_loc_together,pred.val.dolya.all_loc_sep)