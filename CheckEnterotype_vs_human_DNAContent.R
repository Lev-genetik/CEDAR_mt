# AIM: calculate is enterotype depends on proportion of human DNA
# INPUT: can be run as CheckHumanDNA() with the default input files
#  DNA.cont.file = data frame with sample in column 1 and DNA content in column 11(V12),12(V34) and 13(V56), respectively
#  file1, file2 and file3 = files with data frames sample-enterotype obtained by /home/lev-genetik/Desktop/Projects/liege/src/Lev/Script_Get_Enterotypes_and_compare_enterotyping_V12_V34_V56_pairwise_v5.R
#  assumed we work with 3 enterotypes
# OUTPUT: a matrix with each cell corresponding to a p-value of the hypothesis that enterotypes
#  do depend on human DNA contamination (p.value <0.01 => depend!). In our case, V34 enterotypes do depend
#  on human DNA content. Assumed that ent1 is Bacteroides, ent2 = Prevotella, ent3 = other.

CheckHumanDNA <- function(DNA.cont.file='/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/Human_dna_stats_Zhenya_2018_09_12.csv',file1='/media/lev-genetik/980E73270E72FD96/Liege/Results/Clustering/Sample_cluster/2018_09_06_enterotyping/V1V2_3e_sample-enterotype.csv',file2='/media/lev-genetik/980E73270E72FD96/Liege/Results/Clustering/Sample_cluster/2018_09_06_enterotyping/V3V4_3e_sample-enterotype.csv',file3='/media/lev-genetik/980E73270E72FD96/Liege/Results/Clustering/Sample_cluster/2018_09_06_enterotyping/V5V6_3e_sample-enterotype.csv',plots=FALSE, plotdir = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Clustering/Sample_cluster'){
library(stringr)
hum.DNA.table.34 <- read.table(DNA.cont.file,header = T)
ent.V12.3e <- read.table(file1, header = T)
ent.V34.3e <- read.table(file2, header = T)
ent.V56.3e <- read.table(file3, header = T)


# delete replicate samples (they have .1,.2 or .3 at the end)
# j.table <- j.table[-grep('\\.[0-9]',j.table[,1],value = F),]

# substitute .V34 by .V12 and .V56
hum.DNA.table.12 <- hum.DNA.table.34
hum.DNA.table.56 <- hum.DNA.table.34
hum.DNA.table.12[,1] <- str_replace(hum.DNA.table.34[,1],'\\.V34','\\.V12')
hum.DNA.table.56[,1] <- str_replace(hum.DNA.table.34[,1],'\\.V34','\\.V56')

# construct tables for V12, V34 and V56: sample - %human_DNA - enterotype
table.V12 <- ent.V12.3e[,]
match12 <- match(table.V12[,1],hum.DNA.table.12[,1])
# human DNA content for samples in the order as in table.V12[,1]
human.cont.12 <- hum.DNA.table.12[match12,11]
# now add the desired 'human DNA percent' column to the table. Human DNA content really corresponds to the sample - I checked it for 2 samples
table.V12 <- cbind(table.V12,human.DNA.cont=human.cont.12)

table.V34 <- ent.V34.3e[,]
match34 <- match(table.V34[,1],hum.DNA.table.34[,1])
# human DNA content for samples in the order as in table.V34[,1]
human.cont.34 <- hum.DNA.table.34[match34,12]
# now add the desired 'human DNA percent' column to the table. Human DNA content really corresponds to the sample - I checked it for 2 samples
table.V34 <- cbind(table.V34,human.DNA.cont=human.cont.34)



table.V56 <- ent.V56.3e[,]
match56 <- match(table.V56[,1],hum.DNA.table.56[,1])
# human DNA content for samples in the order as in table.V12[,1]
human.cont.56 <- hum.DNA.table.56[match56,13]
# now add the desired 'human DNA percent' column to the table. Human DNA content really corresponds to the sample - I checked it for 2 samples
table.V56 <- cbind(table.V56,human.DNA.cont=human.cont.56)

# tables are okay! checked all the 3 tables table.V12, table.V34, table.V56

# Now get the p-values of the hypothesis that human DNA percent is different between enterotypes.
# levels(table.V12[,2])[2] e.g. is Prevotella.
# wilcox.test = Mann-Whitney U test:)
p_val12.12 <- round((wilcox.test(table.V12[table.V12[,2] == levels(table.V12[,2])[1],3],table.V12[table.V12[,2] == levels(table.V12[,2])[2],3]))$p.value,3)
p_val12.13 <- round((wilcox.test(table.V12[table.V12[,2] == levels(table.V12[,2])[1],3],table.V12[table.V12[,2] == levels(table.V12[,2])[3],3]))$p.value,3)
p_val12.23 <- round((wilcox.test(table.V12[table.V12[,2] == levels(table.V12[,2])[2],3],table.V12[table.V12[,2] == levels(table.V12[,2])[3],3]))$p.value,3)

p_val34.12 <- wilcox.test(table.V34[table.V34[,2] == levels(table.V34[,2])[1],3],table.V34[table.V34[,2] == levels(table.V34[,2])[2],3])$p.value
p_val34.13 <- wilcox.test(table.V34[table.V34[,2] == levels(table.V34[,2])[1],3],table.V34[table.V34[,2] == levels(table.V34[,2])[3],3])$p.value
p_val34.23 <- wilcox.test(table.V34[table.V34[,2] == levels(table.V34[,2])[2],3],table.V34[table.V34[,2] == levels(table.V34[,2])[3],3])$p.value
if(p_val34.12>=0.001){
  p_val34.12 <- round(p_val34.12,3)
}
if(p_val34.13>=0.001){
  p_val34.13 <- round(p_val34.13,3)
}
if(p_val34.23>=0.001){
  p_val34.23 <- round(p_val34.23,3)
}

p_val56.12 <- round((wilcox.test(table.V56[table.V56[,2] == levels(table.V56[,2])[1],3],table.V56[table.V56[,2] == levels(table.V56[,2])[2],3]))$p.value,3)
p_val56.13 <- round((wilcox.test(table.V56[table.V56[,2] == levels(table.V56[,2])[1],3],table.V56[table.V56[,2] == levels(table.V56[,2])[3],3]))$p.value,3)
p_val56.23 <- round((wilcox.test(table.V56[table.V56[,2] == levels(table.V56[,2])[2],3],table.V56[table.V56[,2] == levels(table.V56[,2])[3],3]))$p.value,3)

result <- matrix(c(p_val12.12,p_val12.13,p_val12.23,p_val34.12,p_val34.13,p_val34.23,p_val56.12,p_val56.13,p_val56.23),nrow = 3,ncol=3)
# Here checked that ent1 = Bacteroides and ent2 = Prevotella

if ((length(grep('Bacteroides',levels(table.V12[,2])[1]))<1)|(length(grep('Prevotella',levels(table.V12[,2])[2]))<1)){
  warning('ent1 or ent2 enterotype for V12 has drivers different from Bacteroides and Prevotella respectively')
}
if ((length(grep('Bacteroides',levels(table.V34[,2])[1]))<1)|(length(grep('Prevotella',levels(table.V34[,2])[2]))<1)){
  warning('ent1 or ent2 enterotype for V34 has drivers different from Bacteroides and Prevotella respectively')
}
if ((length(grep('Bacteroides',levels(table.V56[,2])[1]))<1)|(length(grep('Prevotella',levels(table.V56[,2])[2]))<1)){
  warning('ent1 or ent2 enterotype for V56 has drivers different from Bacteroides and Prevotella respectively')
}
  
rownames(result) <- c('ent1(B)-2(P)','ent1(B)-3','ent2(P)-3')
colnames(result) <- c('V1V2','V3V4','V5V6')

if(plots==T){
  setwd(plotdir)
  jpeg(filename = 'V1V2_human_DNA_percent.jpg')
  plotNormalDensity(hum.DNA.table.12[match12,11])
  dev.off()
  jpeg(filename = 'V3V4_human_DNA_percent.jpg')
  plotNormalDensity(hum.DNA.table.34[match34,12])
  dev.off()
  jpeg(filename = 'V5V6_human_DNA_percent.jpg')
  plotNormalDensity(hum.DNA.table.56[match56,13])
  dev.off()
}

return(result)
}