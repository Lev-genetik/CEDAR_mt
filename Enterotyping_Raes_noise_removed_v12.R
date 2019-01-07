# v12: pca and bca results added
# v11: sample coverage filter added
# v10 changes: now can remove outliers (e.g. 1% of most distant samples)
# v9 changes: 
# a. now has an additional argument driv.number (is 10 if not indicated) - number of drivers in drivers_n output element
# b. new output $obs_bet_result contains result of bca function - can be used e.g. to get driver percentages
# by GetDrivers(ent.Raes12.3e$obs_bet_result,n=10) = for first 10 drivers etc.
#input arg 1: the file with taxa abundance.txt (e.g. for genus level /media/lev-genetik/980E73270E72FD96
#/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt) 
# WITH DELETED # Constructed from biom file #OTU ID	!!!
# input arg 2 and 3: output directory, and if noise = taxa with <0.01 average abundance, should be removed
# input arg 4: locat. If exists, only the columns that contain that argument are 
# selected. Can be '.IL.', '.TR.' or '.RE.'
# input arg 5: if plots=F, than arg2 can be '' or anything
# input arg remove.outliers: if not removing outliers, must be = 0 (po umolch). Otherwise if = 0.01 - remove 
# 1% of samples with highest distance to closest 50% of samples. Can be 0 to 1 potentially.
# input arg: sample.coverage.filter = minimum number of reads per samples that are included into enterotyping 
# returns a list:
# $sample.enterotype = data.frame sample-enterotype
# $clusters = a matrix with raws corresponding to clusters and columns to i-th driver of the cluster
#and generates PCA/PCoA results for clustering in the output directory if plots=T

#Uncomment next two lines if R packages are already installed
#install.packages("cluster")
#install.packages("clusterSim")
Enterotyping_Raes <- function (input_file, output_dir, clust.number = 4, noiserem=T, locat = 'all', remove.outliers = 0, 
                               sample.coverage.filter=0, driv.number = 10, plots=F) {
library(cluster)
library(clusterSim)
library(scatterplot3d)
library(knb16slib)
require(ade4)
library(dplyr)
drops.new = NULL

#Download the example data and set the working directory
data=read.table(input_file, header=T, row.names=1, dec=".", sep="\t")
# select specific location - that can be .IL., .TR. or .RE.
if (locat!='all'){
  data <- dplyr::select(data,contains(locat))
}
# data=data[-1,]

# get rid of samples with microbiota coverage < sample.coverage.filter
if (sample.coverage.filter!=0){
  colSummi <- colSums(data)
  data <- data[,colSummi>=sample.coverage.filter]
}
# # data=noise.removal(data, percent=0.01)
noise.removal <- function(dataframe, percent=0.01, top=NULL){
  dataframe->Matrix
  bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent
  Matrix_1 <- Matrix[bigones,]
  # print(percent)
  return(Matrix_1)
}

if(noiserem){data=noise.removal(data, percent=0.01)}

# in order to avoid error if driv.number>number of bacteria
driv.number <- min(driv.number,dim(data)[1])

dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
  
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}

data.dist=dist.JSD(data)

# New in v10: remove outliers before enterotyping
if (remove.outliers!=0){
  library(fitdistrplus)
  dist.m <- as.matrix(data.dist)
  n <- round(nrow(dist.m)/2)
  # the median distance from a sample to 50% of closest
  meds.50 <- apply(dist.m, 1, function(x) {
    x.sort <- sort(x, decreasing = F)[-1]
    median(x.sort[1:n])
  })
  # hist(meds.50, 100, density = T)
  f1 <- fitdist(meds.50,"norm")
  # plotdist(meds.50,"norm",para=list(mean=f1$estimate[1],sd=f1$estimate[2])) #to plot
  p.vals <- pnorm(meds.50, mean = f1$estimate[1], sd = f1$estimate[2], lower.tail = F, log.p = F)
  drops.new <- names(p.vals)[p.vals < remove.outliers]
  # and now delete the outliers
  dist.m.new <- dist.m[!rownames(dist.m) %in% drops.new,!colnames(dist.m) %in% drops.new]
# here we get the new matrix data.dist with outliers removed!
  data.dist <- dist.m.new
  data.dist <- as.dist(data.dist)
  attr(data.dist, "method") <- "dist"
  
  # and remove outliers from original data
  data <- data[,!colnames(data) %in% drops.new]
}

pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}
data.cluster=pam.clustering(data.dist, k=clust.number)

require(clusterSim)
# nclusters = index.G1(t(data), data.cluster, d = data.dist, centrotypes = "medoids")

nclusters=NULL
nsilhouettes=NULL

for (k in 1:20) { 
  if (k==1) {
    nclusters[k]=NA 
    nsilhouettes[k]=NA 
  } else {
    data.cluster_temp=pam.clustering(data.dist, k)
    nclusters[k]=index.G1(t(data),data.cluster_temp,  d = data.dist,
                          centrotypes = "medoids")
    nsilhouettes[k]=mean(silhouette(data.cluster_temp, data.dist)[,3])
    
  }
}
obs.pca.30=dudi.pca(data.frame(t(data)), scannf=F, nf=30) #new in v12
obs.pca=dudi.pca(data.frame(t(data)), scannf=F, nf=10)
obs.bet=bca(obs.pca, fac=as.factor(data.cluster), scannf=F, nf=k-1) 

#making a matrix with raws corresponding to clusters and columns to i-th driver of the cluster
k <- clust.number
clusters <- matrix(0,k,driv.number)
for (i in 1:k){
  for (j in 1:driv.number){
    clusters[i,j] <- GetDrivers(obs.bet[i], n=as.numeric(j))[[2]][[1]][[1]][[as.numeric(j)]]
  }
}

# here in case we have 2 enterotypes with same first drivers, we add second driver after '->'
drivers_unique <- as.vector(clusters[,1])
if(anyDuplicated(drivers_unique)){
  i=2
  while((anyDuplicated(drivers_unique))&(i<=ncol(clusters))){
    for (t in 1:length(drivers_unique)){
      drivers_unique[t] <- paste(drivers_unique[t],clusters[t,i],sep='->')
    }
  }
}

# make a data frame sample-enterotype
zeros <- rep(0, length(colnames(data))*2)
sample.enterotype <- array(zeros, c(length(colnames(data)), 2))
for (i in 1:length(colnames(data))){
  sample.enterotype[i,1] <- colnames(data[i])
  sample.enterotype[i,2] <- drivers_unique[data.cluster[i]]

}
if (plots==T) {
  # print('plot=T')
  setwd(output_dir)
  filename_CH = paste0(word(input_file,start = -2,sep='/'),'.','locat_',locat,'.','CH_pseudo_F.jpg')
  jpeg(filename = filename_CH)
  plot(nclusters, type="h", xlab="Number of clusters", ylab="CH index",main="Calinski-Harabasz pseudo F-statistic", lwd=2, cex.lab=1)
  dev.off()
  filename_Sil = paste0(word(input_file,start = -2,sep='/'),'.','locat_',locat,'.','mean_silhouette.jpg')
  jpeg(filename = filename_Sil)
  plot(nsilhouettes, type="h", xlab="Number of clusters", ylab="Mean silhouette",main="Mean silhouette vs number of clusters",lwd=2, cex.lab=1)
  dev.off()
  
  obs.silhouette=mean(silhouette(data.cluster, data.dist)[,3])
# cat(obs.silhouette) #0.1899451

dir.create(file.path(output_dir, 'BCA_and_PCOA'), showWarnings = FALSE)
setwd(file.path(output_dir, 'BCA_and_PCOA'))
## plot 1
if(dim(obs.bet$ls)[2]>1){
  filename_BCA = paste0(word(input_file,start = -2,sep='/'),'.',clust.number,'e','.','locat_',locat,'.','Between_class_anal_2_eigenvectors_PCA.jpg')
  jpeg(filename = filename_BCA, width=2048, height=2048)
  s.class(obs.bet$ls, fac=as.factor(data.cluster), col=c('red','orange','green','blue'), cpoint = 0.5, grid=F, sub="Between-class analysis for the enterotypes")
  dev.off()
}
# else{
#   print('Note: bca result$li contains only 1 axis coordinates! No BCA for 2 eigen vectors')
# }
#plot 2
obs.pcoa.30=dudi.pco(data.dist, scannf=F, nf=30) # this is new in v12 - done for reporting full PCoA
obs.pcoa=dudi.pco(data.dist, scannf=F, nf=3) 
filename_BC_PCOA = paste0(word(input_file,start = -2,sep='/'),'.',clust.number,'e','.','locat_',locat,'.','Between_class_PCOA.jpg')
jpeg(filename = filename_BC_PCOA, width=2048, height=2048)
s.class(obs.pcoa$li, fac=as.factor(data.cluster), col=c('red','orange','green','blue'), grid=F,sub="Principal coordiante analysis")
dev.off()
filename_BC_PCOA_3D = paste0(word(input_file,start = -2,sep='/'),'.',clust.number,'e','.','locat_',locat,'.','Between_class_PCOA_3D.jpg')
jpeg(filename = filename_BC_PCOA_3D, width=640, height=640)
scatterplot3d(obs.pcoa$li, color=c('red','orange','green','blue')[as.factor(data.cluster)])
dev.off()
}

colnames(sample.enterotype) <- c('sample','enterotype')
sample.enterotype <- as.data.frame(sample.enterotype)




result <- list(clusters,sample.enterotype,drivers_unique,obs.bet,obs.pca.30,obs.pcoa.30,drops.new)
names(result) <- c('drivers_n','sample.enterotype','drivers_unique','obs_bet_result','obs_pca_result',
                   'obs_pcoa_result','outliers')
return(result)
}