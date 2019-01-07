#input arg 1: the file with taxa abundance.txt (e.g. for genus level /media/lev-genetik/980E73270E72FD96
#/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt) 
# WITH DELETED # Constructed from biom file #OTU ID	!!!
#input arg 2 and 3: output directory, and if noise = taxa with <0.01 average abundance, should be removed
#input arg 4: if plots=F, than arg2 can be '' or anything
# returns a data.frame sample-enterotype
#and generates PCA/PCoA results for clustering in the output directory

#Uncomment next two lines if R packages are already installed
#install.packages("cluster")
#install.packages("clusterSim")
Enterotyping_Raes <- function (input_file, output_dir, noiserem=T, plots=F) {
library(cluster)
library(clusterSim)
library(scatterplot3d)
library(knb16slib)
require(ade4)

#Download the example data and set the working directory
data=read.table(input_file, header=T, row.names=1, dec=".", sep="\t")
# data=data[-1,]

# # data=noise.removal(data, percent=0.01)
noise.removal <- function(dataframe, percent=0.01, top=NULL){
  dataframe->Matrix
  bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent
  Matrix_1 <- Matrix[bigones,]
  # print(percent)
  return(Matrix_1)
}

if(noiserem){data=noise.removal(data, percent=0.01)}

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

pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}

data.cluster=pam.clustering(data.dist, k=4)

require(clusterSim)
nclusters = index.G1(t(data), data.cluster, d = data.dist, centrotypes = "medoids")

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

obs.pca=dudi.pca(data.frame(t(data)), scannf=F, nf=10)
obs.bet=bca(obs.pca, fac=as.factor(data.cluster), scannf=F, nf=k-1) 

zeros <- rep(0, length(colnames(data))*2)
sample.enterotype <- array(zeros, c(length(colnames(data)), 2))
for (i in 1:length(colnames(data))){
  sample.enterotype[i,1] <- colnames(data[i])
  if(data.cluster[i]==1){sample.enterotype[i,2] <- GetDrivers(obs.bet[data.cluster[i]], n=1)$drivers$'1'$taxon 
  }else if(data.cluster[i]==2){sample.enterotype[i,2] <- GetDrivers(obs.bet[data.cluster[i]], n=1)$drivers$'2'$taxon 
  }else if(data.cluster[i]==3){sample.enterotype[i,2] <- GetDrivers(obs.bet[data.cluster[i]], n=1)$drivers$'3'$taxon 
  }else if(data.cluster[i]==4){sample.enterotype[i,2] <- GetDrivers(obs.bet[data.cluster[i]], n=1)$drivers$'4'$taxon 
  }
}

if (plots==T) {
print('plot=T')
setwd(output_dir)
jpeg(filename = 'CH_Opt_number_clust_for_L6.jpg')
plot(nclusters, type="h", xlab="Number of clusters", ylab="CH index",main="Calinski-Harabasz pseudo F-statistic", lwd=2, cex.lab=1)
dev.off()
jpeg(filename = 'Silhouettes_Opt_number_clust_for_L6.jpg')
plot(nsilhouettes, type="h", xlab="Number of clusters", ylab="Mean silhouette",main="Mean silhouette vs number of clusters",lwd=2, cex.lab=1)
dev.off()

obs.silhouette=mean(silhouette(data.cluster, data.dist)[,3])
# cat(obs.silhouette) #0.1899451


## plot 1
# dev.new()
#finished - colors
jpeg(filename = 'Between_class_anal_2_eigenvectors_PCA.jpg', width=2048, height=2048)
s.class(obs.bet$ls, fac=as.factor(data.cluster), col=c('red','orange','green','blue'), cpoint = 0.5, grid=F, sub="Between-class analysis for the enterotypes")
dev.off()
#plot 2
obs.pcoa=dudi.pco(data.dist, scannf=F, nf=3)
jpeg(filename = 'Between_class_PCOA.jpg', width=2048, height=2048)
s.class(obs.pcoa$li, fac=as.factor(data.cluster), col=c('red','orange','green','blue'), grid=F,sub="Principal coordiante analysis")
dev.off()
jpeg(filename = 'Between_class_PCOA_3D.jpg', width=640, height=640)
scatterplot3d(bs.pcoa$li, color=c('red','orange','green','blue')[as.factor(data.cluster)])
dev.off()
}

colnames(sample.enterotype) <- c('sample','enterotype')
sample.enterotype <- as.data.frame(sample.enterotype)
return(sample.enterotype)
}