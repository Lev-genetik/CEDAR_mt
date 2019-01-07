
# Read a file of taxa abundance and output a data frame of taxa abundance with:
# 0. location of interest selected if specified
# 1. taxa with abundance < percent_threshold% removed
# 2. Outliers with average distance to closest 50% neighbours >remove.outliers p-value removed
readQiimeSmart <- function(input_file, locat = 'all',percent_threshold = 0.01,remove.outliers = 0,
                           sample.cov.filter=0,to.scale=F){
  library(dplyr)
  data=read.table(input_file, header=T, row.names=1, dec=".", sep="\t")
  # choose the location if interest
  if (locat!='all'){
    data <- dplyr::select(data,contains(locat))
  }
  # filter samples by coverage
  data <- data[,colSums(data)>=sample.cov.filter]
  # remove noise
  noise.removal <- function(dataframe, top=NULL){
    dataframe->Matrix
    bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent_threshold
    Matrix_1 <- Matrix[bigones,]
    return(Matrix_1)
  }
  data=noise.removal(data)
  # remove outliers
  # first, calculate the JSD metric
  if(remove.outliers!=0){
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
    # this to remove outliers per se
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
    # and finally remove outliers from original data
    data <- data[,!colnames(data) %in% drops.new]
  }
  if(to.scale){
    data_col_sum <- colSums(data)
    data_norm <- data
    for (i in 1:dim(data)[2]){
      data_norm[,i] <- data[,i]/data_col_sum[i]
    }
    data <- data_norm
  }
  
  return(data)
}


to.scale <- function(data){
  data_col_sum <- colSums(data)
  data_norm <- data
  for (i in 1:dim(data)[2]){
    data_norm[,i] <- data[,i]/data_col_sum[i]
  }
  return(data_norm)
}
