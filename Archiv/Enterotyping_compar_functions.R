#compare enterotype predictions.
#n number of data frames with 2 columns: sample and enterotype
#make a distance matrix between the samples (distance=number of identical enterotypes)
CompEntPred.matrix <- function(pred.list){
#get rid of .V12, .V34 and .V56
  if (is.list(pred.list)==F){print ('Error: argument should be a list')
    stop()
    geterrmessage('Error: argument should be a list')
  }
for (i in 1:length(pred.list)){
  pred.list[[i]]$sample <- gsub('.V12','',pred.list[[i]]$sample)
  pred.list[[i]]$sample <- gsub('.V34','',pred.list[[i]]$sample)
  pred.list[[i]]$sample <- gsub('.V56','',pred.list[[i]]$sample)
}

# fill in the resul matrixmatrix
result <- matrix(0,length(pred.list),length(pred.list))
for (i in 1:length(pred.list)){
  for (j in 1:length(pred.list)){
    #now get concordance percent for 2 data frames
    result(i,j) <- round(concordance(pred.list[[i]],pred.list[[j]]),3)
  }
}  
return(result)  

}
