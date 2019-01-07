#the function compares 3 data frames : ent1, ent2, ent3 (columns: [,1] = sample and [,2] = enterotype) 
# and yields number of samples that share the same enterotype
#program has bugs
# does not take into account all possible cluster correspondences, used for analyzing Skolkovo clustering data
CompEntPred1 <- function(ent1,ent2,ent3,out_dir){
  #get rid of .V12, .V34 and .V56
  ent1 <- gsub('.V12','',ent1)
  ent2 <- gsub('.V34','',ent2)
  ent3 <- gsub('.V56','',ent3)
  
#defining the result vector
  n_row_ent1 <- nrow(ent1)
  result <- integer(n_row_ent1)

#comparison of enterotypes per se  
    for (n in 1:n_row_ent1){
     current_sample <- ent1[n,1]
     # print(current_sample)
     entV1V2 <- ent1[n,2]
     if((current_sample %in% ent2)&(current_sample %in% ent3)){
      entV3V4 <- ent2[ent2$sample==current_sample]$enterotype
      entV5V6 <- ent3[ent3$sample==current_sample]$enterotype
      if (strcmp(entV1V2,entV3V4)&strcmp(entV1V2,entV5V6)){
        result[n] <- entV1V2
        print(paste0('concordant, n= ',n))
      }else{result[n] <- 'discordant'
       print(paste0('discordant, n= ',n))
      }
    } else{print('some data not available for this sample')}
  
    }
  summary <- table(result)
  Bacteroides <- as.integer(summary[names(summary)=='k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides'])
  Prevotella <- as.integer(summary[names(summary)=='k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotella'])
  Discordant <- as.integer(summary[names(summary)=='discordant'])
  Non_ful_gen <- as.integer(summary[names(summary)=='0'])
  labels <- c('Bacteroides','Prevotella', 'Discordant', 'Not genotyped entirely')
  numbers <- c(Bacteroides, Prevotella, Discordant, Non_ful_gen)
  setwd(out_dir)
  jpeg('Compare_enterotypes_full.jpg')
  pie(numbers,labels,init.angle = 0)
  dev.off()
  jpeg('Compare_enterotypes_short.jpg')
  pie(numbers[1:3],labels[1:3],init.angle = 0)
  dev.off()
  summary_final <- rbind(labels,numbers)
  return(summary_final)
}
