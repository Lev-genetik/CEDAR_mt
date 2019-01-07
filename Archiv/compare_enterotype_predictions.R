#the function compares 3 data frames "data.sample.enterotype": df1 (V12), df2(V34) and df3(V56) (columns: sample and enterotype) 
# and yields number of samples that share the same enterotype 

require('pracma')
CompEntPred <- function(df1, df2, df3){
  df1_mod <- df1
  #add a column to df1_mod$data.sample.enterotype
  df1_mod$data.sample.enterotype$driver <- NA
  # ent1 change to 1 etc.
  df1_mod$data.sample.enterotype$enterotype_number <- gsub ('ent1','1', df1_mod$data.sample.enterotype$enterotype_number) 
  df1_mod$data.sample.enterotype$enterotype_number <- gsub ('ent2','2', df1_mod$data.sample.enterotype$enterotype_number) 
  df1_mod$data.sample.enterotype$enterotype_number <- gsub ('ent3','3', df1_mod$data.sample.enterotype$enterotype_number) 
  df1_mod$data.sample.enterotype$enterotype_number <- as.integer(df1_mod$data.sample.enterotype$enterotype_number)
  #full the 'driver' column of the table
    for (n in 1:nrow(df1_mod$data.sample.enterotype)){
     df1_mod$data.sample.enterotype$driver[n]  <- df1_mod$enterotype.info$driver[df1_mod$data.sample.enterotype$enterotype_number[n]]
    }

#the same operation for df2 and df3    
  df2_mod <- df2
  #add a column to df2_mod$data.sample.enterotype
  df2_mod$data.sample.enterotype$driver <- NA
  # ent1 change to 1 etc.
  df2_mod$data.sample.enterotype$enterotype_number <- gsub ('ent1','1', df2_mod$data.sample.enterotype$enterotype_number) 
  df2_mod$data.sample.enterotype$enterotype_number <- gsub ('ent2','2', df2_mod$data.sample.enterotype$enterotype_number) 
  df2_mod$data.sample.enterotype$enterotype_number <- gsub ('ent3','3', df2_mod$data.sample.enterotype$enterotype_number) 
  df2_mod$data.sample.enterotype$enterotype_number <- as.integer(df2_mod$data.sample.enterotype$enterotype_number)
  #full the 'driver' column of the table
  for (n in 1:nrow(df2_mod$data.sample.enterotype)){
    df2_mod$data.sample.enterotype$driver[n]  <- df2_mod$enterotype.info$driver[df2_mod$data.sample.enterotype$enterotype_number[n]]
  }
  
  df3_mod <- df3
  #add a column to df3_mod$data.sample.enterotype
  df3_mod$data.sample.enterotype$driver <- NA
  # ent1 change to 1 etc.
  df3_mod$data.sample.enterotype$enterotype_number <- gsub ('ent1','1', df3_mod$data.sample.enterotype$enterotype_number) 
  df3_mod$data.sample.enterotype$enterotype_number <- gsub ('ent2','2', df3_mod$data.sample.enterotype$enterotype_number) 
  df3_mod$data.sample.enterotype$enterotype_number <- gsub ('ent3','3', df3_mod$data.sample.enterotype$enterotype_number) 
  df3_mod$data.sample.enterotype$enterotype_number <- as.integer(df3_mod$data.sample.enterotype$enterotype_number)
  #full the 'driver' column of the table
  for (n in 1:nrow(df3_mod$data.sample.enterotype)){
    df3_mod$data.sample.enterotype$driver[n]  <- df3_mod$enterotype.info$driver[df3_mod$data.sample.enterotype$enterotype_number[n]]
  }
  
  #now compare the enterotypes. Based on df1, go row by row
 
  #get rid of .V12, .V34 and .V56
  df1_mod$data.sample.enterotype$sample <- gsub('.V12','',df1_mod$data.sample.enterotype$sample)
  df2_mod$data.sample.enterotype$sample <- gsub('.V34','',df2_mod$data.sample.enterotype$sample)
  df3_mod$data.sample.enterotype$sample <- gsub('.V56','',df3_mod$data.sample.enterotype$sample)
  
#defining the result vector
  n_row_df1 <- as.integer(nrow(df1$data.sample.enterotype))
  result <- integer(n_row_df1)

#comparison of enterotypes per se  
    for (n in 1:n_row_df1){
     current_sample <- df1_mod$data.sample.enterotype$sample[n]
     # print(current_sample)
     entV1V2 <- df1_mod$data.sample.enterotype$driver[n]
     if((current_sample %in% df2_mod$data.sample.enterotype$sample)&(current_sample %in% df3_mod$data.sample.enterotype$sample)){
      entV3V4 <- df2_mod$data.sample.enterotype[df2_mod$data.sample.enterotype$sample==current_sample]$driver
      entV5V6 <- df3_mod$data.sample.enterotype[df3_mod$data.sample.enterotype$sample==current_sample]$driver
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
  jpeg('out/Compare_2enterotypes_full.jpg')
  pie(numbers,labels,init.angle = 0)
  dev.off()
  jpeg('out/Compare_2enterotypes_short.jpg')
  pie(numbers[1:3],labels[1:3],init.angle = 0)
  dev.off()
  summary_final <- rbind(labels,numbers)
  return(summary_final)
   }