# functions_transcriptome.R
# Adapted from: code/CEDAR-master/Microbiome_genome_transcr_comp/Microbiome_transcr_functions_v5.R
# Functions for comparing microbiome with transcriptome data:
#   PC.scatter.plot.create(), corr.trPCs.vs.mPCs(), microb.transcr.make.QQplot(),
#   write.pvalue.vectors.to.file(), corrtest.microb.transcript(),
#   draw.boxPlot.enterotypes()
#
# Dependencies (sourced before this file by 00_setup.ipynb):
#   functions_general.R   — merge.by.rownames.flex(), remove.outliers(),
#                           write.table.smart(), correct4inflation()
#   functions_association.R — check.PC.assoc()
#   functions_figures.R   — qqunif.plot()
#

require(ggplot2)

# --- Scatter plot: microbiome PC vs transcriptome PC ---

PC.scatter.plot.create <- function(microbiome.df, transcriptome.df,
                                   microbiome.PC, transcriptome.PC,
                                   remove.outliers = c(F,F),
                                   remove.outliers.IQR = c(1.5,1.5),
                                   tissue, xlab_ = NULL, ylab_ = NULL,
                                   file.out, width = 800, height = 800,
                                   size.axis.text = 40, size.axis.title = 50,
                                   digits_beta = 3){
    if(is.numeric(microbiome.PC)){
        microbiome.PC_full <- paste0('RS',microbiome.PC)
    } else if (sum((colnames(microbiome.df)==microbiome.PC)==1)){
        microbiome.PC_full <- microbiome.PC
    }else{
        stop('microbiome.PC must be numeric')
    }
    if(is.numeric(transcriptome.PC)){
        transcriptome.PC_full <- paste0('PC',transcriptome.PC)
    } else if (sum((colnames(transcriptome.df)==transcriptome.PC)==1)){
        transcriptome.PC_full <- transcriptome.PC
    } else{
        stop('transcriptome.PC must be numeric')
    }

    df.microb.transcr <- merge.by.rownames.flex(microbiome.df,transcriptome.df,
                                           all.x = F, all.y = F)
    if(sum(colnames(microbiome.df)%in%microbiome.PC_full)!=1){
        stop('microbiome PC not found in the microbiome.df')
    }
    if(sum(colnames(transcriptome.df)%in%transcriptome.PC_full)!=1){
        stop('transcriptome PC not found in the transcriptome.df')
    }
    df.microb.transcr.current.PCs <- df.microb.transcr[,c(microbiome.PC_full,
                                                          transcriptome.PC_full)]
    # remove outliers if required
    if(remove.outliers[1]){
        df.microb.transcr.current.PCs[,1] <- remove.outliers(
            df.microb.transcr.current.PCs[,1], remove.outliers.IQR = remove.outliers.IQR[1])
    }
    if(remove.outliers[2]){
        df.microb.transcr.current.PCs[,2] <- remove.outliers(
            df.microb.transcr.current.PCs[,2], remove.outliers.IQR = remove.outliers.IQR[2])
    }
    pval.current.assoc <- anova(lm(df.microb.transcr.current.PCs[,1]~df.microb.transcr.current.PCs[,2]))[1,'Pr(>F)']
    beta_ <- summary(lm(df.microb.transcr.current.PCs[
        ,1]~df.microb.transcr.current.PCs[,2]))[[4]][2,'Estimate']
    print(paste0('Association p-value (lm anova): ',
                 signif(pval.current.assoc, digits = 3))
          )
    print(paste0('Beta: ',
                 signif(beta_, digits = digits_beta))
    )
    # determine axis labels
    if(is.null(xlab_)){
        xlab_ <- paste0('Microbiome PC',microbiome.PC)
    }
    if(is.null(ylab_)){
        ylab_ <- (paste0(tissue,' transcriptome PC ',transcriptome.PC))
    }
    jpeg(file.out,width=width,height=height)
    print(ggplot(df.microb.transcr.current.PCs, aes(x = df.microb.transcr.current.PCs[, 1],y = df.microb.transcr.current.PCs[, 2])) +
        geom_point(aes(size=3)) +
        stat_smooth(method = "lm", col = "red", aes(size=3)) +
        xlab(xlab_)+
        ylab(ylab_)+
        theme(axis.text=element_text(size=size.axis.text),
              axis.title = element_text(size=size.axis.title),
              legend.position = 'none', plot.margin = margin(2, 0, 1, 1, "cm")))
    dev.off()
}

# --- Core: correlate transcriptome PCs vs microbiome PCs ---

# Calculates correlation/association between microbiome PCs/residuals and
# tissue-specific transcriptome PCs/residuals, optionally makes QQ-plots.
#
# microbiome.PCs: data frame (samples x microbiome PCs)
# transcr.PC.dir: directory containing transcriptome PC files
# tissue.PC.file.name.end: suffix of transcriptome PC files
# expression.PCs.names.vector: e.g. c('IL.PCs','TR.PCs','RE.PCs')
# number.of.PCs.transcriptome: number of PCs to use per tissue
# corr.method: 'spearman', 'pearson', 'lm', or 'kruskal'
# remove.outliers: logical or 2-element logical vector
# transcr.exclude.mode: 'last_columns' or 'first_1_column'
# Returns: named vector of p-values
corr.trPCs.vs.mPCs <- function(microbiome.PCs,
                            transcr.PC.dir = NULL,
                            tissue.PC.file.name.end = '_PCs_after_correction_v2.txt',
                            expression.PCs.names.vector, number.of.PCs.transcriptome = NULL,
                            corr.method = 'spearman',
                            remove.outliers = T,
                            remove.outliers.IQR = 1.5,
                            transcr.exclude.mode = 'last_columns',
                            inflation.correct = F,
                            QQplot.file = NULL, width = 800, height = 800,
                            main = 'Microbiome vs transcriptome PCs (corr. signif)',
                            cex.main = 2.5){
    if (is.null(transcr.PC.dir)) {
        stop('transcr.PC.dir must be specified')
    }
    # in case the number of PCs is not assigned, it must be residual analysis mode.
    if(is.null(number.of.PCs.transcriptome)){
        number.of.PCs.transcriptome <- rep(10^5, times = length(expression.PCs.names.vector))
    }
    if(length(number.of.PCs.transcriptome) != length(expression.PCs.names.vector)){
        stop('Lengths of number.of.PCs.transcriptome and expression.PCs.names.vector
             must be same! expression.PCs.names.vector are names of the
             number.of.PCs.transcriptome PCs.')
    }
    names(number.of.PCs.transcriptome) <- expression.PCs.names.vector
    # Read expression PCs (using file.path instead of setwd)
    for(tissue.PC in expression.PCs.names.vector){
        current.tissue.PC.file <- paste0(word(tissue.PC,sep = '\\.', end = 1),
                                         tissue.PC.file.name.end)
        current.tissue.PC <- read.table(
            file.path(transcr.PC.dir, current.tissue.PC.file),
            row.names = 1, header = T)
        assign(tissue.PC,current.tissue.PC)
        # Quality check
        if((sum(rownames(current.tissue.PC) != current.tissue.PC[,1]))!=0){
            stop('tissue PC file reading error')
        }
    }
    # Check correlation between microbiome PCs and transcriptome PCs
    # Make lists: each element = correlation p-values for a specific tissue
    microb.N_PC.tissue.expr_PCs.corr.coef <- list()
    for(i in seq_along(number.of.PCs.transcriptome)){
        current.tissue.PCs <- names(number.of.PCs.transcriptome)[i]
        if(transcr.exclude.mode == 'last_columns'){
            col.to.exclude = list(NULL,c(1,c((number.of.PCs.transcriptome[i]+2):ncol(
                eval(parse(text = current.tissue.PCs))
            ))))
        } else if (transcr.exclude.mode == 'first_1_column'){
            col.to.exclude = list(NULL,c(1))
        } else{
            stop('transcr.exclude.mode error')
        }
        microb.N_PC.tissue.expr_PCs.corr.coef[[i]] <- check.PC.assoc(
            df1 = microbiome.PCs, df2 = eval(parse(text = current.tissue.PCs)),
            method = corr.method,
            col.to.exclude=col.to.exclude,
            pval.adjust = F, 
            remove.outliers = remove.outliers,
            remove.outliers.IQR = remove.outliers.IQR
        )
    }
    # p-values: make a single named vector out of list of data frames
    p_values_expr_microb.list <- lapply(microb.N_PC.tissue.expr_PCs.corr.coef,as.vector)
    p_values_expr_microb <- numeric()
    for (i in 1:length(p_values_expr_microb.list)){
        pvalues.current.tissue <- p_values_expr_microb.list[[i]]
        names.p2 <- gsub('RS','mPC',rep(rownames(microb.N_PC.tissue.expr_PCs.corr.coef[[i]]),
                                        times = length(colnames(microb.N_PC.tissue.expr_PCs.corr.coef[[i]]))))
        names.p3 <- gsub('PC','trPC',rep(colnames(microb.N_PC.tissue.expr_PCs.corr.coef[[i]]), each =
                                             length(rownames(microb.N_PC.tissue.expr_PCs.corr.coef[[i]]))))
        names.p1 <- word(expression.PCs.names.vector[i], end = 1, sep = '\\.')
        names(pvalues.current.tissue) <- paste(names.p1,names.p2,names.p3, sep = '.')
        p_values_expr_microb <- c(p_values_expr_microb,pvalues.current.tissue)
    }
    # correct p-values for inflation if necessary
    if(inflation.correct){
        p_values_expr_microb <- correct4inflation(p_values_expr_microb)
    }
    if(!is.null(QQplot.file)){
        require(lattice)
        jpeg(QQplot.file,width=width,height=height)
        print(qqunif.plot(p_values_expr_microb,conf.col="lightgray", conf.alpha=.05,
                    main=main,
                    par.strip.text = list(cex = 4, mgp = c(6,6,6,6)),
                    scales=list(tck=c(1,0), x=list(cex=2), y=list(cex=2)),
                    cex.main = cex.main, cex.x = 2, cex.y = 2
        ))
        dev.off()
    }
    
    return(p_values_expr_microb)
}



# Optimized version for permutations: Accepts pre-loaded list of data frames
corr.trPCs.vs.mPCs.fast <- function(microbiome.PCs,
                                   transcr.data.list, # Named list of data frames
                                   number.of.PCs.transcriptome = NULL,
                                   corr.method = 'lm',
                                   remove.outliers = T,
                                   remove.outliers.IQR = 1.5,
                                   transcr.exclude.mode = 'last_columns',
                                   inflation.correct = F,
                                   QQplot.file = NULL, width = 800, height = 800,
                                   main = 'Microbiome vs transcriptome PCs (corr. signif)',
                                   cex.main = 2.5) {
    
    # Pre-checks
    expr_names <- names(transcr.data.list)
    if(is.null(number.of.PCs.transcriptome)){
        number.of.PCs.transcriptome <- rep(10^5, times = length(expr_names))
    }
    names(number.of.PCs.transcriptome) <- expr_names
    
    microb.N_PC.tissue.expr_PCs.corr.coef <- list()
    
    for(i in seq_along(expr_names)){
        tissue_key <- expr_names[i]
        current_df <- transcr.data.list[[tissue_key]]
        
        # Determine columns to exclude without eval(parse)
        if(transcr.exclude.mode == 'last_columns'){
            # Keep up to number.of.PCs.transcriptome (skipping the first column which is IDs)
            n_cols <- ncol(current_df)
            limit <- number.of.PCs.transcriptome[i] + 1
            
            if(n_cols > limit) {
                col.to.exclude = list(NULL, c(1, (limit + 1):n_cols))
            } else {
                col.to.exclude = list(NULL, 1)
            }
        } else if (transcr.exclude.mode == 'first_1_column'){
            col.to.exclude = list(NULL, 1)
        }
        
        # Association calculation
        microb.N_PC.tissue.expr_PCs.corr.coef[[i]] <- check.PC.assoc.fast ( #check.PC.assoc(
            df1 = microbiome.PCs, 
            df2 = current_df,
            method = corr.method,
            col.to.exclude = col.to.exclude,
            #pval.adjust = F, 
            remove.outliers = remove.outliers,
            remove.outliers.IQR = remove.outliers.IQR
        )
    }
    
    # Efficiently flatten and name p-values
    p_values_vec <- unlist(lapply(seq_along(microb.N_PC.tissue.expr_PCs.corr.coef), function(i) {
        res_mat <- microb.N_PC.tissue.expr_PCs.corr.coef[[i]]
        p_vals <- as.vector(res_mat)
        
        # Reconstruct names: Tissue.mPC.trPC
        tissue_prefix <- word(expr_names[i], 1, sep = "\\.")
        m_names <- gsub('RS', 'mPC', rownames(res_mat))
        tr_names <- gsub('PC', 'trPC', colnames(res_mat))
        
        full_names <- paste(tissue_prefix, 
                            rep(m_names, times = length(tr_names)), 
                            rep(tr_names, each = length(m_names)), 
                            sep = ".")
        names(p_vals) <- full_names
        return(p_vals)
    }))
    
    if(inflation.correct){
        p_values_vec <- correct4inflation(p_values_vec)
    }

    if(!is.null(QQplot.file)){
        require(lattice)
        jpeg(QQplot.file,width=width,height=height)
        print(qqunif.plot(p_values_vec,conf.col="lightgray", conf.alpha=.05,
                    main=main,
                    par.strip.text = list(cex = 4, mgp = c(6,6,6,6)),
                    scales=list(tck=c(1,0), x=list(cex=2), y=list(cex=2)),
                    cex.main = cex.main, cex.x = 2, cex.y = 2
        ))
        dev.off()
    }
    
    return(p_values_vec)
}



# Function to estimate Neff using the Sidak Slope method
# p_nominal: The vector of p-values from your actual analysis
# perm_min_p: A vector containing the minimum p-value from each permutation
estimate_neff_sidak <- function(p_nominal, perm_min_p) {
  # Sort nominal p-values to create a range for evaluation
  # We focus on the lower end (e.g., p < 0.1) where the slope is most informative
  eval_p <- sort(p_nominal[p_nominal > 0 & p_nominal < 0.1])
  
  if(length(eval_p) < 10) {
    warning("Too few small p-values to estimate slope accurately. Using all p-values.")
    eval_p <- sort(p_nominal[p_nominal > 0 & p_nominal < 1])
  }

  # Calculate empirical P_expN for each nominal p
  # P_expN = Proportion of permutations where min(p) <= nominal p
  p_exp_n <- sapply(eval_p, function(p) mean(perm_min_p <= p))
  
  # Linearize: log(1 - P_expN) = N * log(1 - p)
  # Y = log(1 - P_expN), X = log(1 - p)
  # We must handle cases where p_exp_n is 1 to avoid -Inf
  valid_idx <- which(p_exp_n < 1 & p_exp_n > 0)
  y_vals <- log(1 - p_exp_n[valid_idx])
  x_vals <- log(1 - eval_p[valid_idx])
  
  # Fit linear model through the origin (intercept = 0)
  # The slope is our Neff  
  fit <- lm(y_vals ~ 0 + x_vals)  
  # output the relative error - as CV
  sum_fit <- summary(fit)
  neff_val <- sum_fit$coefficients[1, "Estimate"]
  neff_se <- sum_fit$coefficients[1, "Std. Error"]
  # Calculate the Coefficient of Variation (CV)
  # This tells you the relative uncertainty (e.g., 0.02 means 2% error)
  neff_cv <- neff_se / neff_val
  cat("Neff coefficient of variation: ", round(neff_cv, 4), "\n")

  neff <- as.numeric(coef(fit))
  
  return(neff)
}




# --- QQ plot helper ---

microb.transcr.make.QQplot <- function(pvalues, QQplot.file,
                                       inflation.correct = F, width = 800,
                                       height = 800, main = '', cex.main = 2){
    if(inflation.correct){
        pvalues <- correct4inflation(pvalues)
    }
    require(lattice)
    jpeg(QQplot.file,width=width,height=height)
    print(qqunif.plot(pvalues,conf.col="lightgray", conf.alpha=.05,
                      main=main,
                      par.strip.text = list(cex = 4, mgp = c(6,6,6,6)),
                      scales=list(tck=c(1,0), x=list(cex=2), y=list(cex=2)),
                      cex.main = cex.main, cex.x = 2, cex.y = 2
    ))
    dev.off()
    return(pvalues)
}

# --- Write p-value table to file ---
#write.pvalue.vectors.to.file <- function(output_file, Original, FDR = NULL,
#                                         Infl.corrected = NULL,
#                                         Infl.corrected.FDR = NULL,
#                                         signif_ = NULL){
#    if(!is.vector(Original)|length(Original)<1){
#        stop('Original must be a vector of non-zero length')
#    }
#    black.labels <- logical(4)
#    if(is.null(FDR)){
#        FDR <- rep(-2, times = length(Original))
#        names(FDR) <- names(Original)
#        black.labels[2] <- T
#    }
#    if(is.null(Infl.corrected)){
#        Infl.corrected <- rep(-2, times = length(Original))
#        names(Infl.corrected) <- names(Original)
#        black.labels[3] <- T
#    }
#    if(is.null(Infl.corrected.FDR)){
#        Infl.corrected.FDR <- rep(-2, times = length(Original))
#        names(Infl.corrected.FDR) <- names(Original)
#        black.labels[4] <- T
#    }
#    if (sum((names(Original)!=names(FDR))|(names(Original)!=names(Infl.corrected.FDR))|
#                                       (names(Original)!=names(Infl.corrected)))!=0){
#        stop('names of the 4 vectors should be same')
#    }
#    df <- matrix(ncol = 4, nrow = length(Original))
#    rownames(df) <- names(Original)
#    colnames(df) <- c('Original', 'FDR', 'Infl.corrected',
#                      'Infl.corrected.FDR')
#    if(!is.null(signif_)){
#        Original <- as.numeric(signif(Original, digits = signif_))
#        FDR <-  as.numeric(signif(FDR, digits = signif_))
#        Infl.corrected <-  as.numeric(signif(Infl.corrected, digits = signif_))
#        Infl.corrected.FDR <-  as.numeric(signif(Infl.corrected.FDR, digits = signif_))
#    }
#    for(i in 1:length(Original)){
#        df[i,] <-c(Original[i],FDR[i],Infl.corrected[i],Infl.corrected.FDR[i])
#    }
#    df <- as.data.frame(df)
#    df <- df[,!black.labels, drop = F]
#    write.table.smart(df,output_file = output_file)
#}


write.pvalue.vectors.to.file <- function(output_file, Original, 
                                         Hochberg = NULL,
                                         Sidak = NULL,
                                         Infl.corrected = NULL,
                                         Infl.corrected.Hochberg = NULL,
                                         signif_ = NULL){
    
    # Use length() instead of is.vector() to avoid the attribute trap
    if(is.null(Original) || length(Original) < 1){
        stop('Original must be a provided and have a non-zero length')
    }
    
    # Convert inputs to standard numeric vectors to strip non-name attributes 
    # while keeping the data and names intact for the dataframe
    clean_vec <- function(v) {
        if(is.null(v)) return(NULL)
        res <- as.numeric(v)
        names(res) <- names(v)
        return(res)
    }

    orig_clean <- clean_vec(Original)
    
    include_col <- rep(TRUE, 5)
    
    fill_null <- function(vec, ref, idx) {
        if(is.null(vec)) {
            include_col[idx] <<- FALSE
            v <- rep(-2, length(ref))
            names(v) <- names(ref)
            return(v)
        }
        return(clean_vec(vec))
    }

    hoch_v <- fill_null(Hochberg, orig_clean, 2)
    sidak_v <- fill_null(Sidak, orig_clean, 3)
    infl_v <- fill_null(Infl.corrected, orig_clean, 4)
    infl_hoch_v <- fill_null(Infl.corrected.Hochberg, orig_clean, 5)

    # Name consistency check
    all_names <- list(names(orig_clean), names(hoch_v), names(sidak_v), 
                      names(infl_v), names(infl_hoch_v))
    
    if (!all(sapply(all_names, function(x) identical(x, names(orig_clean))))) {
        stop('Names of all provided p-value vectors must be identical')
    }

    # Create the dataframe
    df <- data.frame(
        Original = orig_clean,
        Hochberg = hoch_v,
        Sidak = sidak_v,
        Infl.corrected = infl_v,
        Infl.corrected.Hochberg = infl_hoch_v,
        row.names = names(orig_clean)
    )

    # Apply significance rounding
    if(!is.null(signif_)){
        df <- as.data.frame(lapply(df, function(x) {
            # Only round values that aren't the placeholder (-2)
            ifelse(x == -2, x, signif(as.numeric(x), digits = signif_))
        }), row.names = rownames(df))
    }

    # Filter out the NULL columns and save
    df <- df[, include_col, drop = FALSE]
    write.table.smart(df, output_file = output_file)
}

# --- Correlation test: microbe vs transcript across tissues ---

# Tests correlation between a specific microbe and a specific transcript
# across all 9 tissue/cell types.
corrtest.microb.transcript <- function(taxa.abund_t,
    microbe = 'k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Rikenellaceae;g__Rikenellaceae.RC9.gut.group',
    transcr.dir = NULL,
    tissue.Res.file.name.end = '_Corrected4_Covars_PCs.txt',
    transcript = 'X1410114',
    corr.test.method = 'spearman',
    expression.Res.names.vector = c('IL_GE','TR_GE','RE_GE','CD4_GE',
                                     'CD8_GE','CD14_GE','CD15_GE',
                                     'CD19_GE','PLA_GE'),
    scatter.plot = F,
    scatter.plot.dir = NULL,
    scatter.plot.file.suffix = '_PID1expr__Rikenellaceae_RC9_gut_group_abund.jpg',
    ylab_suffix = ': PID1 gene expression',
    xlab_ = 'Rikenellaceae_RC9_gut_group',
    digits_ = 3
    ){
    if (is.null(transcr.dir)) stop('transcr.dir must be specified')
    result <- as.data.frame(matrix(ncol = 2))
    colnames(result) <- c('Corr.coef','p-value')
    for(tissue.Res in expression.Res.names.vector){
        tissue = word(tissue.Res, sep = '_', end = 1)
        print(paste0('Tissue: ', tissue.Res, ' is being processed'))
        current.tissue.Res.file <- paste0(tissue.Res, tissue.Res.file.name.end)
        # Checking that the expression file contains the probe of interest
        line1.current=readLines(con=file.path(transcr.dir, current.tissue.Res.file),
                                n=1)
        if((line1.current%like%paste0(' ',str_sub(transcript, start = 2),' '))|
            (line1.current%like%transcript)){
            current.tissue.Res <- read.table(file.path(transcr.dir, current.tissue.Res.file),
                                             row.names = 1, header = T)
            assign(tissue.Res,current.tissue.Res)
            if((sum(rownames(current.tissue.Res) != current.tissue.Res[,1]))!=0){
                stop('tissue Res file reading error')
            }
            taxa.abund_t.current.tissue.Res.merged <- merge.by.rownames.flex(
                taxa.abund_t, current.tissue.Res, all.x = F, all.y = F)
            taxa.abund_t.current.tissue.Res.merged[,microbe] <- remove.outliers(
                taxa.abund_t.current.tissue.Res.merged[,microbe])
            taxa.abund_t.current.tissue.Res.merged[,transcript] <-
                remove.outliers(taxa.abund_t.current.tissue.Res.merged[
                    ,transcript])
            cor.test.res <- cor.test(taxa.abund_t.current.tissue.Res.merged[
                ,colnames(taxa.abund_t.current.tissue.Res.merged)==microbe],
                taxa.abund_t.current.tissue.Res.merged[
                    ,colnames(taxa.abund_t.current.tissue.Res.merged)==transcript],
                method = corr.test.method)
            result <- rbind(result, c(signif(as.numeric(cor.test.res$estimate),
                                             digits = digits_),
                                      signif(as.numeric(cor.test.res$p.value),
                                                 digits = digits_))
                        )
            rownames(result)[nrow(result)] <- tissue
            if(scatter.plot){
                if (is.null(scatter.plot.dir)) stop('scatter.plot.dir must be specified')
                PC.scatter.plot.create(microbiome.df = taxa.abund_t.current.tissue.Res.merged[
                    colnames(taxa.abund_t.current.tissue.Res.merged)%like%'^k__'],
                    transcriptome.df = taxa.abund_t.current.tissue.Res.merged[
                        !colnames(taxa.abund_t.current.tissue.Res.merged)%like%'^k__'],
                    microbiome.PC = microbe,
                    transcriptome.PC = transcript, tissue = tissue,
                    ylab_ = paste0(tissue, ylab_suffix),
                    xlab_ = xlab_,
                    size.axis.title = 40,
                    file.out = paste0(scatter.plot.dir,'/',tissue,
                                      scatter.plot.file.suffix
                    )
                )
            }
        }
    }
    # removing the first NA row from result
    result <- result[-c(1),]
}

# --- Enterotype boxplot ---

draw.boxPlot.enterotypes <- function(ent.df, df2, variable,
                                     file_out, pdf = F,
                                     title_ = NULL,
                                     x_lab = 'Enterotype',
                                     y_lab = 'Gene expression',
                                     size_axes = 30, width_ = 800, height_ = 800,
                                     color_ = 'black',
                                     enterotypes = c('B-R','P','B','R'),
                                     remove.outliers.df2 = T,
                                     kruskal.test_ = F, signif.kruskal = 3,
                                     wilcox.test_ = F
                                     ){
    df2 <- df2[,colnames(df2)== variable,drop=F]
    if ((ncol(ent.df)!=1)|(ncol(df2)!=1)){
        stop('Input data frame error')
    }
    if(!is.numeric(signif.kruskal)){
        signif.kruskal <- 3
        warning('signif.kruskal must be numeric! Converted to 3')
    }
    ent.df.df2 <- merge.by.rownames.flex(ent.df,df2,all = F)
    if(remove.outliers.df2){
        ent.df.df2[,2]<- remove.outliers(ent.df.df2[,2])
    }
    ent.df.df2 <- ent.df.df2[!is.na(ent.df.df2[,1]),]
    ent.df.df2 <- ent.df.df2[!is.na(ent.df.df2[,2]),]
    colnames(ent.df.df2)[1] <- 'Enterotype'
    levels(ent.df.df2[,1]) <- c('Bact1','Prev','Bact2','Rum')
    colnames(ent.df.df2)[2] <- 'Gene'
    size_ = size_axes
    if(pdf){
        pdf(file_out, width = width_, height = height_)
    }else{
        jpeg(file_out, width = width_, height = height_)
    }
    print(ggplot(ent.df.df2, aes(x=Enterotype, y=Gene))+geom_boxplot()+ geom_jitter(
        shape=20, position=position_jitter(0.2), size = 3) +
        labs(title = title_, x = x_lab, y = y_lab)+
        theme_bw()+
        theme(plot.title = element_text(color = color_, size = size_,
                                        hjust = 0.5),
              axis.text.x = element_text(color = color_, size = size_),
              axis.text.y = element_text(color = color_, size = size_),
              axis.title.x = element_text(color = color_, size = size_),
              axis.title.y = element_text(color = color_, size = size_))
        )
    dev.off()
    if(kruskal.test_){
        print(paste0('Kruskal-Wallis test p-value ',
                     signif(kruskal.test(ent.df.df2[,2]~ent.df.df2[,1])$p.value,
                     digits = signif.kruskal)))
    }
    if(wilcox.test_){
        ent.df.df2.1 <- ent.df.df2
        ent.df.df2.1[,3] <- ifelse(ent.df.df2.1[,1]=='Bact1',1,0)
        ent.df.df2.1[,4] <- ifelse(ent.df.df2.1[,1]=='Prev',1,0)
        ent.df.df2.1[,5] <- ifelse(ent.df.df2.1[,1]=='Bact2',1,0)
        ent.df.df2.1[,6] <- ifelse(ent.df.df2.1[,1]=='Rum',1,0)

        print(paste0('Mann-Whitney U test p-value for ent Bact1 vs others: ',
                     signif(wilcox.test(ent.df.df2.1[,2]~ent.df.df2.1[,3])$p.value,
                            digits = signif.kruskal)))
        print(paste0('Mann-Whitney U test p-value for ent Prev vs others: ',
                     signif(wilcox.test(ent.df.df2.1[,2]~ent.df.df2.1[,4])$p.value,
                            digits = signif.kruskal)))
        print(paste0('Mann-Whitney U test p-value for ent Bact2 vs others: ',
                     signif(wilcox.test(ent.df.df2.1[,2]~ent.df.df2.1[,5])$p.value,
                            digits = signif.kruskal)))
        print(paste0('Mann-Whitney U test p-value for ent Rum vs others: ',
                     signif(wilcox.test(ent.df.df2.1[,2]~ent.df.df2.1[,6])$p.value,
                            digits = signif.kruskal)))
    }
}
