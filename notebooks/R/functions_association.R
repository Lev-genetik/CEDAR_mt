# functions_association.R
# Adapted from: code/CEDAR-master/Microbiome_genome_transcr_comp/Microbiota_QTL_functions_v9.R
# Functions for association testing between data frames (PCs, residuals,
# enterotypes, genotypes, taxa abundances).
#
# Core functions used in the microbiome-transcriptome pipeline:
#   check.PC.assoc()       — iterates over all column pairs, returns p-value matrix
#   test.lm.microb.gt()    — linear model p-values for genotype/PC x taxon/gene
#   df.remove.outliers()   — IQR-based outlier removal on data frame columns
#
# Additional GWAS functions (used in genotype analyses, not in core notebooks):
#   readGenotFile(), add.microbiota.info(), construct.gt.taxon.file(),
#   check.marker.locus(), GWAS.violin.plots(), get.pval.assoc.bact.marker(),
#   PLINK.plot(), manhattan.mod(), prepare_Hardi(),
#   plot.region.of.interest.GWAS(), plot.mGWAS.eQTL.together(),
#   plot.regional.eQTL.GTEx()
#
# Dependencies (sourced before this file by 00_setup.ipynb):
#   functions_general.R — merge.by.rownames.flex(), remove.outliers(),
#                         write.table.smart(), taxonomic.unit()

require(ggplot2)
require(dplyr)
require(tibble)
require(stats)
# Note: library(bedr) is required only for plot.region.of.interest.GWAS()
# and plot.regional.eQTL.GTEx(). Not loaded here to avoid failures in
# environments without bedr. Load it manually before calling those functions.

# ===========================================================================
# CORE FUNCTIONS
# ===========================================================================

# --- check.PC.assoc: pairwise association testing between two data frames ---
#
# Tests all column-pairs between df1 and df2 using correlation or lm or kruskal.
# INPUT:
#   df1, df2: data frames (samples in rows, variables in columns)
#   method: 'spearman', 'pearson', 'lm', or 'kruskal'
#   col.to.exclude: list of 2 (columns to exclude from df1 and df2)
#   remove.outliers: logical vector of 2 (outlier removal for df1 and df2)
#   remove.outliers.IQR: IQR multiplier for outlier detection
#   return: 'p.values' or 'correlation.coef'
# OUTPUT: matrix of p-values (or correlation coefficients)
check.PC.assoc <- function(df1,df2, method = 'spearman', return = 'p.values',
                           col.to.exclude = NULL,
                           pval.adjust = FALSE, remove.outliers = c(F,F),
                           remove.outliers.IQR = 1.5){
    # remove.outliers must be a vector of 2 values
    if(length(remove.outliers)==1){
        remove.outliers <- c(remove.outliers,remove.outliers)
    }
    if(!return%in%c('p.values','correlation.coef')){
        stop('return must be either pvalues or correlations')
    }
    # exclude samples that are not common and sort
    common.samples <- intersect(rownames(df1), rownames(df2))
    if (length(common.samples)<10){
        stop('<10 samples are in common across the 2 input data frames!')
    }
    df1 <- df1[common.samples,,drop=F]
    df2 <- df2[common.samples,]
    df1 <- df1[sort(row.names(df1)),,drop=F]
    df2 <- df2[sort(row.names(df2)),]
    if (!is.null(col.to.exclude[[1]])){df1 <- df1[,-c(col.to.exclude[[1]])]}
    if (!is.null(col.to.exclude[[2]])){df2 <- df2[,-c(col.to.exclude[[2]])]}
    # remove outliers if needed
    if(remove.outliers[1]|remove.outliers[2]){
        if(!is.numeric(remove.outliers.IQR)){
            stop('remove.outliers.IQR must be numeric')
        }
        if(remove.outliers[1]){
            for(i in 1:ncol(df1)){
                lh <- quantile(df1[,i],probs=0.25, na.rm = T)
                uh <- quantile(df1[,i],probs=0.75, na.rm = T)
                step <- remove.outliers.IQR * (uh-lh)
                df1[,i][df1[,i]>uh + step] <- NA
                df1[,i][df1[,i]<lh - step] <- NA
            }
        }
        if(remove.outliers[2]){
         for(i in 1:ncol(df2)){
            lh <- quantile(df2[,i],probs=0.25, na.rm = T)
            uh <- quantile(df2[,i],probs=0.75, na.rm = T)
            step <- remove.outliers.IQR * (uh-lh)
            df2[,i][df2[,i]>uh + step] <- NA
            df2[,i][df2[,i]<lh - step] <- NA
         }
        }
    }
    result <- matrix(NA, nrow = ncol(df1), ncol = ncol(df2))
    rownames(result) <- colnames(df1)
    colnames(result) <- colnames(df2)
    for(i in 1:ncol(df1)){
        for (j in 1:ncol(df2)){
            if (return=='p.values'){
                if(method %in% c('pearson', 'spearman')){
                    result[i,j] <- cor.test(df1[,i], df2[,j], method = method)$p.value
                } else if (method %like% 'kruskal'){
                    if(!is.factor(df1[,i])){
                        df1[,i] <- as.factor(df1[,i])
                        warning(paste0('df1[,',i,'] converted to factor'))
                    }
                    result[i,j] <- kruskal.test(x=df2[,j],g=df1[,i])$p.value
                } else if(method=='lm'){
                    result[i,j] <- anova(
                        lm(formula = df1[,i] ~ df2[,j], data = df1)
                    )[1,'Pr(>F)']
                }
            }
            if (return == 'correlation.coef'){
                if(method %in% c('pearson', 'spearman')){
                    result[i,j] <- as.numeric(cor.test(df1[,i], df2[,j], method = method)[[4]])
                } else{
                    stop('method must be pearson or spearman')
                }
            }
        }
    }
    if(pval.adjust){
        result <- p.adjust(result, method = 'hochberg')
    }
    result
}


# Speed-Optimized check.PC.assoc.fast
#This version vectorizes outlier removal and uses optimized matrix correlations.
check.PC.assoc.fast <- function(df1, df2, method = 'lm', 
                                col.to.exclude = NULL,
                                remove.outliers = c(F, F),
                                remove.outliers.IQR = 1.5) {
    
    # 1. Alignment (Should ideally be done once outside this function)
    common.samples <- intersect(rownames(df1), rownames(df2))
    df1 <- df1[common.samples, , drop = FALSE]
    df2 <- df2[common.samples, , drop = FALSE]
    
    if (!is.null(col.to.exclude[[1]])) df1 <- df1[, -col.to.exclude[[1]], drop = FALSE]
    if (!is.null(col.to.exclude[[2]])) df2 <- df2[, -col.to.exclude[[2]], drop = FALSE]

    # 2. Vectorized Outlier Removal
    # We use apply to avoid the for-loop over columns
    remove_outliers_vec <- function(x, iqr_mult) {
        qs <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
        step <- iqr_mult * (qs[2] - qs[1])
        x[x > qs[2] + step | x < qs[1] - step] <- NA
        return(x)
    }
    
    if(remove.outliers[1]) df1 <- apply(df1, 2, remove_outliers_vec, iqr_mult = remove.outliers.IQR)
    if(remove.outliers[2]) df2 <- apply(df2, 2, remove_outliers_vec, iqr_mult = remove.outliers.IQR)

    # 3. Fast Association Logic
    if (method == 'lm' || method == 'pearson' || method == 'spearman') {
        # Pearson and LM use the same p-value logic
        # Spearman is just Pearson on ranks
        if (method == 'spearman') {
            # use = 'pairwise.complete.obs' handles the NAs created by outlier removal
            r <- cor(df1, df2, method = 'spearman', use = 'pairwise.complete.obs')
        } else {
            r <- cor(df1, df2, method = 'pearson', use = 'pairwise.complete.obs')
        }
        
        # Calculate N for each pair (necessary if there are NAs)
        # crossprod on logical matrix gives the count of non-NA overlaps
        n_matrix <- crossprod(!is.na(df1), !is.na(df2))
        
        # Matrix-based p-value calculation
        t_stat <- r * sqrt((n_matrix - 2) / (1 - r^2))
        result <- 2 * pt(-abs(t_stat), df = n_matrix - 2)
        
        rownames(result) <- colnames(df1)
        colnames(result) <- colnames(df2)
        return(result)
        
    } else if (grepl('kruskal', method)) {
        # Kruskal is harder to vectorize fully, but we can optimize the loop
        result <- matrix(NA, nrow = ncol(df1), ncol = ncol(df2), 
                         dimnames = list(colnames(df1), colnames(df2)))
        for(i in 1:ncol(df1)) {
            g <- as.factor(df1[, i])
            for(j in 1:ncol(df2)) {
                result[i, j] <- kruskal.test(df2[, j] ~ g)$p.value
            }
        }
        return(result)
    }
}



# ultra fast implementation for LM case
check.PC.assoc.lm <- function(df1, df2) {
    # 1. High-speed correlation matrix
    # Pearson correlation p-value is identical to simple linear regression p-value
    r <- cor(df1, df2, use = "pairwise.complete.obs")
    
    # 2. Pairwise sample sizes (only needed if NAs are present)
    n_matrix <- crossprod(!is.na(df1), !is.na(df2))
    
    # 3. Vectorized t-distribution transformation
    # This replaces the double loop entirely
    t_stat <- r * sqrt((n_matrix - 2) / (1 - r^2))
    p_matrix <- 2 * pt(-abs(t_stat), df = n_matrix - 2)
    
    return(p_matrix)
}

## --- 1. Fast Matrix P-value Engine ---
## Mimics simple linear regression p-values using matrix correlation
#get_p_matrix_fast <- function(m_mat, t_mat) {
#  # Calculate Pearson correlation (handles NAs from outlier cleaning)
#  # 'pairwise.complete.obs' is necessary because outliers create NAs
#  r <- cor(m_mat, t_mat, use = "pairwise.complete.obs")
#  
#  # Calculate pairwise sample sizes (N) to account for NAs
#  n_matrix <- crossprod(!is.na(m_mat), !is.na(t_mat))
#  
#  # Vectorized transformation to t-statistic and then to p-value
#  t_stat <- r * sqrt((n_matrix - 2) / (1 - r^2))
#  p_matrix <- 2 * pt(-abs(t_stat), df = n_matrix - 2)
#  
#  return(p_matrix)
#}

# --- Helper: Fast Matrix P-value Engine ---
# Efficiently calculates p-values for all pairs using Pearson/LM equivalence
get_p_matrix_fast <- function(m_mat, t_mat) {
    # Handles NAs from outlier cleaning
    r <- cor(m_mat, t_mat, use = "pairwise.complete.obs")
    n_matrix <- crossprod(!is.na(m_mat), !is.na(t_mat))
    
    # Vectorized p-value calculation
    t_stat <- r * sqrt((n_matrix - 2) / (1 - r^2))
    p_matrix <- 2 * pt(-abs(t_stat), df = n_matrix - 2)
    return(p_matrix)
}


# Fast vectorized Kruskal-Wallis for 1 categorical vs many numeric features
get_p_kruskal_fast <- function(group_factor, data_matrix) {
    # Remove NAs from the group factor
    idx <- !is.na(group_factor)
    g <- group_factor[idx]
    m <- data_matrix[idx, , drop=FALSE]
    
    n <- nrow(m)
    r <- apply(m, 2, rank)
    
    # Kruskal formula: H = (12 / n(n+1)) * sum(T_i^2 / n_i) - 3(n+1)
    group_levels <- unique(g)
    sum_sq_ranks <- 0
    for(lvl in group_levels) {
        group_idx <- g == lvl
        n_i <- sum(group_idx)
        T_i <- colSums(r[group_idx, , drop=FALSE])
        sum_sq_ranks <- sum_sq_ranks + (T_i^2 / n_i)
    }
    
    h_stat <- (12 / (n * (n + 1))) * sum_sq_ranks - 3 * (n + 1)
    df <- length(group_levels) - 1
    p_vals <- pchisq(h_stat, df = df, lower.tail = FALSE)
    
    return(p_vals)
}




prepare_data_universal <- function(microb_df, transcr_list, 
                                   n_pcs = NULL, 
                                   exclude_mode = 'last_columns',
                                   outlier_modes = c(F, T), 
                                   iqr_mult = 5) {
    
    # Internal helper to clean numeric data only
    clean_numeric_only <- function(df, mult) {
        # Identify which columns are actually numeric
        is_num <- sapply(df, is.numeric)
        if(any(is_num)) {
            df[, is_num] <- lapply(df[, is_num, drop=FALSE], function(x) {
                qs <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
                step <- mult * (qs[2] - qs[1])
                x[x > qs[2] + step | x < qs[1] - step] <- NA
                return(x)
            })
        }
        # Return only the numeric parts to ensure matrix safety
        return(df[, is_num, drop=FALSE])
    }

    prepared <- list()
    for (nm in names(transcr_list)) {
        df_t <- transcr_list[[nm]]
        
        # A. Exact logic from your original function for column subsetting
        if (exclude_mode == 'last_columns' && !is.null(n_pcs)) {
            # Your original code kept columns 2 through (N+1)
            # (Excluding col 1 and anything after N)
            limit <- n_pcs[nm]
            # Safety check: don't exceed ncol
            end_col <- min(limit + 1, ncol(df_t))
            if(ncol(df_t) >= 2) {
                df_t <- df_t[, 2:end_col, drop = FALSE]
            }
        } else if (exclude_mode == 'first_1_column') {
            # Exclude only the first column
            if(ncol(df_t) > 1) df_t <- df_t[, -1, drop = FALSE]
        }

        # B. Sample Alignment
        common <- intersect(rownames(microb_df), rownames(df_t))
        if(length(common) < 10) next # Skip if too few samples
        
        m_sub <- microb_df[common, , drop=FALSE]
        t_sub <- df_t[common, , drop=FALSE]
        
        # C. Outlier Cleaning AND Numeric Filtering
        # This step now drops non-numeric columns automatically
        m_clean <- if(outlier_modes[1]) clean_numeric_only(m_sub, iqr_mult) else m_sub[, sapply(m_sub, is.numeric), drop=FALSE]
        t_clean <- if(outlier_modes[2]) clean_numeric_only(t_sub, iqr_mult) else t_sub[, sapply(t_sub, is.numeric), drop=FALSE]
        
        # D. Store as strictly numeric matrices
        prepared[[nm]] <- list(m = as.matrix(m_clean), t = as.matrix(t_clean))
    }
    return(prepared)
}


prepare_data_universal_ent <- function(microb_df, transcr_list, 
                                   n_pcs = NULL, 
                                   exclude_mode = 'last_columns',
                                   outlier_modes = c(F, T), 
                                   iqr_mult = 5) {
    
    clean_numeric_only <- function(df, mult) {
        is_num <- sapply(df, is.numeric)
        if(any(is_num)) {
            df[, is_num] <- lapply(df[, is_num, drop=FALSE], function(x) {
                qs <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
                step <- mult * (qs[2] - qs[1])
                x[x > qs[2] + step | x < qs[1] - step] <- NA
                return(x)
            })
        }
        # For cleaning, we only return numeric. For Enterotypes, we won't trigger this.
        return(df[, is_num, drop=FALSE])
    }

    prepared <- list()
    for (nm in names(transcr_list)) {
        df_t <- transcr_list[[nm]]
        
        # A. Column Selection
        if (exclude_mode == 'last_columns' && !is.null(n_pcs)) {
            limit <- n_pcs[nm]
            end_col <- min(limit + 1, ncol(df_t))
            if(ncol(df_t) >= 2) df_t <- df_t[, 2:end_col, drop = FALSE]
        } else if (exclude_mode == 'first_1_column') {
            if(ncol(df_t) > 1) df_t <- df_t[, -1, drop = FALSE]
        }

        # B. Sample Alignment
        common <- intersect(rownames(microb_df), rownames(df_t))
        if(length(common) < 5) next 
        
        m_sub <- microb_df[common, , drop=FALSE]
        t_sub <- df_t[common, , drop=FALSE]
        
        # C. SMART CLEANING: 
        # If outlier_modes is FALSE, we KEEP the column (even if it's a Factor/Enterotype)
        m_clean <- if(outlier_modes[1]) clean_numeric_only(m_sub, iqr_mult) else m_sub
        t_clean <- if(outlier_modes[2]) clean_numeric_only(t_sub, iqr_mult) else t_sub
        
        # Store. t must be a matrix for math. m can be a data.frame to keep the Factor.
        prepared[[nm]] <- list(m = m_clean, t = as.matrix(t_clean))
    }
    return(prepared)
}


# Association Engine: Preserves your specific naming convention
#run_universal_association <- function(prepared_sets, m_prefix = "mPC", t_prefix = "trPC") {
#    results <- unlist(lapply(names(prepared_sets), function(nm) {
#        s <- prepared_sets[[nm]]
#        p_mat <- get_p_matrix_fast(s$m, s$t)
#        
#        # Naming Logic: Tissue.mPC.trPC
#        row_n <- gsub('RS', m_prefix, rownames(p_mat))
#        col_n <- gsub('PC', t_prefix, colnames(p_mat))
#        
#        tis_prefix <- stringr::word(nm, 1, sep="\\.")
#        full_names <- paste(tis_prefix, 
#                            rep(row_n, times = ncol(p_mat)), 
#                            rep(col_n, each = nrow(p_mat)), sep = ".")
#        
#        p_vec <- as.vector(p_mat)
#        names(p_vec) <- full_names
#        return(p_vec)
#    }))
#    return(results)
#}

#run_universal_association <- function(prepared_sets, m_prefix = "mPC", t_prefix = "trPC", 
#                                      inflation_correct = FALSE) {
#    results <- unlist(lapply(names(prepared_sets), function(nm) {
#        s <- prepared_sets[[nm]]
#        p_mat <- get_p_matrix_fast(s$m, s$t)
#        
#        row_n <- gsub('RS', m_prefix, rownames(p_mat))
#        col_n <- gsub('PC', t_prefix, colnames(p_mat))
#        tis_prefix <- stringr::word(nm, 1, sep="\\.")
#        
#        p_vec <- as.vector(p_mat)
#        
#        # Apply inflation correction BEFORE naming/returning
#        if(inflation_correct) {
#            p_vec <- correct4inflation(p_vec)
#        }
#        
#        names(p_vec) <- paste(tis_prefix, 
#                             rep(row_n, times = ncol(p_mat)), 
#                             rep(col_n, each = nrow(p_mat)), sep = ".")
#        return(p_vec)
#    }))
#    return(results)
#}

#run_universal_association <- function(prepared_sets, m_prefix = "mPC", t_prefix = "trPC", 
#                                      inflation_correct = FALSE) {
#    # 1. Calculate all p-values as before
#    results_list <- lapply(names(prepared_sets), function(nm) {
#        s <- prepared_sets[[nm]]
#        p_mat <- get_p_matrix_fast(s$m, s$t)
#        
#        row_n <- gsub('RS', m_prefix, rownames(p_mat))
#        col_n <- gsub('PC', t_prefix, colnames(p_mat))
#        tis_prefix <- stringr::word(nm, 1, sep="\\.")
#        
#        p_vec <- as.vector(p_mat)
#        names(p_vec) <- paste(tis_prefix, rep(row_n, times = ncol(p_mat)), 
#                             rep(col_n, each = nrow(p_mat)), sep = ".")
#        return(p_vec)
#    })
#    
#    p_all <- unlist(results_list)
#    
#    # 2. Track and apply Lambda
#    current_lambda <- 1.0
#    if (inflation_correct) {
#        current_lambda <- calculate_lambda(p_all)
#        p_all <- apply_fixed_lambda(p_all, current_lambda)
#    }
#    
#    # 3. Attach lambda as an attribute for the Neff pipeline to use
#    attr(p_all, "lambda") <- current_lambda
#    return(p_all)
#}

run_universal_association <- function(prepared_sets, m_prefix = "mPC", t_prefix = "trPC", 
                                      inflation_correct = FALSE,
                                      QQplot.file = NULL, main = "Associations",
                                      width = 10, height = 10) {
    # 1. Calculate P-values
    results_list <- lapply(names(prepared_sets), function(nm) {
        s <- prepared_sets[[nm]]
        p_mat <- get_p_matrix_fast(s$m, s$t)
        
        row_n <- gsub('RS', m_prefix, rownames(p_mat))
        col_n <- gsub('PC', t_prefix, colnames(p_mat))
        tis_prefix <- stringr::word(nm, 1, sep="\\.")
        
        p_vec <- as.vector(p_mat)
        names(p_vec) <- paste(tis_prefix, rep(row_n, times = ncol(p_mat)), 
                             rep(col_n, each = nrow(p_mat)), sep = ".")
        return(p_vec)
    })
    p_all <- unlist(results_list)
    
    # 2. Handle Inflation
    lambda_val <- 1.0
    if (inflation_correct) {
        lambda_val <- calculate_lambda(p_all)
        p_all <- apply_fixed_lambda(p_all, lambda_val)
    }
    
    # 3. Handle QQ-plot
    if (!is.null(QQplot.file)) {
        title_full <- paste0(main, "\n(N=", length(p_all), ", lambda=", round(lambda_val, 3), ")")
        pdf(QQplot.file, width = width, height = height)
        print(qqunif.plot(p_all, main = title_full, cex.main = 1.5))
        dev.off()
    }
    
    if (inflation_correct) {
        attr(p_all, "lambda") <- lambda_val
    }
    return(p_all)
}



run_universal_association_ent <- function(prepared_sets, m_prefix = "Ent", t_prefix = "", 
                                          inflation_correct = FALSE) {
    results_list <- lapply(names(prepared_sets), function(nm) {
        s <- prepared_sets[[nm]]
        
        # Enterotypes are usually a single column in s$m
        # We take the first column (the Enterotype factor)
        p_vec <- get_p_kruskal_fast(s$m[,1], s$t)
        
        tis_prefix <- stringr::word(nm, 1, sep="\\.")
        tr_names <- if(t_prefix == "") colnames(s$t) else paste0(t_prefix, colnames(s$t))
        names(p_vec) <- paste(tis_prefix, m_prefix, tr_names, sep = ".")
        
        return(p_vec)
    })
    
    p_all <- unlist(results_list)
    
    # Handle Inflation and Attributes
    current_lambda <- 1.0
    if (inflation_correct) {
        current_lambda <- calculate_lambda(p_all)
        p_all <- apply_fixed_lambda(p_all, current_lambda)
    }
    
    attr(p_all, "lambda") <- current_lambda
    return(p_all)
}



# Function: run_neff_pipeline
# Input:
#   pval_observed: The vector of p-values from the real data
#   prepared_data: The list of aligned/cleaned [m, t] matrices (from prepare_data_universal)
#   microb_base: The microbiome dataframe to be shuffled (e.g. microbiome.15PCs)
#   n_perm: Number of permutations
#   n_cores: Number of CPUs to use
run_neff_pipeline <- function(pval_observed, prepared_data, microb_base, 
                              n_perm = 1000, n_cores = detectCores() - 1,
                              #inflation_correct = FALSE
                              fixed_lambda = NULL
                              ) {
    
    cat(paste("Starting Neff estimation with", n_perm, "permutations on", n_cores, "cores...\n"))
    t_start <- proc.time()
    
    # --- SET MASTER THREADS TO 1 ---
    # This ensures the main process isn't competing with workers
    blas_set_num_threads(1)
    omp_set_num_threads(1)
    # Setup Parallel Cluster
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    
    # Parallel Permutation Loop
    perm_min_ps <- foreach(i = 1:n_perm, .combine = 'c', 
                           #.export = c('get_p_matrix_fast', 'correct4inflation'),
                           .export = c('get_p_matrix_fast', 'apply_fixed_lambda'),
                           .packages = c('stringr', 'RhpcBLASctl')) %dopar% {
                           #.packages = 'stringr') %dopar% {

        # --- SET WORKER THREADS TO 1 ---
        # Each worker must explicitly set its own BLAS limit
        blas_set_num_threads(1)
        omp_set_num_threads(1)
        
        # Shuffle microbiome indices once per permutation
        shuffled_idx <- sample(nrow(microb_base))
        m_shuffled_full <- microb_base[shuffled_idx, , drop=FALSE]
        rownames(m_shuffled_full) <- rownames(microb_base)
        
        # Calculate p-values across all tissues/cells in this permutation
        #all_p_perm <- unlist(lapply(prepared_data, function(s) {
        #    # Re-align the shuffled microbiome to the samples in this specific tissue/cell set
        #    m_p <- m_shuffled_full[rownames(s$t), , drop=FALSE]
        #    # Fast matrix association
        #    as.vector(get_p_matrix_fast(m_p, s$t))
        #}))
        
        #all_p_perm <- unlist(lapply(prepared_data, function(s) {
        #    m_p <- m_shuffled_full[rownames(s$t), , drop=FALSE]
        #    p_vec <- as.vector(get_p_matrix_fast(m_p, s$t))
        #    # Match the inflation correction used in the real data
        #    if(inflation_correct) p_vec <- correct4inflation(p_vec)
        #    return(p_vec)
        #}))

        all_p_perm <- unlist(lapply(prepared_data, function(s) {
            m_p <- m_shuffled_full[rownames(s$t), , drop=FALSE]
            p_vec <- as.vector(get_p_matrix_fast(m_p, s$t))            
            # Apply the SAME correction factor used in the real data
            if(!is.null(fixed_lambda)) {
                p_vec <- apply_fixed_lambda(p_vec, fixed_lambda)
            }
            return(p_vec)
        }))
        
        # Return the "best" p-value for this permutation
        min(all_p_perm, na.rm = TRUE)
    }
    stopCluster(cl)
    
    # 3. Estimate Neff using the Sidak Slope method
    neff_val <- estimate_neff_sidak(pval_observed, perm_min_ps)
    
    # 4. Apply Sidak adjustment: 1 - (1 - p)^Neff
    pval_sidak <- 1 - (1 - pval_observed)^neff_val

    # add Hochberg adj
    pval_hochberg <- p.adjust(pval_observed, method = 'hochberg')    
    
    # 5. Reporting
    t_end <- proc.time() - t_start
    cat(paste0("Neff Estimation Complete (", round(t_end[3], 1), "s)\n"))
    cat(paste("Total Tests:", length(pval_observed), "\n"))
    cat(paste("Effective Independent Tests (Neff):", round(neff_val, 2), "\n"))
    cat(paste("Reduction in penalty:", round((1 - neff_val/length(pval_observed))*100, 1), "%\n"))
    cat('\np Sidak <= 0.15:', sum(pval_sidak <= 0.15, na.rm = TRUE), '\n')
    print(pval_sidak[pval_sidak <= 0.15])
    cat('\np Hochberg <= 0.15:', sum(pval_hochberg <= 0.15, na.rm = TRUE), '\n')
    print(pval_hochberg[pval_hochberg <= 0.15])
    
    
    return(list(
        pval_sidak = pval_sidak,
        pval_hochberg = pval_hochberg,
        neff = neff_val,
        perm_distribution = perm_min_ps
    ))
}


run_neff_pipeline_ent <- function(pval_observed, prepared_data, microb_base, 
                                  n_perm = 1000, n_cores = detectCores() - 1,
                                  fixed_lambda = NULL) {
    
    cat(paste("Starting Neff estimation (Kruskal-Wallis) with", n_perm, "permutations...\n"))
    
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    
    # We strictly export the Kruskal engine and the lambda helper
    perm_min_ps <- foreach(i = 1:n_perm, .combine = 'c', 
                           .export = c('get_p_kruskal_fast', 'apply_fixed_lambda'),
                           .packages = c('stringr', 'RhpcBLASctl')) %dopar% {
        
        RhpcBLASctl::blas_set_num_threads(1)
        RhpcBLASctl::omp_set_num_threads(1)
        
        # Shuffle enterotype labels
        shuffled_idx <- sample(nrow(microb_base))
        m_shuffled_full <- microb_base[shuffled_idx, , drop = FALSE]
        rownames(m_shuffled_full) <- rownames(microb_base)
        
        all_p_perm <- unlist(lapply(prepared_data, function(s) {
            # Align shuffled Enterotypes to the current transcriptome set
            m_p <- m_shuffled_full[rownames(s$t), , drop = FALSE]
            
            # Use Kruskal-Wallis (m_p[, 1] is the enterotype factor)
            p_vec <- as.vector(get_p_kruskal_fast(m_p[, 1], s$t))
            
            if (!is.null(fixed_lambda)) {
                p_vec <- apply_fixed_lambda(p_vec, fixed_lambda)
            }
            return(p_vec)
        }))
        
        min(all_p_perm, na.rm = TRUE)
    }
    stopCluster(cl)
    
    neff_val <- estimate_neff_sidak(pval_observed, perm_min_ps)
    pval_sidak <- 1 - (1 - pval_observed)^neff_val

    # add Hochberg adj
    pval_hochberg <- p.adjust(pval_observed, method = 'hochberg')    
    
    # 5. Reporting        
    cat(paste("Total Tests:", length(pval_observed), "\n"))
    cat(paste("Effective Independent Tests (Neff):", round(neff_val, 2), "\n"))
    cat(paste("Reduction in penalty:", round((1 - neff_val/length(pval_observed))*100, 1), "%\n"))
    cat('\np Sidak <= 0.15:', sum(pval_sidak <= 0.15, na.rm = TRUE), '\n')
    print(pval_sidak[pval_sidak <= 0.15])
    cat('\np Hochberg <= 0.15:', sum(pval_hochberg <= 0.15, na.rm = TRUE), '\n')
    print(pval_hochberg[pval_hochberg <= 0.15])

    
    return(list(pval_sidak = pval_sidak, pval_hochberg = pval_hochberg, neff = neff_val, perm_distribution = perm_min_ps))
}


# --- Function to calculate the Genomic Inflation Factor (Lambda) ---
calculate_lambda <- function(p_vec) {
    # Remove NAs
    p_vec <- p_vec[!is.na(p_vec)]
    # Convert p-values to Chi-squared statistics (1 degree of freedom)
    chisq <- qchisq(p_vec, df = 1, lower.tail = FALSE)
    # Calculate Lambda: Median(observed) / Median(expected under null)
    lambda <- median(chisq) / qchisq(0.5, 1)
    return(lambda)
}

# Helper to apply a pre-calculated lambda
apply_fixed_lambda <- function(p_vec, lambda) {
    chisq <- qchisq(p_vec, df = 1, lower.tail = FALSE)
    chisq_corr <- chisq / lambda
    p_corr <- pchisq(chisq_corr, df = 1, lower.tail = FALSE)
    return(p_corr)
}


# --- Standard QQ-plot for uniform distribution ---
qqunif.plot <- function(p_vec, conf.col="lightgray", conf.alpha=.05, main="QQ-plot", ...) {
    p_vec <- p_vec[!is.na(p_vec)]
    n <- length(p_vec)
    exp_p <- -log10((1:n) / (n + 1))
    obs_p <- -log10(sort(p_vec))
    
    # Simple lattice-based plot (optimized for many points)
    xyplot(obs_p ~ exp_p, 
           aspect = "iso", 
           xlab = expression(Expected~~-log[10](p)), 
           ylab = expression(Observed~~-log[10](p)),
           main = main,
           panel = function(x, y, ...) {
               panel.abline(0, 1, col = "red")
               panel.xyplot(x, y, pch = 20, ...)
           }, ...)
}

# --- 2. Outlier Cleaning Function ---
# Returns the matrix with outliers replaced by NA
clean_outliers_matrix <- function(mat, iqr_mult = 1.5) {
  apply(mat, 2, function(x) {
    qs <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
    step <- iqr_mult * (qs[2] - qs[1])
    x[x > qs[2] + step | x < qs[1] - step] <- NA
    return(x)
  })
}






# --- df.remove.outliers: IQR-based outlier removal on data frame ---

df.remove.outliers <- function(df, remove.outliers.IQR = 1.5){
    if(!is.numeric(remove.outliers.IQR)){
        stop('remove.outliers.IQR must be numeric')
    }
    for(i in 1:ncol(df)){
        lh <- quantile(df[,i],probs=0.25, na.rm = T)
        uh <- quantile(df[,i],probs=0.75, na.rm = T)
        step <- remove.outliers.IQR * (uh-lh)
        df[,i][df[,i]>uh + step] <- NA
        df[,i][df[,i]<lh - step] <- NA
    }
    return(df)
}

# ===========================================================================
# GENOTYPE-TAXON ASSOCIATION FUNCTIONS
# ===========================================================================

# --- Read genotype file ---

readGenotFile <- function(gt.file, convert.genotype.to.numeric = T,
                          samples.limit = NULL){
    genotypes <- readLines(gt.file)
    genotypes.list <- list()
    result <- as.data.frame(matrix(ncol = 1, nrow = length(
        unlist(strsplit(genotypes[1], split = ' ')))),stringsAsFactors = F)
    for (i in 1:length(genotypes)){
        genotypes.list[[i]] <- unlist(strsplit(genotypes[i], split = ' '))
        result <- cbind(result,unlist(genotypes.list[[i]]), stringsAsFactors = F)
    }
    rownames(result) <- result[,2]
    result <- result[,-c(1,2)]
    for (i in 1:ncol(result)){
        colnames(result)[i] <- as.character(result[1,i])
    }
    result <- result[-1,]
    result <- as.data.frame(result, stringsAsFactors = F)
    if(convert.genotype.to.numeric){
        for (i in 1:ncol(result)){
            result[,i] <- as.numeric(result[,i])
        }
    }
    if(!is.null(samples.limit)){
        result <- result[rownames(result)%in%samples.limit,]
    }
     result
}

# --- Add microbiota info to a genotype data frame ---

add.microbiota.info <- function(input.df, taxa.abundance.file, locaton = 'IL'){
    taxa.abundance <- read.table(taxa.abundance.file, header = T, row.names = 1,
                                 sep = '\t', as.is = T, check.names = F)
    if (locaton!='all'){
        taxa.abundance <- taxa.abundance[,colnames(taxa.abundance)%like%locaton]
    }
    taxa.abundance <- t(taxa.abundance)
    result.df <- input.df
    result.df[,(ncol(input.df)+1):(ncol(input.df)+ncol(taxa.abundance))] <- NA
    colnames(result.df)[(ncol(result.df)-ncol(taxa.abundance)+1):(ncol(result.df))] <-
        colnames(taxa.abundance)
    for (i in 1:nrow(result.df)){
        row.in.taxa.abundance <- which(rownames(taxa.abundance)%like%
                                           rownames(result.df)[i])
        if(length(row.in.taxa.abundance) == 0){
            next()
        }
        for (j in 1:ncol(taxa.abundance)){
            result.df[i,(ncol(input.df)+j)] <- taxa.abundance[
                row.in.taxa.abundance,j]
        }
    }
    result.df
}

# --- Construct genotype + taxon data frame ---

construct.gt.taxon.file <- function(GenotFile, abundance.file,
                                    abundance.file.2 = NULL, add3locations = F,
                                    amplicon = NULL, select.locations = NULL,
                                    samples.in.columns.abundance.files = c(T,F),
                                    convert.genotype.to.numeric = T,
                                    genotype.file.contains.alleles = F,
                                    sep.abundance.file = '\t'
                                    ){
    genotypes.df <- readGenotFile(GenotFile,
                                  convert.genotype.to.numeric = convert.genotype.to.numeric)
    if(genotype.file.contains.alleles){
        genotypes.df <- genotypes.df[-c(1,2),]
    }
    if(add3locations){
        genotypes.df.IL <- genotypes.df
        rownames(genotypes.df.IL) <- paste0(rownames(genotypes.df),'.IL')
        genotypes.df.TR <- genotypes.df
        rownames(genotypes.df.TR) <- paste0(rownames(genotypes.df),'.TR')
        genotypes.df.RE <- genotypes.df
        rownames(genotypes.df.RE) <- paste0(rownames(genotypes.df),'.RE')
        genotypes.df <- rbind(genotypes.df.IL,genotypes.df.TR,genotypes.df.RE)
    }
    if(!is.null(amplicon)){
        rownames(genotypes.df) <- paste0(rownames(genotypes.df),amplicon)
    }
    taxa.abundance <- read.table(file = abundance.file, header = T,
                                 row.names = 1, sep = sep.abundance.file, as.is = T,
                                 check.names = F)
    if(samples.in.columns.abundance.files[1]){
        taxa.abundance <- as.data.frame(t(taxa.abundance), stringsAsFactors = F)
    }
    result <- merge.by.rownames.flex(genotypes.df,taxa.abundance,all.x = T)
    if(!is.null(abundance.file.2)){
        taxa.abundance.2 <- read.table(file = abundance.file.2, header = T,
                                       row.names = 1, sep = '\t', as.is = T,
                                       check.names = F)
        if(samples.in.columns.abundance.files[2]){
            taxa.abundance.2 <- as.data.frame(t(taxa.abundance.2), stringsAsFactors = F)
        }
        result <- merge.by.rownames.flex(result,taxa.abundance.2,all.x = T)
    }
    if(!is.null(select.locations)){
        result_filtered <- NULL
        for(loc in select.locations){
            result_current <- result[rownames(result)%like%loc,]
            result_filtered <- rbind(result_filtered,result_current)
        }
        result <- result_filtered
    }
    result
}

# --- test.lm.microb.gt: linear model test for genotype vs taxon ---

test.lm.microb.gt <- function(SNP.microb.df, microbe_prefix = 'k__|RS[0-9]',
                              SNP_prefix = '^[0-9|X|Y]', return.SNP.lm.coef = F,
                              recessive.model = NULL, check.numeric = F){
    if(check.numeric){
        for(i in (1:ncol(SNP.microb.df))){
            if(!is.numeric(SNP.microb.df[,i])){
                warning(paste0('While checking associations, row ',i,' has
                               non-numeric values'))
            }
        }
    }
    ncol_result <- sum(colnames(SNP.microb.df)%like%microbe_prefix)
    nrow_result <- sum(colnames(SNP.microb.df)%like%SNP_prefix)
    if (min(which(colnames(SNP.microb.df)%like%microbe_prefix))-
        max(which(colnames(SNP.microb.df)%like%SNP_prefix))!=1){
        stop('genotypes should be in the input data frame immediately before bacterial taxa!')
    }
    if(min(which(colnames(SNP.microb.df)%like%SNP_prefix))!=1){
        stop('SNP genotypes in the SNP.microb.df must start from 1st column')
    }
    result <- matrix(2,ncol = ncol_result, nrow = nrow_result)
    colnames(result) <- colnames(
        SNP.microb.df)[(nrow_result+1):(nrow_result+ncol_result)]
    rownames(result) <- colnames(
        SNP.microb.df)[1:nrow_result]
    if(return.SNP.lm.coef){
        SNP.coef <- result
    }
    if(!is.null(recessive.model)){
        for(i in 1:nrow_result){
            SNP.microb.df[,i] <- ifelse(SNP.microb.df[,i]>recessive.model,1,0)
        }
    }
    for (i in 1:nrow_result){
        for (j in (nrow_result+1):(nrow_result+ncol_result)){
            result[i,j-nrow_result] <- anova(lm(formula = SNP.microb.df[,j] ~
                                                    SNP.microb.df[,i],
                                                data = SNP.microb.df))['Pr(>F)'][[1]][1]
            if(return.SNP.lm.coef){
            SNP.coef[i,j-nrow_result] <-lm(formula = SNP.microb.df[,j] ~
                                                      SNP.microb.df[,i],
                                                  data = SNP.microb.df)[[1]][2]
            }
        }
    }
    if(return.SNP.lm.coef){
        pvalues = result
        result <- list(pvalues,SNP.coef)
        names(result) <- c('pvalues','SNP.coef')
        result
    } else{
        result
    }
}

# --- check.marker.locus: check LCT or other marker locus association ---

check.marker.locus <- function(abundance.file, recessive.model = NULL,
                               marker.genotype.file = NULL,
                               SNP = '2_136608646_G_A', violin.plot = F,
                               violin.plot.file = NULL,
                               bacterium = 'k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Bifidobacteriales;f__Bifidobacteriaceae;g__Bifidobacterium',
                               allele1 = 'GG', allele2 = 'AA + AG',
                               width = 1000, height = 600,
                               y_label = "Scaled bacterium abundance",
                               x_label = "Genotype", title = NULL,
                               ...
){
    if (is.null(marker.genotype.file)) stop('marker.genotype.file must be specified')
    genotype.abuncance.df <- construct.gt.taxon.file(GenotFile =
                                                         marker.genotype.file,
                                                     abundance.file = abundance.file,...)
    pval_lm_taxa_SNP <- test.lm.microb.gt(
        SNP.microb.df = genotype.abuncance.df, return.SNP.lm.coef = TRUE,
        recessive.model = recessive.model, check.numeric = T)
    if (sum(rownames(pval_lm_taxa_SNP$SNP.coef)%in%SNP)!=1){
        stop ('SNP not found or found > 1 times in the data')
    }
    if (sum(colnames(pval_lm_taxa_SNP$SNP.coef)%in%bacterium)!=1){
        stop ('bacterium not found or found > 1 times in the data')
    }
    result <- paste0(taxonomic.unit(bacterium),' ~ ',SNP, ifelse(
        recessive.model,' recessive model',''), ' regression coef: ',
        signif(pval_lm_taxa_SNP$SNP.coef[SNP,bacterium], digits = 3),
        ', p-value: ', signif(pval_lm_taxa_SNP$pvalues[SNP,bacterium], digits = 3)
    )
    if(violin.plot == T){
        if(is.null(title)){
            title = taxonomic.unit(bacterium)
        }
        jpeg(violin.plot.file, width = width, height = height)
        par(mar = c(12,6,6,2))
        genotype.abuncance.df_ = genotype.abuncance.df
        GG <- genotype.abuncance.df_[,SNP]<recessive.model
        genotype.abuncance.df_[GG,SNP] <- allele1
        genotype.abuncance.df_[!GG,SNP] <- allele2
        genotype.abuncance.df_[,SNP] <- as.factor(genotype.abuncance.df_[,SNP])
        colnames(genotype.abuncance.df_)[colnames(genotype.abuncance.df_) ==
                                             bacterium] <- 'the_bacterium'
        colnames(genotype.abuncance.df_)[colnames(genotype.abuncance.df_) ==
                                             SNP] <- 'SNP'
        genotype.abuncance.df_ <- genotype.abuncance.df_[
            !is.na(genotype.abuncance.df_$the_bacterium),]
        p <- ggplot(genotype.abuncance.df_, aes(
            x=SNP, y=the_bacterium)) + scale_x_discrete(limits=c(allele1, allele2))+
            geom_violin(trim=FALSE, fill='cadetblue2', color="cadetblue3")
        print(p + geom_boxplot(width=0.6, col = 'darkblue', fill = 'cadetblue2',
                         na.rm = T) +
            geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4, col = 'cornflowerblue',
                         binwidth = 0.16, na.rm = T)+
            labs(title = title, x=x_label,
                 y = y_label)+
            theme(plot.title = element_text(color = "grey20", size = 30,
                                            hjust = 0.5),
                axis.text.x = element_text(color = "grey20", size = 30),
                  axis.text.y = element_text(color = "grey20", size = 30),
                  axis.title.x = element_text(color = "grey20", size = 30),
                  axis.title.y = element_text(color = "grey20", size = 30)))
        dev.off()
    }
    print(result)
}

# --- GWAS.violin.plots: 3-group violin plots for genotype vs taxon ---

GWAS.violin.plots <- function(abundance.file,
                               marker.genotype.file = NULL,
                               SNP = '2_136608646_G_A', violin.plot = T,
                               violin.plot.file = NULL,
                               bacterium = 'k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Bifidobacteriales;f__Bifidobacteriaceae;g__Bifidobacterium',
                               boundaries = list(0.5,1.5,NULL),
                               alleles = c('GG','AG','AA'),
                               y.axis.title = "Scaled bacterium abundance",
                               plot.type = 'clr', check.assoc.pval = T,
                              print.plot = T, color_ = "grey20", size_ = 30,
                              title_ = NULL, dotsize = 0.4, order_alleles = NULL,
                               ...
){
    if (is.null(marker.genotype.file)) stop('marker.genotype.file must be specified')
    if(is.null(order_alleles)){
        order_alleles <- 1:length(alleles)
    }
    genotype.abuncance.df <- construct.gt.taxon.file(GenotFile =
                                                         marker.genotype.file,
                                                     abundance.file = abundance.file,...)
    if(check.assoc.pval==T){
    for (i in 1:length(boundaries)){
         pval_lm_taxa_SNP <- test.lm.microb.gt(
        SNP.microb.df = genotype.abuncance.df, return.SNP.lm.coef = TRUE,
        recessive.model = boundaries[[i]], check.numeric = T)
    if (sum(rownames(pval_lm_taxa_SNP$SNP.coef)%in%SNP)!=1){
        stop ('SNP not found or found > 1 times in the data')
    }
    if (sum(colnames(pval_lm_taxa_SNP$SNP.coef)%in%bacterium)!=1){
        stop ('bacterium not found or found > 1 times in the data')
    }
    if (is.null(boundaries[[i]])){
        model = 'dosage model '
    } else if(boundaries[[i]]==0.5){
        model = 'recessive model '
    } else if (boundaries[[i]]==1.5){
        model = 'dominant model '
    } else {
        model = paste0(boundaries[[i]], ' dosage boundary model')
    }
        result <- paste0(taxonomic.unit(bacterium),' ~ ',SNP, ' regression coef: ',
    signif(pval_lm_taxa_SNP$SNP.coef[SNP,bacterium], digits = 3), ', ', model,
    'p-value: ', signif(pval_lm_taxa_SNP$pvalues[SNP,bacterium], digits = 3)
    )
        print(result)
    }
    }
    if(is.null(title_)){
        title_ <- taxonomic.unit(bacterium)
        }
    if(violin.plot == T){
        number.of.genotype.groups = length(unlist(boundaries))+1
        if(number.of.genotype.groups!=length(alleles)){
            warning('Number of alleles does not correspond to number of numerics
                    in boundaries! Allele captures are hense inaccurate')
        }
        if (sum(sort(unlist(boundaries))!=unlist(boundaries))>0){
            stop('boundaries must be a list with numeric values in increasing order.
                 They also must correspond to alleles')
        }
        if(!is.null(violin.plot.file)){
            jpeg(violin.plot.file, width = 1000, height = 600)
        }
        par(mar = c(12,6,6,2))
        colnames(genotype.abuncance.df)[colnames(genotype.abuncance.df) ==
                                             bacterium] <- 'the_bacterium'
        colnames(genotype.abuncance.df)[colnames(genotype.abuncance.df) ==
                                             SNP] <- 'SNP'
        genotype.abuncance.df <- genotype.abuncance.df[
            !is.na(genotype.abuncance.df$the_bacterium),]
        genotype.abuncance.df <- genotype.abuncance.df[
            !is.na(genotype.abuncance.df[,'SNP']),]
        genotype.abuncance.df_ = genotype.abuncance.df
        genotype.abuncance.df_[,'SNP'] <- alleles[1]
        for(j in seq_along(unlist(boundaries))){
            not.current.genotype.upwards <- genotype.abuncance.df[,'SNP']>=unlist(boundaries)[j]
            genotype.abuncance.df_[not.current.genotype.upwards,'SNP'] <- alleles[j+1]
        }
        genotype.abuncance.df_[,'SNP'] <- as.factor(genotype.abuncance.df_[,'SNP'])
        p <- ggplot(genotype.abuncance.df_, aes(
            x=SNP, y=the_bacterium)) + scale_x_discrete(limits=alleles[order_alleles]) +
            geom_violin(trim=FALSE, fill='cadetblue2', color="cadetblue3")
        if(plot.type == 'fraction'){
            final.plot <- (p + theme_bw() +
                      geom_boxplot(width=0.6, col = 'darkblue', fill = 'cadetblue2',
                                   na.rm = T) +
                      labs(title = title_, x="Genotype",
                           y = y.axis.title)+
                      theme(plot.title = element_text(color = color_, size = size_,
                                                      hjust = 0.5),
                            axis.text.x = element_text(color = color_, size = size_),
                            axis.text.y = element_text(color = color_, size = size_),
                            axis.title.x = element_text(color = color_, size = size_),
                            axis.title.y = element_text(color = color_, size = size_))+
                theme(plot.margin = unit(c(2,2,2,2), "cm")))
        } else if (plot.type == 'clr'){
            final.plot <- (p + theme_bw() +
                      geom_boxplot(width=0.6, col = 'darkblue', fill = 'cadetblue2',
                                   na.rm = T) +
                      geom_dotplot(binaxis='y', stackdir='center', dotsize=dotsize,
                                   col = 'cornflowerblue',
                                   binwidth = 0.16,
                                   na.rm = T)+
                      labs(title = title_, x="Genotype",
                           y = y.axis.title)+
                      theme(plot.title = element_text(color = color_, size = size_,
                                                      hjust = 0.5),
                            axis.text.x = element_text(color = color_, size = size_),
                            axis.text.y = element_text(color = color_, size = size_),
                            axis.title.x = element_text(color = color_, size = size_),
                            axis.title.y = element_text(color = color_, size = size_))+
                theme(plot.margin = unit(c(2,2,2,2), "cm")))
        } else{
            stop('plot type must be \'fractoion\' or \'clr\'')
        }
        if(print.plot==T){
            print(final.plot)
        } else{
            return(final.plot)
        }
            if(!is.null(violin.plot.file)){
            dev.off()
        }
    }
}

# --- get.pval.assoc.bact.marker: one-sided p-value for marker-bacterium association ---

get.pval.assoc.bact.marker <- function(abundance.file,
                           marker.genotype.file,
                           SNP, bacterium,
                           sign.assoc.exp,
                           boundary = NULL,
                           signif_ = 3
){
    genotype.abuncance.df <- construct.gt.taxon.file(GenotFile =
                                                         marker.genotype.file,
                                                     abundance.file = abundance.file)
    samples.nonNA <- Reduce('&',list((!is.na(genotype.abuncance.df[,SNP])),
                                      (!is.na(genotype.abuncance.df[,bacterium]))))
    nonzero.values.no <- sum(samples.nonNA)
    pval_lm_taxa_SNP <- test.lm.microb.gt(
        SNP.microb.df = genotype.abuncance.df, return.SNP.lm.coef = TRUE,
        recessive.model = boundary, check.numeric = T)
    if (sum(rownames(pval_lm_taxa_SNP$SNP.coef)%in%SNP)!=1){
        stop ('SNP not found or found > 1 times in the data')
    }
    if (sum(colnames(pval_lm_taxa_SNP$SNP.coef)%in%bacterium)!=1){
        stop ('bacterium not found or found > 1 times in the data')
    }
    if (is.null(boundary)){
        model = 'dosage model '
    } else if(boundary==0.5){
        model = 'recessive model '
    } else if (boundary==1.5){
        model = 'dominant model '
    } else {
        model = paste0(boundary, ' dosage boundary model')
    }
    if(sign.assoc.exp==sign(pval_lm_taxa_SNP$SNP.coef[SNP,bacterium])){
        p.val.oneSided <- signif(pval_lm_taxa_SNP$pvalues[SNP,bacterium],digits = signif_)
    }else{
        p.val.oneSided <- signif(1-(pval_lm_taxa_SNP$pvalues[SNP,bacterium]),digits = signif_)
    }
    result <- paste0(taxonomic.unit(bacterium),' ~ ',SNP, ' regression coef: ',
                     signif(pval_lm_taxa_SNP$SNP.coef[SNP,bacterium], digits = 3), ', ', model,
                     '2-sided p-value: ', signif(pval_lm_taxa_SNP$pvalues[SNP,bacterium], digits = 3),
                     ' (N = ',nonzero.values.no,')'
                     )
    print(result)
    return(list(p.val.oneSided,nonzero.values.no))
}

# --- prepare_Hardi: get genotype group counts ---

prepare_Hardi <- function(file, allele, thresholds = c(0.5,1.5),
                          samples.limit = NULL){
    genotypes <- readGenotFile(file, samples.limit = samples.limit)
    if(!is.null(samples.limit)){
        genotypes <- genotypes[rownames(genotypes)%in%samples.limit,]
    }
    allele.genotypes <- genotypes[,allele]
    if(!is.numeric(allele.genotypes)){
        allele.genotypes <- as.numeric(allele.genotypes)
    }
    allele.genotypes.RefRef.number <- length(allele.genotypes[allele.genotypes<thresholds[1]])
    allele.genotypes.RefAlt.number <- length(allele.genotypes[Reduce('&', list(
        (allele.genotypes>thresholds[1]),(allele.genotypes<thresholds[2])))])
    allele.genotypes.AltAlt.number <- length(allele.genotypes[allele.genotypes>thresholds[2]])
    return(c('Ref-Ref' = allele.genotypes.RefRef.number,
            'Ref-Alt' = allele.genotypes.RefAlt.number,
            'Alt-Alt' = allele.genotypes.AltAlt.number))
}

# ===========================================================================
# GWAS PLOTTING FUNCTIONS
# (Not used in core microbiome-transcriptome notebooks; included for
#  completeness. These require library(bedr) for bed2index/in.region.)
# ===========================================================================

# --- PLINK Manhattan plot ---

PLINK.plot <- function(file, jpeg.preffix, header_ = T, ylim_ = c(0, 10),
                       cex_ = 2, cex.lab_ = 3,
                       cex.axis_ = 3, lwd_ = 3, annotatePval_ = 0.000001,
                       main_ = T, file.first.N.elements = NULL,
                       file.suffix = '_PLINK_allSNPs.assoc.dosage',
                       output_dir = NULL,
                       CHR = NULL, BP = NULL,
                       ...){
    if (is.null(output_dir)) stop('output_dir must be specified')
    PLINK.result <- read.table(file, header = T, row.names = 1, sep = '', as.is = T,
                               check.names = F)
    if(!is.null(file.first.N.elements)){
        if(!is.numeric(file.first.N.elements)){
            stop('file.first.N.elements must be NULL or numeric')
        }
        PLINK.result <- PLINK.result[1:file.first.N.elements,]
    }
    if (header_){
        PLINK.result <- PLINK.result[-c(1),]
    }
    if((is.null(CHR)|!is.numeric(CHR))){
        CHR <- word(rownames(PLINK.result), end = 1, sep = '_')
        CHR <- as.numeric(CHR)
    }
    if((is.null(BP))|!is.numeric(BP)){
        BP <- word(rownames(PLINK.result), start = 2, end = 2, sep = '_')
        BP <- as.numeric(BP)
    }
    SNP <- rownames(PLINK.result)
    PLINK.result <- add_column(PLINK.result,BP,.before = 1)
    PLINK.result <- add_column(PLINK.result,CHR,.before = 1)
    PLINK.result <- add_column(PLINK.result,SNP,.before = 1)
    PLINK.result$P <- as.numeric(PLINK.result$P)
    if(main_){
        plot_title = word(taxonomic.unit(file),sep = file.suffix)
    } else{
        plot_title = NULL
    }
    jpg.name <- file.path(output_dir,
                          paste0(jpeg.preffix,'_',word(taxonomic.unit(file),sep = file.suffix),'.jpg')
    )
    jpeg(jpg.name, width = 2000, height = 1000)
    par(mar = c(9,9,5,18), mgp = c(5,2,0))
    manhattan.mod(PLINK.result, ylim = ylim_, cex = cex_,
                  cex.lab = cex.lab_, cex.axis = cex.axis_, lwd = lwd_,
                  annotatePval = annotatePval_, main = plot_title,...)
    dev.off()
}

# --- manhattan.mod: modified manhattan from qqman ---

manhattan.mod <- function (x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c("gray10",
                                                                                    "gray60"),
                           chrlabs = NULL, suggestiveline = -log10(1e-05),
                           genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE,
                           annotatePval = NULL, annotateTop = TRUE,
                           col_suggestive = "cadetblue2", col_genomewide = "darkblue",
                           ...)
{
    CHR = BP = P = index = NULL
    if (!(chr %in% names(x)))
        stop(paste("Column", chr, "not found!"))
    if (!(bp %in% names(x)))
        stop(paste("Column", bp, "not found!"))
    if (!(p %in% names(x)))
        stop(paste("Column", p, "not found!"))
    if (!(snp %in% names(x)))
        warning(paste("No SNP column found. OK unless you're trying to highlight."))
    if (!is.numeric(x[[chr]]))
        stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
    if (!is.numeric(x[[bp]]))
        stop(paste(bp, "column should be numeric."))
    if (!is.numeric(x[[p]]))
        stop(paste(p, "column should be numeric."))
    d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]])
    if (!is.null(x[[snp]]))
        d = transform(d, SNP = x[[snp]])
    d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
    d <- d[order(d$CHR, d$BP), ]
    if (logp) {
        d$logp <- -log10(d$P)
    }
    else {
        d$logp <- d$P
    }
    d$pos = NA
    d$index = NA
    ind = 0
    for (i in unique(d$CHR)) {
        ind = ind + 1
        d[d$CHR == i, ]$index = ind
    }
    nchr = length(unique(d$CHR))
    if (nchr == 1) {
        d$pos = d$BP
        ticks = floor(length(d$pos))/2 + 1
        xlabel = paste("Chromosome", unique(d$CHR), "position")
        labs = ticks
    }
    else {
        lastbase = 0
        ticks = NULL
        for (i in unique(d$index)) {
            if (i == 1) {
                d[d$index == i, ]$pos = d[d$index == i, ]$BP
            }
            else {
                lastbase = lastbase + tail(subset(d, index ==
                                                      i - 1)$BP, 1)
                d[d$index == i, ]$pos = d[d$index == i, ]$BP +
                    lastbase
            }
            ticks = c(ticks, (min(d[d$index == i, ]$pos) + max(d[d$index ==
                                                                     i, ]$pos))/2 + 1)
        }
        xlabel = "Chromosome"
        labs <- unique(d$CHR)
    }
    xmax = ceiling(max(d$pos) * 1.03)
    xmin = floor(max(d$pos) * -0.03)
    def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i",
                     las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(0,
                                                                       ceiling(max(d$logp))), xlab = xlabel, ylab = expression(-log[10](italic(p))))
    dotargs <- list(...)
    do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in%
                                                names(dotargs)]))
    if (!is.null(chrlabs)) {
        if (is.character(chrlabs)) {
            if (length(chrlabs) == length(labs)) {
                labs <- chrlabs
            }
            else {
                warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
            }
        }
        else {
            warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
        }
    }
    if (nchr == 1) {
        axis(1, ...)
    }
    else {
        axis(1, at = ticks, labels = labs, ...)
    }
    col = rep(col, max(d$CHR))
    if (nchr == 1) {
        with(d, points(pos, logp, pch = 20, col = col[1], ...))
    }
    else {
        icol = 1
        for (i in unique(d$index)) {
            with(d[d$index == unique(d$index)[i], ], points(pos,
                                                            logp, col = col[icol], pch = 20, ...))
            icol = icol + 1
        }
    }
    if (suggestiveline)
        abline(h = suggestiveline, col = col_suggestive, lwd = 2)
    if (genomewideline)
        abline(h = genomewideline, col = col_genomewide, lwd = 2)
    if (!is.null(highlight)) {
        if (any(!(highlight %in% d$SNP)))
            warning("You're trying to highlight SNPs that don't exist in your results.")
        d.highlight = d[which(d$SNP %in% highlight), ]
        with(d.highlight, points(pos, logp, col = "green3",
                                 pch = 20, ...))
    }
    if (!is.null(annotatePval)) {
        topHits = subset(d, P <= annotatePval)
        par(xpd = TRUE)
        if (annotateTop == FALSE) {
            with(subset(d, P <= annotatePval), textxy(pos, -log10(P),
                                                      offset = 0.625, labs = topHits$SNP),
                 ...)
        }
        else {
            topHits <- topHits[order(topHits$P), ]
            topSNPs <- NULL
            for (i in unique(topHits$CHR)) {
                chrSNPs <- topHits[topHits$CHR == i, ]
                topSNPs <- rbind(topSNPs, chrSNPs[1, ])
            }
            textxy(topSNPs$pos, -log10(topSNPs$P), offset = 0.625,
                   labs = topSNPs$SNP, ...)
        }
    }
    par(xpd = FALSE)
}

# --- plot.region.of.interest.GWAS ---
# Requires: library(bedr) for bed2index() and in.region()

plot.region.of.interest.GWAS <- function(GWAS.file, gene.region.file = NULL,
                                         exon.file = NULL,
                                         plot.file, original.SNP.coord = NULL,
                                         highlight.file = NULL, resolution = 1,
                                         original_SNP_highlight_diameter = 1,
                                         labels_pval_threshold = 1e-12,
                                         main.scale.factor = 1, ylim_ = NULL,
                                         ...){
    if (!requireNamespace('bedr', quietly = TRUE)) {
        stop('plot.region.of.interest.GWAS requires the bedr package')
    }
    Bacterium_50k <- read.table(GWAS.file,
                                sep = '', header = T, row.names = 1)
    Bacterium_50k <- Bacterium_50k[!is.na(Bacterium_50k$P),]
    coordinate <- word(rownames(Bacterium_50k), start = 2, end = 2, sep = '_')
    chromosome <- word(rownames(Bacterium_50k), start = 1, end = 1, sep = '_')[1]
    if(length(unique(word(rownames(Bacterium_50k), start = 1, end = 1, sep = '_')))!=1){
        warning ('>1 chromosome in GWAS.file, should be only one')
    }
    Bacterium.df <- cbind(as.numeric(coordinate), -log10(Bacterium_50k$P))
    SNPs.indexed <- bedr::bed2index(cbind(paste0('chr',word(rownames(Bacterium_50k),
                                                      end = 1, sep = '_')),
                                    Bacterium.df[,1],Bacterium.df[,1]))
    if(!is.null(exon.file)){
        exons.bacterium.region <- read.table(exon.file, sep = '\t')
        exon.bed.indexed <- bedr::bed2index(exons.bacterium.region)
        coding.SNPs <- bedr::in.region(SNPs.indexed,exon.bed.indexed)
        max.ass.coding <- Bacterium_50k[coding.SNPs,][which.max(Bacterium.df[coding.SNPs,2]),]
        print(paste0('Maximum coding SNP associated with bacterium: ',rownames(max.ass.coding),
                     ', p = ', max.ass.coding$P))
    } else{
        coding.SNPs <- NULL
    }
    if(!is.null(gene.region.file)){
        genes.bacterium.region <- read.table(gene.region.file, sep = '\t')
        genes.bed.indexed <- bedr::bed2index(genes.bacterium.region)
        gene.region.SNPs <- bedr::in.region(SNPs.indexed,genes.bed.indexed)
        max.ass.gene <- Bacterium_50k[gene.region.SNPs,][which.max(Bacterium.df[gene.region.SNPs,2]),]
        print(paste0('Maximum gene-region SNP associated with bacterium: ',rownames(max.ass.gene),
                     ', p = ', max.ass.gene$P))
    } else{
        gene.region.SNPs <- NULL
    }
    max.ass <- Bacterium_50k[which.max(Bacterium.df[,2]),]
    print(paste0('Maximum SNP associated with bacterium: ',rownames(max.ass),
             ', p = ', max.ass$P))

    jpeg(plot.file, width = 800*resolution, height = 600*resolution)
    par(mar = c(11+2*resolution,9,9,1), mgp = c(4,1,0))
    plot(Bacterium.df, pch = 20, xlab = '',
         cex = 2*resolution, cex.axis = 2*resolution,
         type="n", xaxt="n", yaxt="n",
         cex.main = 2*resolution, cex.lab = 2*resolution,
         col = 'grey',...)
    my.legend.size <- legend("topleft", c(
        "Region of gene", "Exon +- 5bp", "orignal SNP", "All SNPs"),col=c(
            "black", "red","darkgreen","grey"),
        pch = c(1,20,20,20), cex=c(1*resolution,
                                   1*resolution,
                                   1*resolution,
                                   1*resolution), plot = FALSE
    )
    my.range = range(Bacterium.df[,2])
    my.range[2] <- 1.04*(my.range[2]+my.legend.size$rect$h[1])
    if(is.null(ylim_)){
        ylim_ = my.range
    }
    plot(Bacterium.df, pch = 20, xaxt = 'n',xlab = '',
         ylab = expression('-log'[10]*'(association p-value)'),
         cex = 2*resolution, cex.axis = 2*resolution,
         cex.main = 2*resolution*main.scale.factor, cex.lab = 2*resolution,
         col = 'grey',ylim = ylim_,...)
    if(sum(Bacterium.df[,2]>=-log10(labels_pval_threshold))>0){
        text(x = Bacterium.df[Bacterium.df[,2]>=-log10(labels_pval_threshold),1],
             y = Bacterium.df[Bacterium.df[,2]>=-log10(labels_pval_threshold),2],
             labels = rownames(Bacterium_50k)[Bacterium.df[,2]>=-log10(labels_pval_threshold)],
             cex = 0.7*resolution, adj = c(0,0))
    }
    legend("topleft", legend=c(
        "Region of gene", "Exon +- 5bp", "orignal SNP", "All SNPs"),col=c(
            "black", "red","darkgreen","grey"),
        pch = c(1,20,20,20), cex=c(1*resolution,
                                   1*resolution,
                                   1*resolution,
                                   1*resolution)
    )
    if(!is.null(highlight.file)){
        SNPs2highlight.df <- read.table(highlight.file, sep = '\t', header = F)
        current.chromosome <- unique(word(rownames(Bacterium_50k)
                                          [2:nrow(Bacterium_50k)], end = 1, sep = '_')
                                     )
        current.chromosome <- as.numeric(current.chromosome)
        SNPs2highlight.df$V1 <- as.numeric(SNPs2highlight.df$V1)
        if (length(current.chromosome)!=1){
            warning('All SNPs in highlight.file should be on 1 chromosome! Highlighting will not be performed.')
            SNPs2highlight.df <- NULL
        }
        SNPs2highlight.df <- SNPs2highlight.df[
            SNPs2highlight.df$V1==current.chromosome,]
        points(Bacterium.df[Bacterium.df[,1]%in%SNPs2highlight.df$V4,],
               col = 'bisque4', pch = 20, cex = 2*resolution)
    }
    points(Bacterium.df[gene.region.SNPs,], col = 'black', pch = 1,
           cex = 2*resolution)
    points(Bacterium.df[coding.SNPs,], col = colors()[556], pch = 20,
           cex = 1*resolution)
    if(!is.null(original.SNP.coord)){
        Bacterium.df.orig.SNP <- Bacterium.df[Bacterium.df[,1]%in%original.SNP.coord,]
        if(is.vector(Bacterium.df.orig.SNP)){
            Bacterium.df.orig.SNP <- t(Bacterium.df.orig.SNP)
        }
        points(Bacterium.df.orig.SNP, col = 'darkgreen',
               pch = 20, cex = original_SNP_highlight_diameter*resolution)
    }
    t <- nrow(Bacterium.df)
    x.axis.captures <- c(
        1, round(t), round(t/2), round(t/4), round(3*t/4))
    axis(1, las = 2, labels = round(Bacterium.df[x.axis.captures,1]/1e6, digits = 2),
         at = Bacterium.df[x.axis.captures,1], cex.axis = 2*resolution)
    title(xlab=paste0("Coordinate on chr ",chromosome,", Mb"), line=9+2*resolution,
          cex.lab = 2*resolution)
      dev.off()
}

# --- plot.mGWAS.eQTL.together ---

plot.mGWAS.eQTL.together <- function(GWAS.file, eQTL.file,
                                     plot.file, original.SNP.coord = NULL,
                                     resolution = 1, original_SNP_highlight_diameter = 1,
                                     labels_pval_threshold = 1e-12,
                                     main.scale.factor = 1, ylim_ = NULL,
                                     xlim_ = NULL,
                                     pch_ = c(20,20), color_ = c('black','darkred'),
                                     sub_ = NULL, cex.sub_ = NULL, line.sub_ = NULL,
                                     col.sub_ = NULL,
                                     ...){
    GWAS_pval <- read.table(GWAS.file,
                            sep = '', header = T, row.names = 1)
    GWAS_pval <- GWAS_pval[!is.na(GWAS_pval$P),]
    coordinate <- word(rownames(GWAS_pval), start = 2, end = 2, sep = '_')
    chromosome <- word(rownames(GWAS_pval), start = 1, end = 1, sep = '_')[1]
    if(length(unique(word(rownames(GWAS_pval), start = 1, end = 1, sep = '_')))!=1){
        warning ('>1 chromosome in GWAS.file, should be only one')
    }
    GWAS.df <- cbind(as.numeric(coordinate), -log10(GWAS_pval$P))
    eQTL_pval <- read.table(eQTL.file, sep = '', header = T, row.names = 1)
    eQTL_pval <- eQTL_pval[!is.na(eQTL_pval$P),]
    coordinate <- word(rownames(eQTL_pval), start = 2, end = 2, sep = '_')
    chromosome <- word(rownames(eQTL_pval), start = 1, end = 1, sep = '_')[1]
    if(length(unique(word(rownames(eQTL_pval), start = 1, end = 1, sep = '_')))!=1){
        warning ('>1 chromosome in GWAS.file.2, should be only one')
    }
    eQTL.df <- cbind(as.numeric(coordinate), -log10(eQTL_pval$P))
    my.range.x = range(c(GWAS.df[,1],eQTL.df[,1]))
    if(is.null(xlim_)){
        xlim_ = my.range.x
    }
    GWAS.df = GWAS.df[Reduce('&',list((GWAS.df[,1]>=xlim_[1]),(GWAS.df[,1]<=xlim_[2]))),
                      ]
    eQTL.df = eQTL.df[Reduce('&',list((eQTL.df[,1]>=xlim_[1]),(eQTL.df[,1]<=xlim_[2]))),
                      ]
    jpeg(plot.file, width = 800*resolution, height = 600*resolution)
    par(mar = c(11+2*resolution,9,9,1), mgp = c(4,1,0))
    plot(GWAS.df, pch = pch_[1], xlab = '',
         cex = 2*resolution, cex.axis = 2*resolution,
         type="n", xaxt="n", yaxt="n",
         cex.lab = 2*resolution,
         col = 'grey',xlim = xlim_)
    my.legend.size <- legend("topleft", c(
        "Microbiome GWAS","eQTL","original SNPs"),col=c(color_, 'darkgreen'),
        pch = c(pch_,20), cex=c(1*resolution,
                                1*resolution,
                                original_SNP_highlight_diameter*resolution),
        plot = FALSE
    )
    my.range.y = range(c(GWAS.df[,2],eQTL.df[,2]))
    my.range.y[2] <- 1.04*(my.range.y[2]+my.legend.size$rect$h[1])
    if(is.null(ylim_)){
        ylim_ = my.range.y
    }
    plot(GWAS.df, pch = pch_[1], xaxt = 'n',xlab = '',
         ylab = expression('-log'[10]*'(association p-value)'),
         cex = 2*resolution, cex.axis = 2*resolution,
         cex.main = 2*resolution*main.scale.factor, cex.lab = 2*resolution,
         col = color_[1],ylim = ylim_,xlim = xlim_, ...)
    points(x = eQTL.df[,1],y = eQTL.df[,2], pch = pch_[2],
           col = color_[2], cex = 2*resolution)
    if(sum(GWAS.df[,2]>=-log10(labels_pval_threshold[1]))>0){
        text(x = GWAS.df[GWAS.df[,2]>=-log10(labels_pval_threshold),1],
             y = GWAS.df[GWAS.df[,2]>=-log10(labels_pval_threshold),2],
             labels = rownames(GWAS_pval)[GWAS.df[,2]>=-log10(labels_pval_threshold)],
             cex = 0.7*resolution, adj = c(0,0))
    }
    legend("topleft", c(
        "Microbiome GWAS","eQTL","mGWAS hit SNP"),col=c(color_, 'darkgreen'),
        pch = c(pch_,20), cex=c(1*resolution,
                                1*resolution,
                                original_SNP_highlight_diameter*resolution)
    )
    if(!is.null(sub_)){
        if(is.null(cex.sub_)){
            cex.sub_ =4
        }
        if(is.null(line.sub_)){
            line.sub_ =-30
        }
        if(is.null(col.sub_)){
            col.sub_ ='darkblue'
        }
        title(sub = sub_, cex.sub = cex.sub_, line = line.sub_, col.sub = col.sub_)
    }
    if(!is.null(original.SNP.coord)){
        GWAS.df.orig.SNP <- GWAS.df[GWAS.df[,1]%in%original.SNP.coord,]
        if(is.vector(GWAS.df.orig.SNP)){
            GWAS.df.orig.SNP <- t(GWAS.df.orig.SNP)
        }
        points(GWAS.df.orig.SNP, col = 'darkgreen',
               pch = 20, cex = original_SNP_highlight_diameter*resolution)
        eQTL.df.orig.SNP <- eQTL.df[eQTL.df[,1]%in%original.SNP.coord,]
        if(is.vector(eQTL.df.orig.SNP)){
            eQTL.df.orig.SNP <- t(eQTL.df.orig.SNP)
        }
        points(eQTL.df.orig.SNP, col = 'darkgreen',
               pch = 20, cex = original_SNP_highlight_diameter*resolution)
    }
    t <- nrow(GWAS.df)
    x.axis.captures <- c(
        1, round(t), round(t/2), round(t/4), round(3*t/4))
    axis(1, las = 2, labels = round(GWAS.df[x.axis.captures,1]/1e6, digits = 2),
         at = GWAS.df[x.axis.captures,1], cex.axis = 2*resolution)
    title(xlab=paste0("Coordinate on chr ",chromosome,", Mb"), line=9+2*resolution,
          cex.lab = 2*resolution)
    dev.off()
}

# Note: plot.regional.eQTL.GTEx() is omitted from this adapted version because
# it depends on liftOver binary and hardcoded temporary file paths on the
# original server. It can be restored from Microbiota_QTL_functions_v9.R
# lines 1281-1414 if needed.


