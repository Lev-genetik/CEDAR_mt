# functions_general.R
# Adapted from: code/CEDAR-master/Functions_general_non_enterotyping_v3m.R
# General utility functions used across multiple analysis notebooks.

require(fitdistrplus)
require(dplyr)
require(data.table)
require(stringr)

# --- I/O utilities ---

read.table.smart <- function(input_file, ...) {
    read.table(input_file, header = T, row.names = 1,
               sep = '\t', as.is = T, check.names = F, ...)
}

write.table.smart <- function(data, output_file, sep = "\t", quote = F,
                              col.names = NA, ...) {
    write.table(data, file = output_file, sep = sep, quote = quote,
                col.names = col.names, ...)
}

# --- Merging utilities ---

merge.by.rownames <- function(df1, df2, stringsAsFactors_ = FALSE) {
    result <- merge(df1, df2, by = 'row.names', all = T,
                    stringsAsFactors = stringsAsFactors_)
    rownames(result) <- result[, 1]
    result <- as.data.frame(result[, -1])
}

merge.by.rownames.multiple <- function(df_list, stringsAsFactors_ = FALSE) {
    result.df <- NULL
    for (i in seq_along(df_list)) {
        result.df <- merge.by.rownames(result.df, df_list[[i]])
    }
    result.df
}

merge.by.rownames.flex <- function(df1, df2, stringsAsFactors_ = FALSE,
                                   by = 'row.names', ...) {
    result <- merge(df1, df2, by = by,
                    stringsAsFactors = stringsAsFactors_, ...)
    rownames(result) <- result[, 1]
    result <- as.data.frame(result[, -1])
}

merge.by.rownames.flex.multiple <- function(df_list, stringsAsFactors_ = FALSE,
                                            by = 'row.names', ...) {
    if (length(df_list) < 2) {
        warning('Size of df_list <2')
        return(df_list)
    }
    result.df <- df_list[[1]]
    for (i in 2:length(df_list)) {
        result.df <- merge.by.rownames.flex(result.df, df_list[[i]],
                                            stringsAsFactors_ = FALSE,
                                            by = by, ...)
    }
    return(result.df)
}

merge.by.colnames <- function(df1, df2, stringsAsFactors_ = FALSE) {
    as.data.frame(t(merge.by.rownames(t(df1), t(df2))),
                  stringsAsFactors = stringsAsFactors_)
}

# --- Abundance / preprocessing ---

to.scale <- function(data) {
    data_col_sum <- colSums(data)
    data_norm <- data
    for (i in 1:dim(data)[2]) {
        data_norm[, i] <- data[, i] / data_col_sum[i]
    }
    return(data_norm)
}

# Replace 0 counts by adding pseudocount =
# minimum observed proportion for that taxon * sampleSum for that sample
impute.abundance <- function(df, df_scaled, digits = 2) {
    column_sums <- colSums(df)
    min.fraction <- numeric(nrow(df_scaled))
    for (i in seq_along(rownames(df_scaled))) {
        min.fraction[i] <- min(df_scaled[i, ][df_scaled[i, ] > 0])
        if (min.fraction[i] == Inf) {
            min.fraction[i] <- 0.000001
        }
    }
    for (i in seq_along(rownames(df))) {
        for (j in seq_along(colnames(df))) {
            if (df[i, j] == 0) {
                df[i, j] <- round((min.fraction[i]) * column_sums[j], digits = digits)
            }
        }
    }
    print('Imputation done.')
    df
}

# --- Taxonomy ---

# Get the most precise available taxonomic unit for a taxon
taxonomic.unit <- function(name, sep = ';') {
    start.taxonomic.level <- sapply(strsplit(name, sep), length)
    for (j in start.taxonomic.level:1) {
        taxonomic.levels <- word(name, start = j, end = j, sep = sep)
        if (nchar(taxonomic.levels) > 4) {
            break
        }
    }
    taxonomic.levels
}

refine.taxa.vector <- function(taxa.vector, keep.only = 'g__') {
    taxa.vector <- taxa.vector[taxa.vector %like% keep.only]
    taxa.vector <- taxa.vector[!is.na(taxa.vector)]
    taxa.vector <- taxa.vector[nchar(taxa.vector) > 3]
    taxa.vector <- gsub('g__', '', taxa.vector)
    taxa.vector <- taxa.vector[!taxa.vector %in% c(
        'uncultured', 'metagenome', 'uncultured.bacterium', 'gut.metagenome')]
    return(taxa.vector)
}

# --- Statistical utilities ---

# Correct p-values for genomic inflation (chi-squared lambda method)
correct4inflation <- function(pvalues) {
    p_values_khisq <- qchisq(pvalues, df = 1, low = F)
    p_values_khisq_med <- median(p_values_khisq)
    lambda <- p_values_khisq_med / qchisq(0.5, 1, low = F)
    p_values_khisq_corr <- p_values_khisq / lambda
    p_values_corrected <- pchisq(p_values_khisq_corr, df = 1, low = F)
    return(p_values_corrected)
}

# Remove outliers from a vector: values beyond Q1/Q3 +/- IQR*threshold -> NA
remove.outliers <- function(x, remove.outliers.IQR = 1.5) {
    lh <- quantile(x, probs = 0.25, na.rm = T)
    uh <- quantile(x, probs = 0.75, na.rm = T)
    step <- remove.outliers.IQR * (uh - lh)
    x[x > uh + step] <- NA
    x[x < lh - step] <- NA
    return(x)
}

# JSD distance matrix
dist.JSD.df <- function(inMatrix, pseudocount = 0.000001, ...) {
    KLD <- function(x, y) sum(x * log(x / y))
    JSD <- function(x, y) sqrt(0.5 * KLD(x, (x + y) / 2) + 0.5 * KLD(y, (x + y) / 2))
    matrixColSize <- length(colnames(inMatrix))
    colnames <- colnames(inMatrix)
    resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
    inMatrix <- apply(inMatrix, 1:2, function(x) ifelse(x == 0, pseudocount, x))
    for (i in 1:matrixColSize) {
        for (j in 1:matrixColSize) {
            resultsMatrix[i, j] <- JSD(as.vector(inMatrix[, i]),
                                       as.vector(inMatrix[, j]))
        }
    }
    colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
    resultsMatrix <- as.data.frame(resultsMatrix)
    return(resultsMatrix)
}

# Average abundance across locations
get.average.abundance.loc <- function(input_df, cols.to.corr = 8,
                                      number_of_letters_rownames = 6,
                                      round_ = 1) {
    samples <- unique(str_sub(rownames(input_df), end = number_of_letters_rownames))
    result <- matrix(NA, ncol = cols.to.corr, nrow = length(samples))
    colnames(result) <- colnames(input_df)[1:cols.to.corr]
    rownames(result) <- samples
    black.labels <- logical(nrow(input_df))
    for (sample in samples) {
        rownumbers <- which(rownames(input_df) %like% sample)
        for (i in 1:cols.to.corr) {
            result[which(samples %in% sample), i] <-
                mean(as.numeric(input_df[c(rownumbers), i]), na.rm = TRUE)
        }
        if (length(rownumbers) > 1) {
            black.labels[rownumbers[2:length(rownumbers)]] <- T
        }
    }
    if (cols.to.corr < ncol(input_df)) {
        result_ <- cbind(result, input_df[!black.labels, (cols.to.corr + 1):ncol(input_df)])
    } else {
        result_ <- result
    }
    if (!is.null(round_)) {
        result_ <- round(result_, digits = round_)
    }
    return(result_)
}

# Weighted z-score combination of p-values
get.combined.pval <- function(pval1.right, pval2.right, N1, N2) {
    z_combined <- (qnorm(1 - pval1.right) * sqrt(N1) + qnorm(1 - pval2.right) * sqrt(N2)) / sqrt(N1 + N2)
    return(1 - pnorm(z_combined))
}

# Read QIIME output with filtering, outlier removal, and optional scaling
# (Full function from v3m — included for completeness, not needed for core analysis
#  since corrected tables are provided)
readQiimeSmart <- function(input_file, locat = 'all', percent_threshold = 0.01,
                           remove.outliers = 0, sample.cov.filter = 0, to.scale = F,
                           rare.taxa.advocates = NULL, advocate.coverage = c(10000, 10000),
                           cov.filter.at.end.of.cleaning = F) {
    data <- read.table(input_file, header = T, row.names = 1, dec = ".", sep = ",",
                       check.names = F) %>%
        dplyr::select(-Location, -MiSeq_sample_number, -MiSeq_lane, -MiSeq_run_date_YYMMDD)

    if (!is.null(rare.taxa.advocates)) {
        if (length(rare.taxa.advocates) != 2)
            stop('rare.taxa.advocates must be a character vector of 2 files')
        if (length(advocate.coverage) != 2)
            stop('advocate.coverage must be a numeric vector of 2 values')
        advocate1 <- read.table(rare.taxa.advocates[1], header = T, row.names = 1, dec = ".", sep = ",")
        advocate2 <- read.table(rare.taxa.advocates[2], header = T, row.names = 1, dec = ".", sep = ",")
        advocate1 <- advocate1[, colSums(advocate1) >= advocate.coverage[1]]
        advocate2 <- advocate2[, colSums(advocate2) >= advocate.coverage[2]]
    }
    if (locat != 'all') {
        data <- dplyr::select(data, contains(locat))
        if (!is.null(rare.taxa.advocates)) {
            advocate1 <- dplyr::select(advocate1, contains(locat))
            advocate2 <- dplyr::select(advocate2, contains(locat))
        }
    }
    if (!cov.filter.at.end.of.cleaning) {
        data <- data[, colSums(data) >= sample.cov.filter]
    }
    noise.removal <- function(dataframe, top = NULL) {
        dataframe -> Matrix
        bigones.data <- rowSums(Matrix) * 100 / (sum(rowSums(Matrix))) > percent_threshold
        bigones <- bigones.data
        if (!is.null(rare.taxa.advocates)) {
            bigones.advocate1 <- rowSums(advocate1) * 100 / (sum(rowSums(advocate1))) > percent_threshold
            bigones.advocate2 <- rowSums(advocate2) * 100 / (sum(rowSums(advocate2))) > percent_threshold
            for (taxon in names(bigones)) {
                bigones[taxon] <- as.logical(max(bigones.data[taxon],
                                                 bigones.advocate1[taxon],
                                                 bigones.advocate2[taxon], na.rm = T))
            }
        }
        Matrix_1 <- Matrix[bigones, ]
        return(Matrix_1)
    }
    data <- noise.removal(data)
    if (remove.outliers != 0) {
        dist.JSD <- function(inMatrix, pseudocount = 0.000001, ...) {
            KLD <- function(x, y) sum(x * log(x / y))
            JSD <- function(x, y) sqrt(0.5 * KLD(x, (x + y) / 2) + 0.5 * KLD(y, (x + y) / 2))
            matrixColSize <- length(colnames(inMatrix))
            colnames <- colnames(inMatrix)
            resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
            inMatrix <- apply(inMatrix, 1:2, function(x) ifelse(x == 0, pseudocount, x))
            for (i in 1:matrixColSize) {
                for (j in 1:matrixColSize) {
                    resultsMatrix[i, j] <- JSD(as.vector(inMatrix[, i]),
                                               as.vector(inMatrix[, j]))
                }
            }
            colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
            as.dist(resultsMatrix) -> resultsMatrix
            attr(resultsMatrix, "method") <- "dist"
            return(resultsMatrix)
        }
        data.dist <- dist.JSD(data)
        dist.m <- as.matrix(data.dist)
        n <- round(nrow(dist.m) / 2)
        meds.50 <- apply(dist.m, 1, function(x) {
            x.sort <- sort(x, decreasing = F)[-1]
            median(x.sort[1:n])
        })
        f1 <- fitdist(meds.50, "norm")
        p.vals <- pnorm(meds.50, mean = f1$estimate[1], sd = f1$estimate[2],
                        lower.tail = F, log.p = F)
        drops.new <- names(p.vals)[p.vals < remove.outliers]
        dist.m.new <- dist.m[!rownames(dist.m) %in% drops.new, !colnames(dist.m) %in% drops.new]
        data.dist <- as.dist(dist.m.new)
        attr(data.dist, "method") <- "dist"
        data <- data[, !colnames(data) %in% drops.new]
    }
    if (cov.filter.at.end.of.cleaning) {
        data <- data[, colSums(data) >= sample.cov.filter]
    }
    if (to.scale) {
        data_col_sum <- colSums(data)
        data_norm <- data
        for (i in 1:dim(data)[2]) {
            data_norm[, i] <- data[, i] / data_col_sum[i]
        }
        data <- data_norm
    }
    return(data)
}
