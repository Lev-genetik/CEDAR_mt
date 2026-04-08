# functions_figures.R
# Adapted from: code/CEDAR-master/Scripts/Figures_for_paper_mt/constants.R
#               code/data_for_plots/general_func.R
# Constants and plotting utilities for manuscript figures.
# Note: the two copies of constants.R and general_func.R are identical.

require(data.table)
require(lattice)

# --- Constants ---

loc.full <- list(RE = 'rectum',
                 IL = 'ileum',
                 TR = 'transverse colon')

immune.cells.types <- list(CD4 = 'CD4',
                           CD8 = 'CD8',
                           CD14 = 'CD14',
                           CD15 = 'CD15',
                           CD19 = 'CD19',
                           PLA = 'PLA')

loc.plus.immune <- c(loc.full, immune.cells.types)

regions <- c('V1V2', 'V3V4', 'V5V6')

# Color palette: named vector for all 9 tissue/cell types
colors <- c('darkgoldenrod1', 'darkcyan', 'deeppink4', 'deepskyblue', 'orangered1',
            'palevioletred2', 'chartreuse4', 'firebrick4', 'black')
names(colors) <- c('IL', 'TR', 'RE', 'CD4', 'CD8', 'CD14', 'CD15', 'CD19', 'PLA')

# Transparent color variants (alpha = 50)
.make_transparent_colors <- function(colors, alpha_val) {
    colors.rgb.dt <- as.data.table(col2rgb(colors))
    result <- lapply(names(colors), function(nm) {
        rgb(colors.rgb.dt[1, nm, with = FALSE][[1]],
            colors.rgb.dt[2, nm, with = FALSE][[1]],
            colors.rgb.dt[3, nm, with = FALSE][[1]],
            alpha = alpha_val, maxColorValue = 255)
    })
    names(result) <- names(colors)
    return(result)
}

colors.transparent <- .make_transparent_colors(colors, 50)
colors.transparent.80 <- .make_transparent_colors(colors, 80)

p.val.thresh <- 0.05

# --- Plotting functions ---

# Simple QQ plot for p-values
p.qqplot <- function(pvector, main = NULL, ...) {
    par(pty = "s")
    o <- -log10(sort(pvector, decreasing = F))
    e <- -log10(1:length(o) / length(o))
    limits <- c(0, max(c(e, o)))
    plot(x = e, y = o, pch = 19, cex = 1, main = main, ...,
         xlab = expression(Expected ~~ -log[10](italic(p))),
         ylab = expression(Observed ~~ -log[10](italic(p))),
         xlim = limits, ylim = limits, col = 'cadetblue')
    lines(e, e, col = "red")
}

# QQ plot with optional pretty scatter mode
qq <- function(dt, main, loc, colors.transparent, num.dots.to.label = 2,
               mode = 'simple', cex.lab = 1.5) {
    setorder(dt, pval)
    dt$logp <- -log10(dt$pval)
    o <- dt$logp
    e <- -log10(((1:length(o))) / length(o))
    dt$e <- e
    par(pty = "s", pch = 20, cex = 0.7)
    if (mode == 'simple') {
        plot(e, o,
             xlab = expression(Expected ~~ -log[10](italic(p))),
             ylab = expression(Observed ~~ -log[10](italic(p))),
             xlim = c(0, max(e)),
             ylim = c(0, max(o)),
             main = main,
             cex.lab = 1.5,
             col = 'royalblue')
    } else if (mode == 'pretty') {
        PrettyScatter(x = e, y = o,
                      xlab = expression(Expected ~~ -log[10](italic(p))),
                      ylab = expression(Observed ~~ -log[10](italic(p))),
                      main = main,
                      bg = colors.transparent[[loc]])
    }
    labels <- paste0(dt[1:num.dots.to.label, ]$bPC, '&', dt[1:num.dots.to.label, ]$gPC)
    lines(e, e, col = "black")
    return(dt)
}

PrettyScatter <- function(x, y, main, bg, panel.first.step = 1, cex.lab = 1.5,
                          xlim = F, ylim = F, xlab = '', ylab = '', abline = F) {
    if (identical(xlim, FALSE)) xlim <- c(min(x), max(x))
    if (identical(ylim, FALSE)) ylim <- c(min(y), max(y))
    plot(x, y, xlim = xlim, ylim = ylim,
         xlab = xlab, ylab = ylab,
         cex.lab = cex.lab,
         bg = bg, col = bg, pch = 20,
         axes = FALSE, frame.plot = FALSE, main = main)
    if (abline == T) {
        abline(h = seq(round(min(x)), round(max(y)), panel.first.step), col = 'grey60')
    }
    at <- pretty(x)
    at <- at[-c(1, length(at))]
    mtext(side = 1, text = at, at = at, col = "grey20", line = 1, cex = 0.9)
    at <- pretty(y)
    at <- at[-c(1, length(at))]
    mtext(side = 2, text = at, at = at, col = "grey20", line = 1, cex = 0.9)
}

# Variance explained helpers
SumSubvector <- function(vector, num.elements) {
    sum(vector[1:num.elements])
}

FindNumberOfPCsExplainingVarianceThreshold <- function(variance, variance.sum.thresh) {
    subsum <- SumSubvector(variance, num.elements = 1)
    num.elements <- 1
    while (subsum < variance.sum.thresh) {
        num.elements <- num.elements + 1
        subsum <- SumSubvector(variance, num.elements)
    }
    return(num.elements)
}

# Outlier detection (IQR-based, from general_func.R)
GetOutliers <- function(dt, feature) {
    dt <- dt[, c('id', feature), with = F]
    x <- as.numeric(unlist(dt[, feature, with = F]))
    names(x) <- dt$id
    bottom.thresh <- quantile(x)['25%'] - 5 * iqr(x)
    up.thresh <- quantile(x)['75%'] + 5 * iqr(x)
    names(x[x < bottom.thresh | x > up.thresh])
}

# Recursive merge
MergeRecurse <- function(list.dts, by.col) {
    Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = by.col, all.x = TRUE, all.y = TRUE),
           list.dts)
}

# Lattice-based QQ plot for uniform p-values (from general_func.R)
qqunif.plot <- function(pvalues,
                        should.thin = T, thin.obs.places = 2, thin.exp.places = 2,
                        xlab = expression(paste("Expected (", -log[10], " p-value)")),
                        ylab = expression(paste("Observed (", -log[10], " p-value)")),
                        draw.conf = TRUE, conf.points = 1000, conf.col = "lightgray", conf.alpha = .05,
                        already.transformed = FALSE, pch = 20, aspect = "iso",
                        prepanel = prepanel.qqunif,
                        par.settings = list(superpose.symbol = list(pch = pch)), ...) {
    if (length(pvalues) == 0) stop("pvalue vector is empty, can't draw plot")
    if (!(class(pvalues) == "numeric" ||
          (class(pvalues) == "list" && all(sapply(pvalues, class) == "numeric"))))
        stop("pvalue vector is not numeric, can't draw plot")
    if (any(is.na(unlist(pvalues)))) stop("pvalue vector contains NA values, can't draw plot")
    if (already.transformed == FALSE) {
        if (any(unlist(pvalues) == 0)) stop("pvalue vector contains zeros, can't draw plot")
    } else {
        if (any(unlist(pvalues) < 0)) stop("-log10 pvalue vector contains negative values, can't draw plot")
    }
    grp <- NULL
    n <- 1
    exp.x <- c()
    if (is.list(pvalues)) {
        nn <- sapply(pvalues, length)
        rs <- cumsum(nn)
        re <- rs - nn + 1
        n <- min(nn)
        if (!is.null(names(pvalues))) {
            grp <- factor(rep(names(pvalues), nn), levels = names(pvalues))
            names(pvalues) <- NULL
        } else {
            grp <- factor(rep(1:length(pvalues), nn))
        }
        pvo <- pvalues
        pvalues <- numeric(sum(nn))
        exp.x <- numeric(sum(nn))
        for (i in 1:length(pvo)) {
            if (!already.transformed) {
                pvalues[rs[i]:re[i]] <- -log10(pvo[[i]])
                exp.x[rs[i]:re[i]] <- -log10((rank(pvo[[i]], ties.method = "first") - .5) / nn[i])
            } else {
                pvalues[rs[i]:re[i]] <- pvo[[i]]
                exp.x[rs[i]:re[i]] <- -log10((nn[i] + 1 - rank(pvo[[i]], ties.method = "first") - .5) / (nn[i] + 1))
            }
        }
    } else {
        n <- length(pvalues) + 1
        if (!already.transformed) {
            exp.x <- -log10((rank(pvalues, ties.method = "first") - .5) / n)
            pvalues <- -log10(pvalues)
        } else {
            exp.x <- -log10((n - rank(pvalues, ties.method = "first") - .5) / n)
        }
    }
    panel.qqconf <- function(n, conf.points = 1000, conf.col = "gray", conf.alpha = .05, ...) {
        require(grid)
        conf.points <- min(conf.points, n - 1)
        mpts <- matrix(nrow = conf.points * 2, ncol = 2)
        for (i in seq(from = 1, to = conf.points)) {
            mpts[i, 1] <- -log10((i - .5) / n)
            mpts[i, 2] <- -log10(qbeta(1 - conf.alpha / 2, i, n - i))
            mpts[conf.points * 2 + 1 - i, 1] <- -log10((i - .5) / n)
            mpts[conf.points * 2 + 1 - i, 2] <- -log10(qbeta(conf.alpha / 2, i, n - i))
        }
        grid.polygon(x = mpts[, 1], y = mpts[, 2], gp = gpar(fill = conf.col, lty = 0), default.units = "native")
    }
    if (should.thin == T) {
        if (!is.null(grp)) {
            thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                                      exp.x = round(exp.x, thin.exp.places), grp = grp))
            grp <- thin$grp
        } else {
            thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                                      exp.x = round(exp.x, thin.exp.places)))
        }
        pvalues <- thin$pvalues
        exp.x <- thin$exp.x
    }
    gc()
    prepanel.qqunif <- function(x, y, ...) {
        A <- list()
        A$xlim <- range(x, y) * 1.02
        A$xlim[1] <- 0
        A$ylim <- A$xlim
        return(A)
    }
    xyplot(pvalues ~ exp.x, groups = grp, xlab = xlab, ylab = ylab, aspect = aspect,
           prepanel = prepanel, 
           
           #scales = list(axs = "i"), 
           
           pch = pch,
           panel = function(x, y, ...) {
               if (draw.conf) {
                   panel.qqconf(n, conf.points = conf.points,
                                conf.col = conf.col, conf.alpha = conf.alpha)
               }
               panel.xyplot(x, y, ...)
               panel.abline(0, 1)
           }, par.settings = par.settings, ...)
}
