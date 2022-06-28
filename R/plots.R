computeScaleFactor<- function(object, sizes) {
    if (is.null(object)) {
        length <- sizes[[2]] - sizes[[1]]
    }
    else if (is(object, "msscafExp")) {
        length <- sum(sizes)
    }
    else if (is(object, "msscaf2RefExp")) {
        length <- min(object@size1, object@size2)
    }
    else if (is(object, "msscafRefExp")) {
        length <- object@size
    }
    else {
        stop("Do not know what to do with object.")
    }
    scaleFactor <- ceiling(log10(length))
    if (scaleFactor <= 3) {
        return(1)
    }
    return(10^(scaleFactor - 3))
}

rescale <- function(data, scale) {
    if (scale == 1) {
        return(data)
    }
    data %<>%
        dplyr::mutate(bin1 = as.integer(round(bin1 / scale) * scale)) %>%
        dplyr::mutate(bin2 = as.integer(round(bin2 / scale) * scale)) %>%
        dplyr::mutate(count = as.numeric(count))
    if ("ref1" %in% colnames(data)) {
        data %<>%
            dplyr::group_by(ref1, ref2, bin1, bin2) %>%
            dplyr::summarise(count = mean(count)) %>%
            dplyr::ungroup()
        return(data)
    }
    data %>%
        dplyr::group_by(bin1, bin2) %>%
        dplyr::summarise(count = mean(count)) %>%
        dplyr::ungroup()
}

rescaleValue <- function(value, scale) {
    as.integer(round(value / scale) * scale)
}

plotMoleculeSizeDistribution <- function(object) {
    object@interactionMatrix %>%
        filter(ref1 == ref2) %>%
        mutate(distance = abs(bin1 - bin2)) %>%
        dplyr::select(distance, count) %>%
        sample_n(min(object@parameters@sampleSize, nrow(.))) %>%
        mutate(loess = predict(loess(count ~ distance, data = ., span = 0.01))) %>%
            ggplot(aes(distance, count)) +
            geom_point(color = "grey50") +
            geom_line(aes(x = distance, y = loess)) + 
            geom_hline(yintercept = object@parameters@minCount, linetype = "dashed") +
            geom_vline(xintercept = object@parameters@minLinkRange, linetype = "dashed") +
            xlim(1, max(30, 2 * object@parameters@minLinkRange)) +
            scale_y_log10()
}

.plotDiagonalStrengthFit <- function(object) {
    if (! is(object, "msscafExp")) {
        stop("Object should be a 'msscafExp'.")
    }
    tmp <- object@breaks@data %>%
        dplyr::filter(nCells >= object@parameters@breakNCells) %>%
        dplyr::filter(fcMeanCount >= 0) %>%
        dplyr::pull(fcMeanCount)
    standardDev <- sd(c(tmp, -tmp))
    object@breaks@data %>%
        dplyr::filter(nCells >= object@parameters@breakNCells) %>%
        ggplot(aes(fcMeanCount)) +
        geom_histogram(aes(y = stat(density))) +
        stat_function(fun = dnorm, args = list(mean = 0.0, sd = standardDev)) +
        ggtitle(object@name)
}

# This plot shows how triangles (near the diagonal) are fitted.
plotDiagonalStrengthFit <- function(object) {
    if (! is(object, "msscafClass")) {
        stop("Object should be a 'msscafClass'.")
    }
    plots <- purrr::map(object@data, .plotDiagonalStrengthFit)
    return(do.call(cowplot::plot_grid, c(plots, ncol = length(plots))))
}

.plotBreakStats <- function(object, refName, bins = NULL, minBin = NULL, maxBin = NULL) {
    if (! is(object, "msscafExp")) {
        stop("Object should be a 'msscafExp'.")
    }
    p <- object@breaks@data %>%
        dplyr::filter(ref == refName) %>%
        tidyr::gather(key = "type", value = "value", fcMeanCount, nCells) %>%
        ggplot(aes(x = bin, y = value)) +
            geom_line() +
            facet_grid(rows = vars(type), scales = "free_y") +
            ggtitle(object@name) +
            xlab(refName)
    if (!is.null(bins)) {
        for (bin in bins) {
          p <- p + geom_vline(xintercept = bin,
                              colour = "red",
                              linetype = "longdash")
        }
    }
    if (!is.null(minBin)) {
        if (is.null(maxBin)) {
            stop("'minBin' and 'maxBin' should be either both set, or both unset.")
        }
        p <- p + xlim(minBin, maxBin)
    }
    p
}

# This plot give the mean triangle values, and the #cells used for each triangle
plotBreakStats <- function(object, refName, bins = NULL, minBin = NULL, maxBin = NULL) {
    if (! is(object, "msscafClass")) {
        stop("Object should be a 'msscafClass'.")
    }
    plots <- purrr::map(object@data, .plotBreakStats, refName = refName, bins = bins, minBin = minBin, maxBin = maxBin)
    return(do.call(cowplot::plot_grid, c(plots, ncol = length(plots))))
}


plotBackgroundCountDistribution <- function(object) {
    object@interactionMatrix %>%
        filter(ref1 != ref2) %>%
        dplyr::select(count) %>%
        sample_n(min(object@parameters@sampleSize, nrow(.))) %>%
        ggplot(aes(count)) +
            geom_freqpoly(binwidth=1) +
            geom_vline(xintercept = object@parameters@minCount, linetype = "dashed") +
            xlim(1, max(30, 2 * object@parameters@minCount)) +
            scale_y_log10()
}

plotRowCountDistribution <- function(object) {
    object@interactionMatrix %>%
        makeSymmetric() %>%
        group_by(ref1, bin1) %>%
        summarise(countSum = sum(count)) %>%
        ungroup() %>%
        dplyr::select(countSum) %>%
        ggplot(aes(countSum)) +
#            geom_freqpoly(binwidth=1) +
            geom_freqpoly() +
            geom_vline(xintercept = object@parameters@minRowCount, linetype = "dashed") +
            xlim(1, max(30, 2 * object@parameters@minRowCount))
}

plotScaffoldSizeDistribution <- function(object, log) {
    p <- object@interactionMatrix %>%
        filter(ref1 == ref2) %>%
        makeSymmetric() %>%
        group_by(ref1) %>%
        summarise(size = max(bin1)) %>%
        ungroup() %>%
        ggplot(aes(x = size)) + 
            geom_histogram()
    if (log) {
        p <- p + scale_y_log10()
    }
    return(p)
}

plotCountDistribution <- function(object, log) {
    p <- ggplot(object@interactionMatrix, aes(x = count)) + 
        geom_histogram(binwidth = 5)
    if (log) {
        p <- p + scale_y_log10()
    }
    return(p)
}

.plotDistanceCount <- function(object, log) {
    d <- object@interactionMatrix %>%
        dplyr::filter(ref1 == ref2) %>%
        dplyr::sample_n(min(object@parameters@sampleSize, nrow(.))) %>%
        dplyr::mutate(distance = abs(bin1 - bin2)) %>%
        dplyr::mutate(distance = as.numeric(round(distance / object@parameters@metaSize))) %>%
        dplyr::group_by(ref1, bin1, distance) %>%
        dplyr::summarize(count = sum(count), .groups = "drop")
    maxC <- unname(quantile(d$count, 0.99))
    p <- d %>%
        dplyr::filter(count <= maxC) %>%
        ggplot(aes(x = distance, y = count)) + 
            xlim(1, 4 * object@parameters@minLinkRange) +
            geom_bin2d(binwidth = c(1, 1)) +
            geom_smooth(method = "loess", formula = y ~ x) +
            ggtitle(object@name)
    if (log) {
        p <- p + scale_y_log10()
    }
    return(p)
        
#   samples    <- purrr::map_chr(object@data, "name")
#   sampleSize <- sum(purrr::map_dbl(purrr::map(object@data, "parameters"), "sampleSize"))
#   metaSizes  <- purrr::map_int(purrr::map(object@data, "parameters"), "metaSize")
#   p <- purrr::map_dfr(object@data, "interactionMatrix", .id = "sample") %>%
#       dplyr::filter(ref1 == ref2) %>%
#       dplyr::mutate(sample = factor(samples[as.numeric(sample)], levels = samples)) %>%
#       dplyr::mutate(metaSize = metaSizes[as.numeric(sample)]) %>%
#       dplyr::mutate(distance = abs(bin1 - bin2)) %>%
#       dplyr::mutate(distance = distance / metaSize) %>%
#       dplyr::group_by(sample, ref1, bin1, distance) %>%
#       dplyr::summarize(count = sum(count), .groups = "drop") %>%
#       dplyr::sample_n(min(sampleSize, nrow(.))) %>%
#       ggplot(aes(x = distance, y = count)) + 
#           xlim(0, 100) + ylim(0, 50) +
#           geom_bin2d(binwidth = c(1, 1)) +
#           geom_smooth(method = "loess", formula = y ~ x) +
#           facet_grid(cols = vars(sample), scale = "free", space = "free")
#   if (log) {
#       p <- p + scale_y_log10()
#   }
    return(p)
}

plotDistanceCount <- function(object, log) {
    plots <- purrr::map(object@data, .plotDistanceCount, log)
    return(do.call(cowplot::plot_grid, c(plots, ncol = length(plots))))
}

plotDiagonalStrength <- function(object) {
    tmp <- object@interactionMatrix %>%
        dplyr::filter(ref1 == ref2) %>%
        makeSymmetric() %>%
        dplyr::mutate(distance = abs(bin1 - bin2)) %>%
        dplyr::select(ref1, bin1, distance, count) %>%
        dplyr::arrange(ref1, bin1, distance)
    tmp1 <- tmp %>%
        dplyr::group_by(ref1, bin1) %>%
        dplyr::mutate(diagonalSum = cumsum(count)) %>%
        dplyr::ungroup() %>%
        dplyr::select(ref1, bin1, distance, diagonalSum)
    p <- tmp %>%
        dplyr::group_by(ref1, bin1) %>%
        dplyr::summarise(countSum = sum(count)) %>%
        dplyr::ungroup() %>%
        dplyr::right_join(tmp1, by = c("ref1", "bin1")) %>%
        dplyr::mutate(diagPc = diagonalSum / countSum) %>%
        dplyr::filter(distance <= 2 * object@parameters@minLinkRange) %>%
        dplyr::select(ref1, bin1, distance, diagPc) %>%
        dplyr::rename(ref = ref1) %>%
        dplyr::rename(bin = bin1) %>%
        dplyr::group_by(ref)
    n <- p %>% dplyr::group_keys()
    p <- p %>%
        dplyr::group_split(keep = FALSE) %>%
        purrr::map(~ ggplot(., aes(x = bin, y = distance, fill = diagPc)) + 
                geom_tile())
    names(p) <- n$ref
    return(p)
}

.plotBreakPvalueFit <- function(object) {
    tmp <- object@breaks@data %>%
        dplyr::filter(nCells >= object@parameters@breakNCells) %>%
        dplyr::filter(fcMeanCount >= 0) %>%
        dplyr::pull(fcMeanCount)
    standardDev <- sd(c(tmp, -tmp))
    object@breaks@data %>%
        dplyr::filter(nCells >= object@parameters@breakNCells) %>%
        ggplot(aes(fcMeanCount))+
            geom_histogram(aes(y = ..density..)) +
            stat_function(fun = dnorm, args = list(mean = 0, sd = standardDev)) +
            ggtitle(object@name)
}

plotBreakPvalueFit <- function(object, pvalue) {
    plots <- purrr::map(object@data, .plotBreakPvalueFit)
    return(do.call(cowplot::plot_grid, c(plots, ncol = length(plots))))
}

.plotCornerFit <- function(object) {
    nDistances <- nrow(object@parameters@cornerScores)
    colors <- rainbow(nDistances)
    distances <- seq.int(nDistances)
    p <- data.frame(x = seq.int(nDistances),
    		c = as.factor(as.character(distances))) %>%
    	 ggplot(aes(x))
    for (i in distances) {
        p <- p + stat_function(aes_(color = factor(as.character(i))), fun = dgamma, args = list(shape = object@parameters@cornerScores[[i, "shape"]], rate = object@parameters@cornerScores[[i, "rate"]]),)
    }
    p <- p + scale_colour_manual("Distance", values = colors, breaks = as.character(seq.int(nDistances)))
    return(p)
}

# Plot the corner count/distance distributions
plotCornerFit <- function(object) {
    if (! is(object, "msscafClass")) {
        stop("Object should be a 'msscafClass'.")
    }
    plots <- purrr::map(object@data, .plotCornerFit)
    return(do.call(cowplot::plot_grid, c(plots, ncol = length(plots))))
}

plotInsertions1 <- function(table, ref) {
    table %>% ggplot(aes(x = bin, y = distance, fill = difference)) + 
        geom_tile() +
        scale_fill_gradient2() +
        coord_fixed() +
        ggtitle(ref)
}

plotInsertions2 <- function(object, bin, distance) {
    df <- tibble::tibble(xmin = bin - distance,
                 xmax = bin + distance,
                 ymin = bin - distance,
                 ymax = bin + distance)
    plot.msscafRef(object, logColor = TRUE) +
        #annotate("rect",
        geom_rect(data = df,
                  aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                  color = "black",
                  inherit.aes = FALSE,
                  alpha = 0.01)
}


.plotRowCountFit <- function(object, sizes) {
    if (object@parameters@metaSize > 1) {
        colSums <- object@interactionMatrix %>%
            computeSymmetricColSumMeta(sizes, object@parameters@metaSize) %>%
            tibble::as_tibble() %>%
            # do not use the last bin of the ref
            dplyr::mutate(size = sizes[ref]) %>%
            dplyr::filter(size - bin >= object@parameters@metaSize)
    }
    else {
        colSums <- object@interactionMatrix %>%
            computeSymmetricColSum(sizes) %>%
            tibble::as_tibble()
    }
    tmp <- colSums %>%
        dplyr::filter(sum > 0) %>%
        dplyr::pull(sum)
    quartiles    <- quantile(tmp, prob = c(.25, .75))
    iqr          <- quartiles[[2]] - quartiles[[1]]
    firstOutlier <- quartiles[[1]] - 1.5 * iqr
    lastOutlier  <- quartiles[[2]] + 1.5 * iqr
    ggplot(colSums, aes(x = sum)) + 
        geom_density(aes(y = stat(density))) +
        stat_function(fun = function(x){dnbinom(round(x), size = object@parameters@rowCountSize, mu = object@parameters@rowCountMu)}, col = "red") + 
        xlim(0, max(lastOutlier, object@parameters@maxRowCount)) +
        xlab(object@name) +
        geom_vline(xintercept = object@parameters@minRowCount, linetype = "dashed") +
        geom_vline(xintercept = object@parameters@maxRowCount, linetype = "dotted") +
        geom_vline(xintercept = firstOutlier, linetype = "dashed", color = "red") +
        geom_vline(xintercept = lastOutlier,  linetype = "dotted", color = "red")
}

# This plot shows how column sums are fitted
plotRowCountFit <- function(object) {
    plots <- purrr::map(object@data, .plotRowCountFit, object@sizes)
    return(do.call(cowplot::plot_grid, c(plots, ncol = length(plots))))
}

.computeRowDensity <- function(chr, data = data) {
    data %>%
        makeSymmetric() %>%
        filter(ref1 == chr) %>%
        group_by(ref1, bin1) %>%
        summarise(countSum = sum(count)) %>%
        ungroup()
}

.plotOneRowCountDensity <- function(data) {
    meanValue = mean(data$countSum)
    varValue = var(data$countSum)
    ggplot(data, aes(x = countSum)) + 
        geom_density(aes(y = stat(density))) +
        geom_line(aes(y = dnbinom(countSum, mu = meanValue, size = meanValue ^ 2 / (varValue - meanValue)), color = "red")) +
        scale_color_discrete(labels = c("expected"))
}

.plotRowCountDensity <- function(chr, data = data) {
    .plotOneRowCountDensity(.computeRowDensity(chr, data = data))
}

plotRowCountDensity <- function(object) {
    lapply(object@chromosomes, .plotRowCountDensity,
           data = object@interactionMatrix)
}

.plot.msscafJoin <- function(object, r1, r2, after1, after2, size1, size2, logColor, zoom, outliers = outliers, metaSize = metaSize) {
    data <- object@interactionMatrix
    outlierBins <- object@outlierBins
    newRefs <- c(r1, r2)
    sizes   <- c(size1, size2)
    if (! after1) {
        newRefs <- c(r2, r1)
        sizes   <- c(size2, size1)
    }
    data <- data %>%
        dplyr::mutate(ref1 = forcats::fct_relevel(ref1, newRefs)) %>%
        dplyr::mutate(ref2 = forcats::fct_relevel(ref2, newRefs))
    outlierBins <- outlierBins %>%
        dplyr::mutate(ref = forcats::fct_relevel(ref, newRefs))
    if (metaSize) {
        data <- transformMeta(data, object@parameters@metaSize)
    }
    data <- data %>%
        makeSymmetric()
    # Removing outliers requires the matrix to be symmetric
    # It also inserts NAs on the full chromosomes
    if (! outliers) {
        data %<>% removeOutliersCpp(outlierBins, sizes) %>% tibble::as_tibble()
    }
    if (after1) {
        data <- data %>%
            dplyr::filter((ref1 != r1) | (bin1 >= size1 - zoom)) %>%
            dplyr::filter((ref2 != r1) | (bin2 >= size1 - zoom))
    }
    else {
        data <- data %>%
            dplyr::filter((ref1 != r1) | (bin1 <= zoom)) %>%
            dplyr::filter((ref2 != r1) | (bin2 <= zoom))
    }
    if (after2) {
        data <- data %>%
            dplyr::filter((ref1 != r2) | (bin1 >= size2 - zoom)) %>%
            dplyr::filter((ref2 != r2) | (bin2 >= size2 - zoom))
    }
    else {
        data <- data %>%
            dplyr::filter((ref1 != r2) | (bin1 <= zoom)) %>%
            dplyr::filter((ref2 != r2) | (bin2 <= zoom))
    }
    if (logColor) {
        data %<>% dplyr::filter(is.na(count) | (count > 0))
    }
    if (after1 == after2) {
        data <- data %>%
            dplyr::mutate(bin1 = dplyr::if_else(ref1 == r2, -bin1, bin1)) %>%
            dplyr::mutate(bin2 = dplyr::if_else(ref2 == r2, -bin2, bin2))
    }
    if (nrow(data) == 0) return(NULL)
    p <- data %>% 
            ggplot(aes(x = bin1, y = bin2)) + 
            geom_raster(aes(fill = count)) + 
            facet_grid(cols = vars(ref1), rows = vars(ref2), scale = "free", space = "free") +
            scale_x_continuous(expand = c(0, 0)) +
            scale_y_reverse(expand = c(0, 0)) +
            theme_bw() +
            theme(panel.spacing = unit(0, "lines")) +
            ggtitle(object@name)
    if (logColor) {
        p <- p + scale_fill_gradient(low = "grey90", high = "red", trans = "log")
    } else {
        p <- p + scale_fill_gradient2()
    }
    return(p)
}

# This will plot the region near a possible join
plot.msscafJoin <- function(object, ref1, ref2, after1, after2, logColor = TRUE, outliers = FALSE, metaSize = FALSE, nBinZoom = NULL) {
    if (! is(object, "msscafClass")) {
        stop("Object should be a 'msscafClass'.")
    }
    if (ref1 == ref2) {
        stop("Reference names should be different.")
    }
    if (is.null(nBinZoom)) {
        nBinZoom <- max(purrr::map_dbl(purrr::map(object@data, "parameters"), "nBinZoom"))
    }
    object <- keepScaffolds(object, c(ref1, ref2))
    plots <- purrr::map(object@data, .plot.msscafJoin, ref1, ref2, after1, after2, object@sizes[[ref1]], object@sizes[[ref2]], logColor, nBinZoom, outliers, metaSize) %>%
        purrr::discard(purrr::is_null)
    return(do.call(cowplot::plot_grid, c(plots, ncol = length(plots))))
}

# This will plot the region near a possible break
plot.msscafBreak <- function(object, ref, bin, outliers = FALSE, metaSize = FALSE, nBinZoom = NULL) {
    if (is.null(nBinZoom)) {
        nBinZoom <- max(purrr::map_dbl(purrr::map(object@data, "parameters"), "nBinZoom"))
    }
    plot.msscaf(object = object, ref1 = ref, bin1 = bin - nBinZoom, bin2 = bin + nBinZoom, outliers = outliers, highlightedBins = bin, meta = metaSize)
}

# Agregate bins to meta bins
transformMeta <- function(data, metaSize) {
    if (metaSize == 1) return(data)
    data <- data %>%
        dplyr::mutate(metaBin1 = as.integer(trunc(bin1 / metaSize))) %>%
        dplyr::mutate(metaBin2 = as.integer(trunc(bin2 / metaSize)))
    if ("ref1" %in% colnames(data)) {
        if ("ref2" %in% colnames(data)) {
            data <- data %>%
                dplyr::group_by(ref1, ref2, metaBin1, metaBin2)
        }
        else {
            data <- data %>%
                dplyr::group_by(ref1, metaBin1, metaBin2)
        }
    }
    else {
        data <- data %>%
            dplyr::group_by(metaBin1, metaBin2)
    }
    data %>%
        dplyr::summarise(count = sum(count), .groups = "drop") %>%
        dplyr::mutate(metaBin1 = metaBin1 * metaSize) %>%
        dplyr::mutate(metaBin2 = metaBin2 * metaSize) %>%
        dplyr::rename(bin1 = metaBin1, bin2 = metaBin2)
}

plot.msscafRef <- function(object, logColor = TRUE, binMin = NULL, binMax = NULL, bins = NULL, outliers = TRUE, meta = FALSE, linkRanges = FALSE) {
    if ((! is.null(binMin)) & (! is.null(binMax))) {
        if (binMax < binMin) {
            stop(paste0("Second bin (", binMax, ") should be less than first bin (", binMin, ")."))
        }
    }
    if (length(bins) == 0) {
        bins <- NULL
    }
    if (is.null(binMin)) {
        scaleFactor <- computeScaleFactor(object)
        binMin <- 0
        binMax <- object@size
    }
    else {
        if (binMin < 0) {
            binMin <- 0
        }
        if (binMax > object@size) {
            binMax <- object@size
        }
        scaleFactor <- computeScaleFactor(NULL, c(binMin, binMax))
    }
    data     <- object@interactionMatrix %>%
                    dplyr::filter(bin1 >= binMin, bin1 <= binMax) %>%
                    dplyr::filter(bin2 >= binMin, bin2 <= binMax) %>%
                    rescale(scaleFactor)
    minCount <- data %>% dplyr::pull(count) %>% min()
    if (meta) {
        data <- transformMeta(data, object@parameters@metaSize)
    }
    if ((minCount < 0) & (logColor)) {
        stop("Trying to plot a map with negative count in log scale.")
    }
    if (logColor) {
        data %<>% dplyr::filter(is.na(count) | (count > 0))
    }
    data %<>% makeSymmetric()
    if (! outliers) {
        data %<>% removeOutliersRef(object@outlierBins, object@size, binMin, binMax)
    }
    p <- data %>% 
        ggplot(aes(x = bin1, y = bin2)) + 
            geom_raster(aes(fill = count)) +
	    scale_x_continuous(lim = c(binMin, binMax), expand = c(0, 0)) +
	    scale_y_reverse(lim = c(binMax, binMin), expand = c(0, 0)) +
            theme_bw() +
            theme(panel.spacing = unit(0, "lines")) +
            ggtitle(object@name) +
            xlab(object@chromosome) +
            ylab(object@chromosome)
    if (logColor) {
        p <- p + scale_fill_gradient(low = "grey90", high = "red", trans = "log")
    } else {
        p <- p + scale_fill_gradient2()
    }
    if (length(bins) > 100) {
        message(paste0("Cannot display ", length(bins), " lines in the plots: too much RAM needed."))
        bins <- NULL
    }
    if (! rlang::is_null(bins)) {
        p <- p +
            geom_hline(yintercept = bins, linetype = "dotted") +
            geom_vline(xintercept = bins, linetype = "dotted")
    }
    if (linkRanges) {
        p <- p +
            geom_abline(slope = -1, intercept = object@parameters@minLinkRange, linetype = "dashed") +
            geom_abline(slope = -1, intercept = object@parameters@maxLinkRange, linetype = "dotted") +
            geom_abline(slope = -1, intercept = - object@parameters@minLinkRange, linetype = "dashed") +
            geom_abline(slope = -1, intercept = - object@parameters@maxLinkRange, linetype = "dotted")
    }
    return(p)
}

plot.msscaf2Ref <- function(object,
			 logColor = TRUE,
			 circles  = FALSE,
			 b1       = NULL,
			 b2       = NULL,
			 radius   = NULL,
                         outliers = TRUE) {
    if (is.null(radius)) {
        scaleFactor <- computeScaleFactor(object)
    }
    else {
        scaleFactor <- computeScaleFactor(NULL, c(2 * radius, 2 * radius))
    }
    data        <- object@interactionMatrix %>% rescale(scaleFactor)
    x1          <- 0
    x2          <- rescaleValue(object@size1, scaleFactor)
    y1          <- 0
    y2          <- rescaleValue(object@size2, scaleFactor)
    if (! is.null(b1)) {
	x1 <- max(x1, rescaleValue(b1 - radius, scaleFactor))
	x2 <- min(x2, rescaleValue(b1 + radius, scaleFactor))
    }
    if (! is.null(b2)) {
	y1 <- max(y1, rescaleValue(b2 - radius, scaleFactor))
	y2 <- min(y2, rescaleValue(b2 + radius, scaleFactor))
    }
    p <- data %>%
	ggplot(aes(x = bin1, y = bin2)) + 
	    geom_raster(aes(fill = count)) + 
	    scale_x_continuous(lim = c(x1, x2), expand = c(0, 0)) +
	    scale_y_reverse(lim = c(y2, y1), expand = c(0, 0)) +
	    xlab(object@chromosome1) +
	    ylab(object@chromosome2) +
	    theme_bw() +
	    theme(panel.spacing = unit(0, "lines")) +
	    coord_fixed() +
            ggtitle(object@name) +
            xlab(object@chromosome1) +
            ylab(object@chromosome2)
    if (circles) {
	circles <- tibble::tibble(
	    x      = c(1, 1, object@size1, object@size1),
	    y      = c(1, object@size2, 1, object@size2),
	    radius = rep.int(object@parameters@minLinkRange, 4)) %>%
	    transpose()
	addCircle <- function (plot, parameters) {
	    nPoints <- 1000
	    lim1 <- object@size1 / 2
	    lim2 <- object@size2 / 2
	    circle <- bind_rows(tibble::tibble(x = parameters$x - parameters$radius + seq(0, parameters$radius),
				       y = parameters$y + seq(0, parameters$radius)),
				tibble::tibble(x = parameters$x - parameters$radius + seq(0, parameters$radius),
				       y = parameters$y - seq(0, parameters$radius)),
				tibble::tibble(x = parameters$x + seq(0, parameters$radius),
				       y = parameters$y + parameters$radius - seq(0, parameters$radius)),
				tibble::tibble(x = parameters$x + seq(0, parameters$radius),
				       y = parameters$y - parameters$radius + seq(0, parameters$radius))) %>%
		mutate(x = dplyr::if_else((parameters$x < lim1) == (x < lim1), x, lim1)) %>%
		mutate(y = dplyr::if_else((parameters$y < lim2) == (y < lim2), y, lim2)) %>%
		filter(x >= x1) %>%
		filter(y >= y1) %>%
		filter(x <= x2) %>%
		filter(y <= y2)
	    plot + annotate(geom = "point", x = circle$x, y = circle$y, size = 0.1, colour = "grey", alpha = 0.5)
	}
	p <- purrr::reduce(circles, addCircle, .init = p)
    }
#   if (!is.null(x2)) {
#       p <- p + geom_vline(xintercept = x2, linetype = "dashed", color = "red")
#   }
#   if (!is.null(y1)) {
#       p <- p + geom_hline(yintercept = y1, linetype = "dotted", color = "red")
#   }
#   if (!is.null(y2)) {
#       p <- p + geom_hline(yintercept = y2, linetype = "dashed", color = "red")
#   }
    if (logColor) {
	p <- p + scale_fill_gradient(low = "grey90", high = "red", trans = "log")
    }
    else {
	p <- p + scale_fill_gradient(low = "blue", high = "red")
    }
    return(p)
}

plot.msscafDataset <- function(object, sizes, logColor = TRUE, ref1 = NULL, ref2 = NULL, bin1 = NULL, bin2 = NULL, outliers = TRUE, highlightedBins = c(), meta = FALSE, radius = NULL, linkRanges = FALSE) {
    if (! is(object, "msscafExp")) {
        stop("Object should be a 'msscafExp'.")
    }
    if (is.null(bin1) != is.null(bin2)) {
        stop("None, or both bins should be set.")
    }
    if (is.null(ref1)) {
        scaleFactor <- computeScaleFactor(object, sizes)
        # message(paste0("Scale factor: ", scaleFactor))
        data        <- object@interactionMatrix %>% rescale(scaleFactor)
        p <- data %>%
            makeSymmetric() %>%
            ggplot(aes(x = bin1, y = bin2)) + 
            geom_raster(aes(fill = count)) + 
            facet_grid(cols = vars(ref1), rows = vars(ref2), scale = "free", space = "free") +
            scale_x_continuous(expand = c(0, 0)) +
            scale_y_reverse(expand = c(0, 0)) +
            theme_bw() +
            theme(panel.spacing = unit(0, "lines")) +
            ggtitle(object@name) +
            xlab("") +
            ylab("")
        if (logColor) {
            p <- p + scale_fill_gradient(low = "grey90", high = "red", trans = "log")
        }
        else {
            p <- p + scale_fill_gradient(low = "blue", high = "red")
        }
        return(p)
    }
    if (! ref1 %in% levels(object@interactionMatrix$ref1)) {
        message(paste0("Reference #1 '", ref1, "', is not a known reference in dataset '", object@name, "'."))
    }
    if (is.null(ref2)) {
        ref2 <- ref1
    }
    if (ref1 == ref2) {
        object <- extractRef(object, ref1, sizes[[ref1]])
        return(plot.msscafRef(object, logColor = logColor, binMin = bin1, binMax = bin2, outliers = outliers, bins = highlightedBins, meta = meta, linkRanges = linkRanges))
    }
    if (! ref2 %in% levels(object@interactionMatrix$ref1)) {
        message(paste0("Reference #2 '", ref2, "', is not a known reference in dataset '", object@name, "'."))
    }
    if (match(ref1, levels(object@interactionMatrix$ref1)) < match(ref2, levels(object@interactionMatrix$ref1))) {
       tmp <- ref1
       ref1 <- ref2
       ref2 <- tmp
       tmp <- bin1
       bin1 <- bin2
       bin2 <- tmp
    }
    object <- extract2Ref(object, ref1, ref2, sizes[[ref1]], sizes[[ref2]])
    return(plot.msscaf2Ref(object, logColor = logColor, outliers = outliers, b1 = bin1, b2 = bin2, radius = radius))
}

plot.msscaf <- function(object, sizes = NULL, logColor = TRUE, datasetName = NULL, ref1 = NULL, ref2 = NULL, bin1 = NULL, bin2 = NULL, outliers = TRUE, highlightedBins = c(), meta = FALSE, radius = NULL, linkRanges = FALSE) {
    if (is(object, "msscafClass")) {
        sizes  <- object@sizes
        if (is.null(datasetName)) {
             plots <- purrr::map(object@data, plot.msscafDataset, sizes, logColor, ref1, ref2, bin1, bin2, outliers, highlightedBins, meta, radius, linkRanges)
             return(do.call(cowplot::plot_grid, c(plots, ncol = length(plots))))
        }
        else {
            datasetNames <- purrr::map(object@data, "name")
            if (datasetName %in% datasetNames) {
                dataset <- object@data[datasetName == datasetNames][[1]]
                return(plot.msscafDataset(dataset, object@sizes, logColor, ref1, ref2, bin1, bin2, outliers, highlightedBins, meta, radius, linkRanges))
            }
            else {
                stop(paste0("Dataset name '", datasetName, "' is not known."))
            }
        }
    }
    else if (! is(object, "msscafExp")) {
        stop("Object should be a 'msscafClass' or 'msscafExp'.")
    }
    return(plot.msscafDataset(object, sizes, logColor, ref1, ref2, bin1, bin2, outliers, highlightedBins, meta, radius, linkRanges))
}
