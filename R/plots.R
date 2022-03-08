computeScaleFactor<- function(object, sizes) {
    if (is.null(object)) {
        length <- sizes[[2]] - sizes[[1]]
    }
    else if (is(object, "tenxcheckerExp")) {
        length <- sum(sizes)
    }
    else if (is(object, "tenxchecker2RefExp")) {
        length <- min(object@size1, object@size2)
    }
    else if (is(object, "tenxcheckerRefExp")) {
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
            geom_vline(xintercept = object@parameters@maxLinkRange, linetype = "dashed") +
            xlim(1, max(30, 2 * object@parameters@maxLinkRange)) +
            scale_y_log10()
}

.plotDiagonalStrengthFit <- function(object) {
    if (! is(object, "tenxcheckerExp")) {
        stop("Object should be a 'tenxcheckerExp'.")
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
    if (! is(object, "tenxcheckerClass")) {
        stop("Object should be a 'tenxcheckerClass'.")
    }
    plots <- purrr::map(object@data, .plotDiagonalStrengthFit)
    return(do.call("plot_grid", c(plots, ncol = length(plots))))
}

plotBreakStats <- function(object, dataName, refName, bins = NULL) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Object should be a 'tenxcheckerClass'.")
    }
    p <- getDataset(object, dataName)@breaks@data %>%
        dplyr::filter(ref == refName) %>%
        tidyr::gather(key = "type", value = "value", fcMeanCount, nCells) %>%
        ggplot(aes(x = bin, y = value)) +
            geom_line() +
            facet_grid(rows = vars(type), scales = "free_y")
    if (!is.null(bins)) {
        for (bin in bins) {
          p <- p + geom_vline(xintercept = bin,
                              colour = "red",
                              linetype = "longdash")
        }
    }
    p
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

plotMD <- function(object, log) {
    samples    <- map_chr(object@data, "name")
    sampleSize <- sum(map_dbl(map(object@data, "parameters"), "sampleSize"))
    p <- map_dfr(object@data, "interactionMatrix", .id = "sample") %>%
        filter(ref1 == ref2) %>%
        sample_n(min(sampleSize, nrow(.))) %>%
        mutate(sample = factor(samples[as.numeric(sample)], levels = samples)) %>%
        mutate(distance = abs(bin1 - bin2)) %>%
        ggplot(aes(x = distance, y = count)) + 
            xlim(0, 100) + ylim(0, 50) +
            geom_bin2d(binwidth = c(1, 1)) +
            geom_smooth(method = "loess", formula = y ~ x) +
	    facet_grid(cols = vars(sample), scale = "free", space = "free")
    if (log) {
        p <- p + scale_y_log10()
    }
    return(p)
}

plotDiagonalStrength <- function(object) {
    tmp <- object@interactionMatrix %>%
        filter(ref1 == ref2) %>%
        makeSymmetric() %>%
        mutate(distance = abs(bin1 - bin2)) %>%
        dplyr::select(ref1, bin1, distance, count) %>%
        arrange(ref1, bin1, distance)
    tmp1 <- tmp %>%
        group_by(ref1, bin1) %>%
        mutate(diagonalSum = cumsum(count)) %>%
        ungroup() %>%
        dplyr::select(ref1, bin1, distance, diagonalSum)
    p <- tmp %>%
        group_by(ref1, bin1) %>%
        summarise(countSum = sum(count)) %>%
        ungroup() %>%
        right_join(tmp1, by = c("ref1", "bin1")) %>%
        mutate(diagPc = diagonalSum / countSum) %>%
        filter(distance <= 2 * object@parameters@maxLinkRange) %>%
        dplyr::select(ref1, bin1, distance, diagPc) %>%
        dplyr::rename(ref = ref1) %>%
        dplyr::rename(bin = bin1) %>%
        group_by(ref)
    n <- p %>% group_keys()
    p <- p %>%
        group_split(keep = FALSE) %>%
        map(~ ggplot(., aes(x = bin, y = distance, fill = diagPc)) + 
                geom_tile())
    names(p) <- n$ref
    return(p)
}

plotInsertions1 <- function(table, ref) {
    table %>% ggplot(aes(x = bin, y = distance, fill = difference)) + 
        geom_tile() +
        scale_fill_gradient2() +
        coord_fixed() +
        ggtitle(ref)
}

plotInsertions2 <- function(object, bin, distance) {
    df <- tibble(xmin = bin - distance,
                 xmax = bin + distance,
                 ymin = bin - distance,
                 ymax = bin + distance)
    plot.10XRef(object, logColor = TRUE) +
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
            as_tibble() %>%                                                                                                                                                                                        
            # do not use the last bin of the ref                                                                                                                                                                   
            dplyr::mutate(size = sizes[ref]) %>%                                                                                                                                                                   
            dplyr::filter(size - bin >= object@parameters@metaSize)                                                                                                                                                
    }                                                                                                                                                                                                              
    else {
        colSums <- object@interactionMatrix %>%
            computeSymmetricColSum(sizes) %>%
            as_tibble()
    }
    tmp <- colSums %>%
        dplyr::filter(sum > 0) %>%
        dplyr::pull(sum)
    quartiles        <- quantile(tmp, prob = c(.25, .75))
    iqr       <- quartiles[[2]] - quartiles[[1]]                                                                                                                                                                   
    firstOutlier <- quartiles[[1]] - 1.5 * iqr                                                                                                                                                                     
    lastOutlier  <- quartiles[[2]] + 1.5 * iqr                                                                                                                                                                     
    tmp <- tmp %>%                                                                                                                                                                                                 
        keep(~ .x >= firstOutlier) %>%                                                                                                                                                                             
        keep(~ .x <= lastOutlier)                                                                                                                                                                                  
    fitNB            <- fitdistr(tmp, "negative binomial")
    size             <- fitNB$estimate[[1]]
    mu               <- fitNB$estimate[[2]]
    ggplot(colSums, aes(x = sum)) + 
        geom_density(aes(y = stat(density))) +
        stat_function(fun = function(x){dnbinom(round(x), size = size, mu = mu)}, col = "red") + 
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
    return(do.call("plot_grid", c(plots, ncol = length(plots))))
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

plot.10XRef <- function(object, logColor = TRUE, binMin = NULL, binMax = NULL, bins = NULL, outliers = TRUE, meta = FALSE) {
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
    minCount <- data %>% pull(count) %>% min()
    if (meta) {
        data <- transformMeta(data, object@parameters@metaSize)
    }
    if ((minCount < 0) & (logColor)) {
        stop("Trying to plot a map with negative count in log scale.")
    }
    if (logColor) {
        data %<>% dplyr::filter(count > 0)
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
    if (! is_null(bins)) {
        p <- p +
            geom_hline(yintercept = bins, linetype = "dotted") +
            geom_vline(xintercept = bins, linetype = "dotted")
    }    
    return(p)
}

plot.10X2Ref <- function(object,
			 logColor = TRUE,
			 circles = FALSE,
			 center  = NULL,
			 radius  = NULL,
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
    if (! is.null(center)) {
	if (length(center) != 2) {
	    stop(paste0("'center' parameter should have size 2 (size ", length(center), " found)."))
	}
	x1 <- max(x1, rescaleValue(center[[1]] - radius, scaleFactor))
	x2 <- min(x2, rescaleValue(center[[1]] + radius, scaleFactor))
	y1 <- max(y1, rescaleValue(center[[2]] - radius, scaleFactor))
	y2 <- min(y2, rescaleValue(center[[2]] + radius, scaleFactor))
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
	circles <- tibble(
	    x      = c(1, 1, object@size1, object@size1),
	    y      = c(1, object@size2, 1, object@size2),
	    radius = rep.int(object@parameters@maxLinkRange, 4)) %>%
	    transpose()
	addCircle <- function (plot, parameters) {
	    nPoints <- 1000
	    lim1 <- object@size1 / 2
	    lim2 <- object@size2 / 2
	    circle <- bind_rows(tibble(x = parameters$x - parameters$radius + seq(0, parameters$radius),
				       y = parameters$y + seq(0, parameters$radius)),
				tibble(x = parameters$x - parameters$radius + seq(0, parameters$radius),
				       y = parameters$y - seq(0, parameters$radius)),
				tibble(x = parameters$x + seq(0, parameters$radius),
				       y = parameters$y + parameters$radius - seq(0, parameters$radius)),
				tibble(x = parameters$x + seq(0, parameters$radius),
				       y = parameters$y - parameters$radius + seq(0, parameters$radius))) %>%
		mutate(x = if_else((parameters$x < lim1) == (x < lim1), x, lim1)) %>%
		mutate(y = if_else((parameters$y < lim2) == (y < lim2), y, lim2)) %>%
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

plot.10XDataset <- function(object, sizes, logColor = TRUE, ref1 = NULL, ref2 = NULL, bin1 = NULL, bin2 = NULL, outliers = TRUE, highlightedBins = c(), meta = FALSE) {
    if (! is(object, "tenxcheckerExp")) {
        stop("Object should be a 'tenxcheckerExp'.")
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
        return(plot.10XRef(object, logColor = logColor, binMin = bin1, binMax = bin2, outliers = outliers, bins = highlightedBins, meta = meta))
    }
    if (! ref2 %in% levels(object@interactionMatrix$ref1)) {
        message(paste0("Reference #2 '", ref2, "', is not a known reference in dataset '", object@name, "'."))
    }
    if (match(ref1, object@chromosomes) < match(ref2, object@chromosomes)) {
       tmp <- ref1
       ref1 <- ref2
       ref2 <- tmp
       tmp <- bin1
       bin1 <- bin2
       bin2 <- tmp
    }
    object <- extract2Ref(object, ref1, ref2, sizes[[ref1]], sizes[[ref2]])
    return(plot.10X2Ref(object, logColor = logColor, outliers = outliers))
}

plot.10X <- function(object, sizes = NULL, logColor = TRUE, datasetName = NULL, ref1 = NULL, ref2 = NULL, bin1 = NULL, bin2 = NULL, outliers = TRUE, highlightedBins = c(), meta = FALSE) {
    if (is(object, "tenxcheckerClass")) {
        sizes  <- object@sizes
        if (is.null(datasetName)) {
             plots <- purrr::map(object@data, plot.10XDataset, sizes, logColor, ref1, ref2, bin1, bin2, outliers, highlightedBins, meta)
             return(do.call("plot_grid", c(plots, ncol = length(plots))))
        }
        else {
            datasetNames <- map(object@data, "name")
            if (datasetName %in% datasetNames) {
                dataset <- object@data[datasetName == datasetNames][[1]]
                return(plot.10XDataset(dataset, object@sizes, logColor, ref1, ref2, bin1, bin2, outliers, highlightedBins, meta))
            }
            else {
                stop(paste0("Dataset name '", datasetName, "' is not known."))
            }
        }
    }
    else if (! is(object, "tenxcheckerExp")) {
        stop("Object should be a 'tenxcheckerClass' or 'tenxcheckerExp'.")
    }
    return(plot.10XDataset(object, sizes, logColor, ref1, ref2, bin1, bin2, outliers, highlightedBins, meta))
}

# This will plot the region near a possible break
plot.10XBreak <- function(object, ref, bin) {
    nBinZoom <- max(map_dbl(map(object@data, "parameters"), "nBinZoom"))
    plot.10X(object = object, ref1 = ref, bin1 = bin - nBinZoom, bin2 = bin + nBinZoom, outliers = FALSE, highlightedBins = bin, meta = TRUE)
}
