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
    p <- object@interactionMatrix %>%
        filter(ref1 == ref2) %>%
        mutate(distance = abs(bin1 - bin2)) %>%
        ggplot(aes(x = distance, y = count)) + 
            geom_bin2d(binwidth = c(1, 10)) +
            geom_smooth()
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


plotRowCounts <- function(object) {
    p <- object@interactionMatrix %>%
        makeSymmetric() %>%
        group_by(ref1, bin1) %>%
        summarise(countSum = sum(count)) %>%
        ungroup() %>%
        group_by(ref1)
    n <- p %>% group_keys()
    p <- p %>%
        group_split(keep = FALSE) %>%
        map(~ ggplot(., aes(x = bin1, y = countSum)) + 
                geom_point() +
                geom_smooth())
                #facet_grid(rows = vars(ref1)) #+
    names(p) <- n$ref1
    p
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

plot.10XRef <- function(object, logColor = TRUE, bins = NULL, lim = NULL) {
    minLim <- 1
    maxLim <- object@size
    if (length(bins) == 0) {
        bins <- NULL
    }
    if (length(lim) == 2) {
        scaleFactor <- computeScaleFactor(NULL, lim)
    }
    else {
        scaleFactor <- computeScaleFactor(object)
    }
    if (length(lim) == 2) {
        minLim <- max(lim[[1]], minLim)
        maxLim <- min(lim[[2]], maxLim)
    }
    else if ((length(lim) != 1) & (! is_null(lim))) {
        stop("'lim' should be a vector of size 2, or NULL.")
    }
    data     <- object@interactionMatrix %>%
                    dplyr::filter(bin1 >= minLim, bin1 <= maxLim) %>%
                    dplyr::filter(bin2 >= minLim, bin2 <= maxLim) %>%
                    rescale(scaleFactor)
    minCount <- data %>% pull(count) %>% min()
    if ((minCount < 0) & (logColor)) {
        stop("Trying to plot a map with negative count in log scale.")
    }
    if (logColor) {
        data %<>% dplyr::filter(count > 0)
    }
    data %<>% makeSymmetric()
    # if (!is.na(bin)) {
    #     xmin <- max(0, bin - object@parameters@nBinZoom)
    #     xmax <- min(object@size, bin + object@parameters@nBinZoom)
    #     tmp %<>% filter(bin1 >= xmin,
    #                     bin1 <= xmax,
    #                     bin2 >= xmin,
    #                     bin2 <= xmax)
    # }
    p <- data %>% 
        ggplot(aes(x = bin1, y = bin2)) + 
            geom_raster(aes(fill = count)) +
	    scale_x_continuous(lim = c(minLim, maxLim), expand = c(0, 0)) +
	    scale_y_reverse(lim = c(maxLim, minLim), expand = c(0, 0)) +
            theme_bw() +
            theme(panel.spacing = unit(0, "lines")) +
            ggtitle(object@chromosome)
    if (logColor) {
        p <- p + scale_fill_gradient(low = "grey90", high = "red", trans = "log")
    } else {
        p <- p + scale_fill_gradient2()
        #p <- p + scale_fill_gradient(low = "blue", high = "red")
    }
    if (length(bins) > 100) {
        bins <- NULL
    }
    if (! is_null(bins)) {
        for (bin in bins) {
            p <- p +
                geom_hline(yintercept = bin, linetype = "dotted") +
                geom_vline(xintercept = bin, linetype = "dotted")
        }
    }    
    return(p)
}

plot.10X2Ref <- function(object,
			 logColor = TRUE,
			 circles = FALSE,
			 center  = NULL,
			 radius  = NULL) {
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
	    coord_fixed()
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

plot.10XDataset <- function(object, sizes, logColor = TRUE, ref1 = NULL, ref2 = NULL) {
    if (! is(object, "tenxcheckerExp")) {
        stop("Object should be a 'tenxcheckerExp'.")
    }
    if (! is.null(ref1)) {
        if (! is.null(ref2)) {
            object <- extract2Ref(object, ref1, ref2, sizes[[ref1]], sizes[[ref2]])
            return(plot.10X2Ref(object))
        }
        object <- extractRef(object, ref1, sizes[[ref1]])
        return(plot.10XRef(object))
    }
    scaleFactor <- computeScaleFactor(object, sizes)
    message(paste0("Scale factor: ", scaleFactor))
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
        ggtitle(object@name)
    if (logColor) {
	p <- p + scale_fill_gradient(low = "grey90", high = "red", trans = "log")
    }
    else {
	p <- p + scale_fill_gradient(low = "blue", high = "red")
    }
    return(p)
}

plot.10X <- function(object, sizes = NULL, logColor = TRUE, dataset = NULL, ref1 = NULL, ref2 = NULL) {
    if (is(object, "tenxcheckerClass")) {
        if (! is.null(ref1)) {
            if (! ref1 %in% object@chromosomes) {
                stop(paste0("Reference #1 '", ref1, "', is not a known reference."))
            }
            if (! is.null(ref2)) {
                if (! ref2 %in% object@chromosomes) {
                    stop(paste0("Reference #2 '", ref2, "', is not a known reference."))
                }
                if (match(ref1, object@chromosomes) < match(ref2, object@chromosomes)) {
                   tmp <- ref1
                   ref1 <- ref2
                   ref2 <- tmp
                }
            }
        }
        sizes  <- object@sizes
        if (is.null(dataset)) {
             plots <- purrr::map(object@data, plot.10XDataset, sizes, logColor, ref1, ref2)
             return(do.call("plot_grid", c(plots, ncol = length(plots))))
        }
        else {
            object <- getDataset(object, dataset)
        }
    }
    else if (! is(object, "tenxcheckerExp")) {
        stop("Object should be a 'tenxcheckerClass' or 'tenxcheckerExp'.")
    }
    return(plot.10XDataset(object, sizes, logColor, ref1, ref2))
}
