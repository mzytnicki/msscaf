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
        select(ref1, bin1, distance, count) %>%
        arrange(ref1, bin1, distance)
    tmp1 <- tmp %>%
        group_by(ref1, bin1) %>%
        mutate(diagonalSum = cumsum(count)) %>%
        ungroup() %>%
        select(ref1, bin1, distance, diagonalSum)
    p <- tmp %>%
        group_by(ref1, bin1) %>%
        summarise(countSum = sum(count)) %>%
        ungroup() %>%
        right_join(tmp1, by = c("ref1", "bin1")) %>%
        mutate(diagPc = diagonalSum / countSum) %>%
        filter(distance <= 2 * object@parameters@maxLinkRange) %>%
        select(ref1, bin1, distance, diagPc) %>%
        rename(ref = ref1) %>%
        rename(bin = bin1) %>%
        group_by(ref)
    n <- p %>% group_keys()
    p <- p %>%
        group_split(keep = FALSE) %>%
        map(~ ggplot(., aes(x = bin, y = distance, fill = diagPc)) + 
                geom_tile())
    names(p) <- n$ref
    return(p)
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

plot.10XRef <- function(object, logColor = TRUE, bin = NA) {
    tmp <- object@interactionMatrix %>%
        makeSymmetric()
    if (!is.na(bin)) {
        xmin <- max(0, bin - object@parameters@nBinZoom)
        xmax <- min(object@size, bin + object@parameters@nBinZoom)
        tmp %<>% filter(bin1 >= xmin,
                        bin1 <= xmax,
                        bin2 >= xmin,
                        bin2 <= xmax)
    }
    p <- tmp %>% 
        ggplot(aes(x = bin1, y = bin2)) + 
            geom_raster(aes(fill = count)) +
            scale_x_continuous(expand = c(0, 0)) +
            scale_y_reverse(expand = c(0, 0)) +
            theme_bw() +
            theme(panel.spacing = unit(0, "lines")) +
            ggtitle(object@chromosome)
    if (logColor) {
        p <- p + scale_fill_gradient(low = "grey90", high = "red", trans = "log")
    } else {
        p <- p + scale_fill_gradient2()
        #p <- p + scale_fill_gradient(low = "blue", high = "red")
    }
    if (!is.na(bin)) {
        p <- p +
            geom_hline(yintercept = bin, linetype = "dotted") +
            geom_vline(xintercept = bin, linetype = "dotted")
    }
    return(p)
}

plot.10X2Ref <- function(object,
                         logColor = TRUE,
                         x1 = NULL,
                         x2 = NULL,
                         y1 = NULL,
                         y2 = NULL) {
    p <- object@interactionMatrix %>%
        ggplot(aes(x = bin1, y = bin2)) + 
            geom_raster(aes(fill = count)) + 
            scale_x_continuous(expand = c(0, 0)) +
            scale_y_reverse(expand = c(0, 0)) +
            xlab(object@chromosome1) +
            ylab(object@chromosome2) +
            theme_bw() +
            theme(panel.spacing = unit(0, "lines"))
    if (!is.null(x1)) {
        p <- p + geom_vline(xintercept = x1, linetype = "dotted", color = "red")
    }
    if (!is.null(x2)) {
        p <- p + geom_vline(xintercept = x2, linetype = "dashed", color = "red")
    }
    if (!is.null(y1)) {
        p <- p + geom_hline(yintercept = y1, linetype = "dotted", color = "red")
    }
    if (!is.null(y2)) {
        p <- p + geom_hline(yintercept = y2, linetype = "dashed", color = "red")
    }
    if (logColor) {
        p <- p + scale_fill_gradient(low = "grey90", high = "red", trans = "log")
    }
    else {
        p <- p + scale_fill_gradient(low = "blue", high = "red")
    }
    return(p)
}


plot.10X <- function(object, logColor = TRUE, ref = NULL) {
    if (!is.null(ref)) {
        p <- object@interactionMatrix %>%
            filter(ref1 == ref) %>%
            filter(ref2 == ref) %>%
            tenxcheckerRefExp(ref, object@parameters) %>%
            plot.10XRef()
        return(p)
    }
    p <- object@interactionMatrix %>%
        makeSymmetric() %>%
        ggplot(aes(x = bin1, y = bin2)) + 
        geom_raster(aes(fill = count)) + 
        facet_grid(cols = vars(ref1), rows = vars(ref2), scale = "free", space = "free") +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_reverse(expand = c(0, 0)) +
        theme_bw() +
        theme(panel.spacing = unit(0, "lines"))
    if (logColor) {
        p <- p + scale_fill_gradient(low = "grey90", high = "red", trans = "log")
    }
    else {
        p <- p + scale_fill_gradient(low = "blue", high = "red")
    }
    return(p)
}
