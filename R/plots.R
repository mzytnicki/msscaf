plotScaffoldSizeDistribution <- function(object, log) {
    p <- object@interactionMatrix %>%
        makeSymmetric() %>%
        group_by(ref1) %>%
        summarise(size = max(bin1)) %>%
        ungroup() %>%
        ggplot(aes(x = ref1, y = size)) + 
            geom_col()
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
            geom_bin2d(binwidth = c(1, 10))
    if (log) {
        p <- p + scale_x_log10() + scale_y_log10()
    }
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

plot.10X <- function(object, logColor = TRUE) {
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

  
plot.10XRef <- function(object, logColor = TRUE, split = FALSE) {
    p <- object@interactionMatrix %>%
        makeSymmetric() %>%
        ggplot(aes(x = bin1, y = bin2)) + 
            geom_raster(aes(fill = count)) + 
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

  