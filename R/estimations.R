estimateMoleculeSize <- function(object) {
    #backgroundTopFrac = 0.75
    d <- object@interactionMatrix %>%
        filter(ref1 == ref2) %>%
        mutate(distance = abs(bin1 - bin2)) %>%
        select(distance, count) %>%
        sample_n(min(object@parameters@sampleSize, nrow(.))) %>%
        mutate(loess = predict(loess(count ~ distance, data = ., span = 0.1)))
    background <- object@parameters@minCount
    # background <- object@interactionMatrix %>%
    #     filter(ref1 != ref2) %>%
    #     select(count) %>%
    #     sample_n(min(sampleSize, nrow(.))) %>%
    #     top_frac(backgroundTopFrac, count) %>%
    #     arrange(count) %>%
    #     head(n = 1) %>%
    #     pull(count)
    # backgroundPlot <- object@interactionMatrix %>%
    #     filter(ref1 != ref2) %>%
    #     select(count) %>%
    #     ggplot(aes(x = count)) +
    #     geom_freqpoly() +
    #     geom_vline(xintercept = background, linetype = "dashed") +
    #     scale_y_log10()
    distance <- d %>%
        select(distance, loess) %>%
        distinct() %>%
        arrange(distance) %>%
        filter(loess > background) %>%
        tail(n = 1) %>%
        pull(distance)
    p <- ggplot(d, aes(distance, count)) +
        geom_point(color = "grey50") +
        geom_line(aes(x = distance, y = loess)) + 
        geom_hline(yintercept = background, linetype = "dashed") +
        geom_vline(xintercept = distance, linetype = "dashed") +
        scale_y_log10()
    message(paste0("Estimated molecule size: ", distance, "."))
    return(list(size = distance, plot = p))
    #return(list(size = distance, plot = p, backgroundPlot = backgroundPlot))
}

estimateBackgroundCounts <- function(object) {
    d <- object@interactionMatrix %>%
        filter(ref1 != ref2) %>%
        #filter((ref1 != ref2) | (abs(bin1 - bin2) > object@parameters@maxLinkRange)) %>%
        select(count) %>%
        sample_n(min(object@parameters@sampleSize, nrow(.)))
    t <- transform(table(d), cum_freq = cumsum(Freq)) %>%
        mutate(relative = cum_freq / object@parameters@sampleSize) %>%
        filter(relative >= 0.5) %>%
        head(1) %>%
        pull(d)
    t <- as.integer(levels(t))[t] + 1
    p <- ggplot(d, aes(count)) +
        # geom_freqpoly() +
        geom_freqpoly(binwidth=1) +
        geom_vline(xintercept = t, linetype = "dashed") +
        xlim(1, max(30, 2 * t)) +
        scale_y_log10()
    message(paste0("Estimated background count: ", t, "."))
    return(list(count = t, plot = p))
}

estimateMinRowCount <- function(object) {
    xWidth <- 100000
    tmp <- object@interactionMatrix %>%
        makeSymmetric() %>%
        group_by(ref1, bin1) %>%
        summarise(countSum = sum(count)) %>%
        ungroup() %>%
        select(countSum)
    # f <- fitdistr(tmp$countSum, densfun="negative binomial")
    # p <- ggplot(tmp, aes(countSum)) +
    #     geom_density() + xlim(0, xWidth) +
    #     stat_function(n = xWidth+1,
    #                   fun = dnbinom,
    #                   args = list(size = f$estimate["size"], mu = f$estimate["mu"]),
    #                   xlim = c(0, xWidth),
    #                   color = "red") +
    #     geom_vline(xintercept = t, linetype = "dashed")
    # p
    t <- object@parameters@minCount
    p <- NULL
    tryCatch({
            f <- fitdistr(tmp$countSum, densfun="logistic")
            t <- max(object@parameters@minCount,
                     f$estimate["location"] + f$estimate["scale"] * qlogis(0.01))
            p <- ggplot(tmp, aes(countSum)) +
                geom_density() + xlim(0, xWidth) +
                stat_function(n = xWidth+1,
                              fun = dlogis,
                              args = list(location = f$estimate["location"], scale = f$estimate["scale"]),
                              xlim = c(0, xWidth),
                              color = "red") +
                geom_vline(xintercept = t, linetype = "dashed")
            message(paste0("Estimated min. row count: ", t, "."))
           return(list(count = t, plot = p))
        },
        error = function(e) {
            p <- ggplot(tmp, aes(countSum)) +
                geom_density() + xlim(0, xWidth) +
                geom_vline(xintercept = t, linetype = "dashed")
            message(paste0("Cannot estimate min. row count, keeping default ", t, "."))
           return(list(count = t, plot = p))
        }
    )
}

estimateDistributions <- function(object) {
    l <- estimateBackgroundCounts(object)
    object@parameters@minCount <- l$count
    if (is.null(object@parameters@maxLinkRange)) {
        l <- estimateMoleculeSize(object)
        object@parameters@maxLinkRange <- l$size
    }
    object@parameters@breakNCells <- 0.75 * ((object@parameters@maxLinkRange * (object@parameters@maxLinkRange + 1)) / 2)
    l <- estimateMinRowCount(object)
    object@parameters@minRowCount <- l$count
    return(object)
}
