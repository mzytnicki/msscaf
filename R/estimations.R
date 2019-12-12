estimateMoleculeSize <- function(object) {
    sampleSize = 10000
    d <- object@interactionMatrix %>%
        filter(ref1 == ref2) %>%
        mutate(distance = abs(bin1 - bin2)) %>%
        select(distance, count) %>%
        sample_n(min(sampleSize, nrow(.))) %>%
        mutate(loess = predict(loess(count ~ distance, data = ., span = 0.1)))
    distance <- d %>%
        select(distance, loess) %>%
        distinct() %>%
        arrange(distance) %>%
        filter(loess >= object@parameters@minCount) %>%
        tail(n = 1) %>%
        pull(distance)
    p <- ggplot(d, aes(distance, count)) +
        geom_point(color = "grey50") +
        geom_line(aes(x = distance, y = loess)) + xlim(0, 50) +
        geom_vline(xintercept = distance, linetype = "dashed")
    message(paste0("Estimated molecule size: ", distance, "."))
    return(list(size = distance, plot = p))
}

estimateBackgroundCounts <- function(object) {
    sampleSize <- 100000
    d <- object@interactionMatrix %>%
        filter((ref1 != ref2) | (abs(bin1 - bin2) > object@parameters@maxLinkRange)) %>%
        select(count)
    sampleSize <- min(sampleSize, nrow(d))
    d %<>% sample_n(sampleSize)
    t <- transform(table(d), cum_freq = cumsum(Freq)) %>%
        mutate(relative = cum_freq / sampleSize) %>%
        filter(relative >= 0.9) %>%
        head(1) %>%
        pull(d)
    t <- as.integer(levels(t))[t] + 1
    p <- ggplot(d, aes(count)) +
        geom_histogram(binwidth=1) + xlim(0, 10) +
        geom_vline(xintercept = t, linetype = "dashed")
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
    f <- fitdistr(tmp$countSum, densfun="logistic")
    t <- f$estimate["location"] + f$estimate["scale"] * qlogis(0.01)
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
}

estimateDistributions <- function(object) {
    l <- estimateMoleculeSize(object)
    object@parameters@maxLinkRange <- l$size
    l <- estimateBackgroundCounts(object)
    object@parameters@minCount <- l$count
    l <- estimateMinRowCount(object)
    object@parameters@minRowCount <- l$count
    return(object)
}