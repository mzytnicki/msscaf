estimateMoleculeSize <- function(object) {
    object@parameters@maxLinkRange <- object@interactionMatrix %>%
        filter(ref1 == ref2) %>%
        mutate(distance = abs(bin1 - bin2)) %>%
        dplyr::select(distance, count) %>%
        sample_n(min(object@parameters@sampleSize, nrow(.))) %>%
        mutate(loess = predict(loess(count ~ distance, data = ., span = 0.1))) %>%
        dplyr::select(distance, loess) %>%
        distinct() %>%
        arrange(distance) %>%
        filter(loess > object@parameters@minCount) %>%
        tail(n = 1) %>%
        pull(distance)
    message(paste0("Estimated molecule size: ", object@parameters@maxLinkRange, "."))
    return(object)
}

estimateBackgroundCounts <- function(object) {
    d <- object@interactionMatrix %>%
        filter(ref1 != ref2) %>%
        dplyr::select(count) %>%
        sample_n(min(object@parameters@sampleSize, nrow(.)))
    t <- transform(table(d), cum_freq = cumsum(Freq)) %>%
        mutate(relative = cum_freq / object@parameters@sampleSize) %>%
        filter(relative >= 0.5) %>%
        head(1) %>%
        pull(d)
    object@parameters@minCount <- as.integer(levels(t))[t] + 1
    message(paste0("Estimated background count: ", t, "."))
    return(object)
}

estimateMinRowCount <- function(object) {
    tmp <- object@interactionMatrix %>%
        makeSymmetric() %>%
        group_by(ref1, bin1) %>%
        summarise(countSum = sum(count)) %>%
        ungroup() %>%
        pull(countSum)
    tryCatch({
            f <- fitdistr(tmp, densfun="logistic")
            t <- max(object@parameters@minCount,
                     f$estimate["location"] + f$estimate["scale"] * qlogis(0.01))
            message(paste0("Estimated min. row count: ", t, "."))
            object@parameters@minRowCount <- t
            return(object)
        },
        error = function(e) {
            message(paste0("Cannot estimate min. row count, keeping default ", object@parameters@minCount, "."))
            object@parameters@minRowCount <- object@parameters@minCount
            return(object)
        }
    )
}

estimateDistributions <- function(object) {
    object <- estimateBackgroundCounts(object)
    if (is.null(object@parameters@maxLinkRange)) {
        object <- estimateMoleculeSize(object)
    }
    object@parameters@breakNCells <- 0.75 * ((object@parameters@maxLinkRange * (object@parameters@maxLinkRange + 1)) / 2)
    object <- estimateMinRowCount(object)
    return(object)
}
