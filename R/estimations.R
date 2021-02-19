.estimateBackgroundCounts <- function(object) {
    if (! is(object, "tenxcheckerExp")) {
        stop("Parameter should be a tenxcheckerExp.")
    }
    d <- object@interactionMatrix %>%
        filter(ref1 != ref2) %>%
        dplyr::select(count)
    if (nrow(d) == 0) {
        message(paste0("Dataset '", object@name, "': Cannot estimate background count: there is no count outside diagonal matrices.\n\tSetting it as 1."))
        object@parameters@minCount <- 1
        return(object)
    }
    sampleSize <- min(object@parameters@sampleSize, nrow(d))
    d <- d %>%
        sample_n(sampleSize)
    threshold <- transform(table(d), cum_freq = cumsum(Freq)) %>%
        mutate(relative = cum_freq / sampleSize) %>%
        filter(relative >= 0.5) %>%
        head(1) %>%
        pull(d)
    threshold <- as.integer(levels(threshold))[threshold]
    object@parameters@minCount <- threshold
    message(paste0("Dataset '", object@name, "': Estimated background count: ", threshold, "."))
    return(object)
}

estimateBackgroundCounts <- function(object) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    object@data <- purrr::map(object@data, .estimateBackgroundCounts)
    return(object)
}

.estimateMoleculeSize <- function(object) {
    if (! is(object, "tenxcheckerExp")) {
        stop("Parameter should be a tenxcheckerExp.")
    }
    if (is.null(object@parameters@maxLinkRange)) {
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
        message(paste0("Dataset '", object@name, "': Estimated molecule size: ", object@parameters@maxLinkRange, "."))
    }
    return(object)
}

estimateMoleculeSize <- function(object) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    object@data <- purrr::map(object@data, .estimateMoleculeSize)
    return(object)
}

.estimateMinRowCount <- function(object, sizes) {
    if (! is(object, "tenxcheckerExp")) {
        stop("Parameter should be a tenxcheckerExp.")
    }
    tmp <- computeSymmetricColSum(object@interactionMatrix, sizes)
    tryCatch({
            f <- fitdistr(tmp, densfun="logistic")
            t <- max(object@parameters@minCount,
                     f$estimate["location"] + f$estimate["scale"] * qlogis(0.01))
            message(paste0("Dataset '", object@name, "': Estimated min. row count: ", t, "."))
            object@parameters@minRowCount <- t
            return(object)
        },
        error = function(e) {
            message(paste0("Dataset '", object@name, "': Cannot estimate min. row count, keeping default ", object@parameters@minCount, "."))
            object@parameters@minRowCount <- object@parameters@minCount
            return(object)
        }
    )
}

estimateMinRowCount <- function(object, sizes) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    object@data <- purrr::map(object@data, .estimateMinRowCount, sizes = sizes)
    return(object)
}

estimateDistributions <- function(object) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    object <- estimateBackgroundCounts(object)
    object <- estimateMoleculeSize(object)
    return(invisible(object))
}
