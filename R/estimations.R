.estimateBackgroundCounts <- function(object) {
    if (! is(object, "tenxcheckerExp")) {
        stop("Parameter should be a tenxcheckerExp.")
    }
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
    message(paste0("Dataset '", object@name, "': Estimated background count: ", t, "."))
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
    object@parameters@breakNCells <- 0.75 * ((object@parameters@maxLinkRange * (object@parameters@maxLinkRange + 1)) / 2)
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
    object <- estimateMinRowCount(object, object@sizes)
    return(invisible(object))
}
