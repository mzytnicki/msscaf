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
            group_by(distance) %>%
            summarise(count = sum(count)) %>%
            arrange(distance) %>%
            mutate(cumulated = cumsum(count)) %>%
            mutate(cumulated = cumulated / sum(.$count)) %>%
            filter(cumulated <= 0.5) %>%
#           mutate(loess = predict(loess(count ~ distance, data = ., span = 0.1))) %>%
#           dplyr::select(distance, loess) %>%
#           distinct() %>%
#           arrange(distance) %>%
#           filter(loess > object@parameters@minCount) %>%
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

estimateDistanceCount <- function(object) {
    if (! is(object, "tenxcheckerExp")) {
        stop("Parameter should be a tenxcheckerExp.")
    }
    message("\t\tEstimating distance/count distribution.")
    # Do not forget to add empty cells: zero count
    diagonal <- object@interactionMatrix %>%
        dplyr::filter(ref1 == ref2) %>%
        dplyr::mutate(distance = bin1 - bin2) %>%
        dplyr::select(ref1, bin1, distance, count) %>%
        dplyr::filter(distance <= object@parameters@maxLinkRange)
    object@parameters@distanceCount <- diagonal %>%
        tidyr::expand(nesting(ref1, bin1), distance = seq.int(from = 0, to = object@parameters@maxLinkRange, by = 1)) %>%
        dplyr::left_join(diagonal, by = c("ref1", "bin1", "distance")) %>%
        tidyr::replace_na(list(count = 0)) %>%
        dplyr::group_by(distance) %>%
        dplyr::summarise(meanCount = mean(count))
    return(object)
}

compareCornerOffset <- function(offset, corner, background) {
    cornerEuclideanDistance(corner     %>% dplyr::filter(distance >= offset),
                            background %>% dplyr::filter(distance >= offset))
}

compareCorner <- function(corner, background, distance, pb) {
    pb$tick()
    purrr::map(seq.int(from = 0, to = distance, by = 1), compareCornerOffset, corner, background)
}

extractCornerFromPoint <- function(parameters, object, pb) {
    pb$tick()
    corner <- object@interactionMatrix %>%
        dplyr::filter(ref1 == parameters$ref1) %>%
        dplyr::filter(ref2 == parameters$ref2) %>%
        dplyr::mutate(bin1 = bin1 - parameters$bin1) %>%
        dplyr::mutate(bin2 = bin2 - parameters$bin2) %>%
        dplyr::filter(bin1 >= 0, bin1 <= object@parameters@maxLinkRange) %>%
        dplyr::filter(bin2 >= 0, bin2 <= object@parameters@maxLinkRange) %>%
        dplyr::filter(bin1 >= bin2)
    computeDistanceCount(corner, object@parameters@maxLinkRange)
}

findRandomCornerPoints <- function(object, sizes, nSamples) {
    # Extract all observed pairs of refs
    object@interactionMatrix %>%
        dplyr::select(ref1, ref2) %>%
        dplyr::filter(ref1 != ref2) %>%
        dplyr::distinct() %>%
    # Compute ref sizes
        dplyr::mutate(size1 = sizes[ref1]) %>%
        dplyr::mutate(size2 = sizes[ref2]) %>%
        dplyr::mutate(size  = size1 * size2) %>%
    # Sample, weighted by size
        dplyr::slice_sample(n = nSamples, weight_by = size, replace = TRUE) %>%
    # Get random point
        dplyr::mutate(size1 = size1 - object@parameters@maxLinkRange) %>%
        dplyr::mutate(size2 = size2 - object@parameters@maxLinkRange) %>%
        dplyr::mutate(bin1 = as.integer(round(runif(nrow(.)) * size1))) %>%
        dplyr::mutate(bin2 = as.integer(round(runif(nrow(.)) * size2))) %>%
        dplyr::select(ref1, bin1, ref2, bin2)
}

# Find offset where corner detection is not significant
estimateCornerLimits <- function(object, minNBins) {
    bins <- seq.int(from = 0, to = object@parameters@maxLinkRange, by = 1)
    emptyCorner <- tibble(bin1 = bins, bin2 = bins) %>%
        tidyr::expand(bin1, bin2) %>%
        dplyr::filter(bin1 >= bin2) %>%
        dplyr::mutate(distance = bin1 - bin2) %>%
        dplyr::filter(distance <= object@parameters@maxLinkRange) %>%
        dplyr::filter(bin1 <= minNBins) %>%
        dplyr::filter(bin2 <= minNBins) %>%
        dplyr::mutate(count = 0) %>%
        dplyr::select(distance, count)
    values <- purrr::map(bins, computeCornerDistanceOffset, object@parameters@distanceCount, emptyCorner, object@parameters@maxLinkRange) %>% unlist()
message(str(object@parameters@maxLinkRange))
message(str(object@parameters@cornerScores))
message(str(values))
    object@parameters@cornerScores %>%
        dplyr::mutate(observed = values) %>%
        dplyr::mutate(difference = observed - score) %>%
        dplyr::filter(difference < 0.0) %>%
        dplyr::slice_min(order_by = distance, n = 1) %>%
        dplyr::pull(distance)
}

.estimateCornerDistanceThreshold <- function(object, sizes, pvalueThreshold, minNBins) {
    if (! is(object, "tenxcheckerExp")) {
        stop("Parameter should be a tenxcheckerExp.")
    }
    nSamples <- 10 / pvalueThreshold
    message(paste0("\tDataset '", object@name, "':"))
    object   <- estimateDistanceCount(object)
    points   <- findRandomCornerPoints(object, sizes, nSamples) %>%
        purrr::transpose()
    message("\t\tExtracting random matrices.")
    pb       <- progress_bar$new(total = length(points))
    corners  <- purrr::map(points, extractCornerFromPoint, object, pb)
    message("\t\tComputing stats.")
    pb       <- progress_bar$new(total = length(points))
    object@parameters@cornerScores <- purrr::map(corners, compareCorner, object@parameters@distanceCount, object@parameters@maxLinkRange, pb) %>%
        purrr::transpose() %>%
        purrr::map(unlist) %>%
        purrr::map(sort) %>%
        purrr::map(nSamples * pvalueThreshold) %>%
        unlist() %>%
        tibble::enframe(name = "distance", value = "score") %>%
        dplyr::mutate(distance = distance - 1)
    object@parameters@cornerLimit <- estimateCornerLimits(object, minNBins)
    return(object)
}

estimateCornerDistanceThreshold <- function(object, pvalueThreshold) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    message("Join shape estimations.")
    object@data <- purrr::map(object@data, .estimateCornerDistanceThreshold, object@sizes, pvalueThreshold, object@minNBins)
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
