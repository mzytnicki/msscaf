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
        # Rearrange into bin1/distance matrix
        diagMatrices <- object@interactionMatrix %>%
            dplyr::filter(ref2 == ref1) %>%
            dplyr::mutate(distance = bin1 - bin2) %>%
            dplyr::arrange(ref1, bin1, distance)
        # A block is a contiguous stretch of non-zero cells
        blocks <- diagMatrices %>%
            dplyr::group_by(ref1, bin1) %>%
            # compute groups of contiguous values
            dplyr::mutate(diffDist = distance - dplyr::lag(distance) - 1) %>%
            tidyr::replace_na(list(diffDist = 0)) %>%
            dplyr::mutate(diffDist = if_else(diffDist > 1, 1, diffDist)) %>%
            dplyr::mutate(group = cumsum(diffDist)) %>%
            # find start/end/size
            dplyr::group_by(ref1, bin1, group) %>%
            dplyr::summarise(groupStart = min(distance), groupEnd = max(distance), .groups = "drop") %>%
            dplyr::select(-group) %>%
            dplyr::mutate(groupSize = as.integer(groupEnd - groupStart + 1))
        # A hole is a contiguous stretch of zero cells
        holes <- diagMatrices %>%
            dplyr::group_by(ref1, bin1) %>%
            dplyr::mutate(diffDist = as.integer(distance - lag(distance))) %>%
            dplyr::slice(-1) %>%
            dplyr::ungroup() %>%
            dplyr::filter(diffDist != 1) %>%
            dplyr::mutate(groupEnd = distance - diffDist, groupStart = distance) %>%
            dplyr::mutate(holeSize = as.integer(diffDist - 1)) %>%
            # compute size, and end/start positions of flanking groups
            dplyr::select(ref1, bin1, groupEnd, groupStart, holeSize)
        # An empty diagonal is a hole which starts at the diagonal, and is larger than the first block
        notEmptyDiag <- blocks %>%
            dplyr::group_by(ref1, bin1) %>%
            dplyr::arrange(groupStart) %>%
            dplyr::summarise(groupStart = dplyr::first(groupStart), groupSize = dplyr::first(groupSize), .groups = "drop") %>%
            dplyr::filter(groupSize >= groupStart) %>%
            dplyr::select(ref1, bin1)
        # A true hole should have size greater than flanking groups
        firstHoles <- holes %>%
            dplyr::left_join(blocks, by = c("ref1", "bin1", "groupEnd"), suffix = c("", "_group")) %>%
            dplyr::rename(groupSizeLeft = groupSize) %>%
            dplyr::select(ref1, bin1, groupEnd, groupStart, holeSize, groupSizeLeft) %>%
            dplyr::left_join(blocks, by = c("ref1", "bin1", "groupStart"), suffix = c("", "_group")) %>%
            dplyr::rename(groupSizeRight = groupSize) %>%
            dplyr::select(ref1, bin1, groupEnd, groupSizeLeft, holeSize, groupSizeRight) %>%
            dplyr::filter(holeSize > groupSizeLeft & holeSize > groupSizeRight) %>%
            # take the first true hole
            dplyr::group_by(ref1, bin1) %>%
            dplyr::summarise(groupEnd = min(groupEnd), .groups = "drop") 
        # Compute the most distant point
        lastPoints <- diagMatrices %>%
            dplyr::group_by(ref1, bin1) %>%
            dplyr::summarise(last = max(distance), .groups = "drop")
        # The molecule should be at the beginning of the first true hole.
        #    If it does not exist, it is the most distant point.
        #    Empty diagonals are discarded.
        object@parameters@maxLinkRange <- notEmptyDiag %>%
            dplyr::left_join(lastPoints, by = c("ref1", "bin1")) %>%
            dplyr::left_join(firstHoles, by = c("ref1", "bin1")) %>%
            dplyr::mutate(size = dplyr::if_else(is.na(groupEnd), last, groupEnd)) %>%
            dplyr::slice_max(size, prop = 0.1) %>%
            dplyr::arrange(size) %>%
            dplyr::slice(1) %>%
            dplyr::pull(size)
#       object@parameters@maxLinkRange <- object@interactionMatrix %>%
#           filter(ref1 == ref2) %>%
#           mutate(distance = abs(bin1 - bin2)) %>%
#           dplyr::select(distance, count) %>%
#           sample_n(min(object@parameters@sampleSize, nrow(.))) %>%
#           group_by(distance) %>%
#           summarise(count = sum(count)) %>%
#           arrange(distance) %>%
#           mutate(cumulated = cumsum(count)) %>%
#           mutate(cumulated = cumulated / sum(.$count)) %>%
#           filter(cumulated <= 0.5) %>%
#           mutate(loess = predict(loess(count ~ distance, data = ., span = 0.1))) %>%
#           dplyr::select(distance, loess) %>%
#           distinct() %>%
#           arrange(distance) %>%
#           filter(loess > object@parameters@minCount) %>%
#           tail(n = 1) %>%
#           pull(distance)
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

estimateDistanceCount <- function(object, sizes) {
    if (! is(object, "tenxcheckerExp")) {
        stop("Parameter should be a tenxcheckerExp.")
    }
    message("\t\tEstimating distance/count distribution.")
    # Set genome to (ref1, bin1, distance), and transform missing values to 0
    if (sum(sizes) <= object@parameters@sampleSize) {
        distanceCount <- object@interactionMatrix %>%
            dplyr::filter(ref1 == ref2) %>%
            dplyr::mutate(distance = bin1 - bin2) %>%
            dplyr::filter(distance <= object@parameters@maxLinkRange) %>%
            dplyr::select(ref1, bin1, distance, count) %>%
            tidyr::complete(nesting(ref1, bin1), distance, fill = list(count = 0)) %>%
            dplyr::rename(ref = ref1, bin = bin1)
    }
    # Cannot do the loess on the whole genome (too big): sample a few positions
    else {
        sampleSize <- max(100, object@parameters@sampleSize / (object@parameters@maxLinkRange + 1))
        lines <- sizes %>%
            tibble::enframe(name = "ref", value = "size") %>%
            dplyr::mutate(size = size - object@parameters@maxLinkRange) %>%
            dplyr::filter(size > 0) %>%
            dplyr::sample_n(sampleSize, replace = TRUE, weigth = size) %>%
            dplyr::mutate(pos = runif(nrow(.))) %>%
            dplyr::mutate(bin = as.integer(pos * size)) %>%
            dplyr::select(ref, bin) %>%
            dplyr::mutate(ref = factor(ref, levels = names(sizes)))
        distanceCount <- extractLines(object@interactionMatrix, lines, object@parameters@maxLinkRange) %>%
            as_tibble()
#message(str(distanceCount))
    }
    # Perform the loess
    object@parameters@distanceCount <- distanceCount %>%
        # replace loess with mean?
        dplyr::mutate(count = predict(loess(count ~ distance, data = .))) %>%
        dplyr::select(distance, count) %>%
        # dplyr::group_by(distance) %>%
        # dplyr::summarise(count = median(count), .groups = "drop") %>%
        dplyr::arrange(distance) %>%
        dplyr::distinct()
    return(object)
}

compareCornerOffset <- function(offset, corner, background, maxDistance) {
    cornerDistance(corner     %>%
                       dplyr::mutate(distance = distance - offset) %>%
                       dplyr::filter(distance >= 0),
                   background %>%
                       dplyr::mutate(distance = distance - offset) %>%
                       dplyr::filter(distance >= 0))
#   cornerEuclideanDistance(corner     %>% dplyr::filter(distance >= offset),
#                           background %>% dplyr::filter(distance >= offset))
}

compareCorner <- function(corner, background, distance, pb) {
# message("corner")
# message(str(corner))
# message("background")
# message(str(background))
    pb$tick()
    corner <- corner %>%
        # replace loess with mean?
        dplyr::mutate(count = predict(loess(count ~ distance, data = .))) %>%
        dplyr::select(distance, count) %>%
        # dplyr::group_by(distance) %>%
        # dplyr::summarise(count = median(count), .groups = "drop") %>%
        dplyr::arrange(distance) %>%
        dplyr::distinct()
    purrr::map(seq.int(from = 0, to = distance - 1, by = 1), compareCornerOffset, corner, background, distance)
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
    fillCorner(corner, object@parameters@maxLinkRange, FALSE)
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

estimateCornerVariance <- function(object, sizes, pvalueThreshold) {
    if (! is(object, "tenxcheckerExp")) {
        stop("Parameter should be a tenxcheckerExp.")
    }
    nSamples <- 10 / pvalueThreshold
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
    return(object)
}

# Find offset where corner detection is not significant
estimateCornerLimits <- function(object, minNBins) {
    return(object@parameters@maxLinkRange)
#   bins <- seq.int(from = 0, to = object@parameters@maxLinkRange, by = 1)
#   emptyCorner <- tibble(bin1 = bins, bin2 = bins) %>%
#       tidyr::expand(bin1, bin2) %>%
#       dplyr::filter(bin1 >= bin2) %>%
#       dplyr::mutate(distance = bin1 - bin2) %>%
#       dplyr::filter(distance <= object@parameters@maxLinkRange) %>%
#       dplyr::filter(bin1 <= minNBins) %>%
#       dplyr::filter(bin2 <= minNBins) %>%
#       dplyr::mutate(count = 0) %>%
#       dplyr::select(distance, count)
#   values <- purrr::map(bins, computeCornerDistanceOffset, object@parameters@distanceCount, emptyCorner, object@parameters@maxLinkRange) %>% unlist()
#   object@parameters@cornerScores %>%
#       dplyr::mutate(observed = values) %>%
#       dplyr::mutate(difference = observed - score) %>%
#       dplyr::filter(difference < 0.0) %>%
#       dplyr::slice_min(order_by = distance, n = 1) %>%
#       dplyr::pull(distance)
}

.estimateCorners <- function(object, sizes, pvalueThreshold, minNBins) {
    if (! is(object, "tenxcheckerExp")) {
        stop("Parameter should be a tenxcheckerExp.")
    }
    message(paste0("\tDataset '", object@name, "':"))
    object   <- estimateDistanceCount(object, sizes)
    object   <- estimateCornerVariance(object, sizes, pvalueThreshold)
    object@parameters@cornerLimit <- estimateCornerLimits(object, minNBins)
    return(object)
}

estimateCorners <- function(object, pvalueThreshold) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    message("Join shape estimations.")
    object@data <- purrr::map(object@data, .estimateCorners, object@sizes, pvalueThreshold, object@minNBins)
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
