.estimateBackgroundCounts <- function(object) {
    if (! is(object, "msscafExp")) {
        stop("Parameter should be a msscafExp.")
    }
    d <- object@interactionMatrix %>%
        dplyr::filter(ref1 != ref2) %>%
        dplyr::select(count)
    if (nrow(d) == 0) {
        message(paste0("Dataset '", object@name, "': Cannot estimate background count: there is no count outside diagonal matrices.\n\tSetting it as 1."))
        object@parameters@minCount <- 1
        return(invisible(object))
    }
    sampleSize <- min(object@parameters@sampleSize, nrow(d))
    object@parameters@minCount <- d
        dplyr::sample_n(sampleSize) %>%
        table() %>%
        tibble::as_tibble() %>%
        dplyr::rename("count" = ".") %>%
        dplyr::mutate(count = as.integer(count)) %>%
        dplyr::mutate(cumCount = cumsum(n)) %>%
        dplyr::mutate(relative = cumCount / sampleSize) %>%
        dplyr::filter(relative >= 0.5) %>%
        dplyr::slice(1:1) %>%
        dplyr::pull(count)
    message(paste0("Dataset '", object@name, "': Estimated background count: ", threshold, "."))
    return(invisible(object))
}

estimateBackgroundCounts <- function(object) {
    if (! is(object, "msscafClass")) {
        stop("Parameter should be a msscafClass.")
    }
    object@data <- purrr::map(object@data, .estimateBackgroundCounts)
    return(invisible(object))
}

.estimateMoleculeSize <- function(object, sizes) {
    if (! is(object, "msscafExp")) {
        stop("Parameter should be a msscafExp.")
    }
    if (! is.null(object@parameters@maxLinkRange)) {
        return(invisible(object))
    }
    message(paste0(object@name, ":"))
    object@parameters@maxLinkRange <- estimateMoleculeSizeCpp(object@interactionMatrix, object@outlierBins, sizes, object@parameters@metaSize, object@parameters@minCount)
#   # Rearrange into bin1/distance matrix
#   diagMatrices <- object@interactionMatrix %>%
#       dplyr::filter(ref2 == ref1) %>%
#       dplyr::mutate(distance = bin1 - bin2) %>%
#       dplyr::arrange(ref1, bin1, distance)
#   # A block is a contiguous stretch of non-zero cells
#   blocks <- diagMatrices %>%
#       dplyr::group_by(ref1, bin1) %>%
#       # compute groups of contiguous values
#       dplyr::mutate(diffDist = distance - dplyr::lag(distance) - 1) %>%
#       tidyr::replace_na(list(diffDist = 0)) %>%
#       dplyr::mutate(diffDist = dplyr::if_else(diffDist > 1, 1, diffDist)) %>%
#       dplyr::mutate(group = cumsum(diffDist)) %>%
#       # find start/end/size
#       dplyr::group_by(ref1, bin1, group) %>%
#       dplyr::summarise(groupStart = min(distance), groupEnd = max(distance), .groups = "drop") %>%
#       dplyr::select(-group) %>%
#       dplyr::mutate(groupSize = as.integer(groupEnd - groupStart + 1))
#   # A hole is a contiguous stretch of zero cells
#   holes <- diagMatrices %>%
#       dplyr::group_by(ref1, bin1) %>%
#       dplyr::mutate(diffDist = as.integer(distance - lag(distance))) %>%
#       dplyr::slice(-1) %>%
#       dplyr::ungroup() %>%
#       dplyr::filter(diffDist != 1) %>%
#       dplyr::mutate(groupEnd = distance - diffDist, groupStart = distance) %>%
#       dplyr::mutate(holeSize = as.integer(diffDist - 1)) %>%
#       # compute size, and end/start positions of flanking groups
#       dplyr::select(ref1, bin1, groupEnd, groupStart, holeSize)
#   # An empty diagonal is a hole which starts at the diagonal, and is larger than the first block
#   notEmptyDiag <- blocks %>%
#       dplyr::group_by(ref1, bin1) %>%
#       dplyr::arrange(groupStart) %>%
#       dplyr::summarise(groupStart = dplyr::first(groupStart), groupSize = dplyr::first(groupSize), .groups = "drop") %>%
#       dplyr::filter(groupSize >= groupStart) %>%
#       dplyr::select(ref1, bin1)
#   # A true hole should have size greater than flanking groups
#   firstHoles <- holes %>%
#       dplyr::left_join(blocks, by = c("ref1", "bin1", "groupEnd"), suffix = c("", "_group")) %>%
#       dplyr::rename(groupSizeLeft = groupSize) %>%
#       dplyr::select(ref1, bin1, groupEnd, groupStart, holeSize, groupSizeLeft) %>%
#       dplyr::left_join(blocks, by = c("ref1", "bin1", "groupStart"), suffix = c("", "_group")) %>%
#       dplyr::rename(groupSizeRight = groupSize) %>%
#       dplyr::select(ref1, bin1, groupEnd, groupSizeLeft, holeSize, groupSizeRight) %>%
#       dplyr::filter(holeSize > groupSizeLeft & holeSize > groupSizeRight) %>%
#       # take the first true hole
#       dplyr::group_by(ref1, bin1) %>%
#       dplyr::summarise(groupEnd = min(groupEnd), .groups = "drop") 
#   # Compute the most distant point
#   lastPoints <- diagMatrices %>%
#       dplyr::group_by(ref1, bin1) %>%
#       dplyr::summarise(last = max(distance), .groups = "drop")
#   # The molecule should be at the beginning of the first true hole.
#   #    If it does not exist, it is the most distant point.
#   #    Empty diagonals are discarded.
#   object@parameters@maxLinkRange <- notEmptyDiag %>%
#       dplyr::left_join(lastPoints, by = c("ref1", "bin1")) %>%
#       dplyr::left_join(firstHoles, by = c("ref1", "bin1")) %>%
#       dplyr::mutate(size = dplyr::if_else(is.na(groupEnd), last, groupEnd)) %>%
#       dplyr::slice_max(size, prop = 0.1) %>%
#       dplyr::arrange(size) %>%
#       dplyr::slice(1) %>%
#       dplyr::pull(size)


#       selectedBins <- object@interactionMatrix %>%
#           dplyr::select(ref1, bin1) %>%
#           dplyr::distinct() %>%
#           dplyr::sample_n(min(object@parameters@sampleSize, nrow(.)))
#       selectedCells <- object@interactionMatrix %>%
#           dplyr::filter(ref1 == ref2) %>%
#           dplyr::right_join(selectedBins, by = c("ref1", "bin1")) %>%
#           dplyr::mutate(distance = abs(bin1 - bin2)) %>%
#           dplyr::select(ref1, bin1, distance, count) %>%
#           tidyr::expand(nesting(ref1, bin1), tidyr::full_seq(distance, 1)) %>%
#           dplyr::rename(distance = `tidyr::full_seq(distance, 1)`) %>%
#           dplyr::mutate(bin2 = bin1 - distance) %>%
#           dplyr::filter(bin2 >= 0)
#       distanceCount <- object@interactionMatrix %>%
#           dplyr::right_join(selectedCells, by = c("ref1", "bin1", "bin2")) %>%
#           dplyr::select(distance, count) %>%
#           tidyr::replace_na(list(count = 0)) %>%
#           dplyr::group_by(distance) %>%
#           dplyr::summarise(medCount = mean(count), .groups = "drop")
#       tmp <- selectedBins %>%
#           dplyr::mutate(ref2 = ref1) %>%
#           dplyr::mutate(bin2 = bin1)
#           tidyr::expand(tidyr::nesting(ref1, bin1), tidyr::nesting(ref2, bin2)) %>%
#           dplyr::filter(ref1 != ref2) %>%
#           dplyr::right_join(object@interactionMatrix, by = c("ref1", "bin1", "ref2", "bin2")) %>%
#           tidyr::replace_na(list(count = 0)) %>%
#           dplyr::summarise(medCount = mean(count))
#           

#           dplyr::mutate(loess = predict(loess(count ~ distance, data = ., span = 0.1))) %>%
#           dplyr::select(distance, loess) %>%
#           dplyr::distinct() %>%
#           dplyr::arrange(distance)
        
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
    message(paste0("Dataset '", object@name, "': Estimated molecule size: ", object@parameters@maxLinkRange, " (meta-bin size: ", object@parameters@metaSize, ")."))
    return(invisible(object))
}

estimateMoleculeSize <- function(object) {
    if (! is(object, "msscafClass")) {
        stop("Parameter should be a msscafClass.")
    }
    object@data <- purrr::map(object@data, .estimateMoleculeSize, sizes = object@sizes)
    return(invisible(object))
}

.estimateMetaBinsMoleculeSize <- function(object, sizes) {
    if (! is(object, "msscafExp")) {
        stop("Parameter should be a msscafExp.")
    }
    nMeta <- 10
    minCount <- 5
    moleculeSize <- (! is.null(object@parameters@maxLinkRange))
    l <- estimateMetaBinsMoleculeSizeCpp(object@interactionMatrix, sizes, minCount, nMeta, moleculeSize)
    if (l$metaSize == 0) {
        stop(paste0("Error!  Data '", object@name, "' is almost empty.  Cannot analyze it."))
    }
    object@parameters@metaSize <- l$metaSize
    if (object@parameters@metaSize > 1) {
        message(paste0("Dataset '", object@name, "': Estimated metabin size: ", object@parameters@metaSize, "."))
    }
    if (! moleculeSize) {
        object@parameters@maxLinkRange <- l$maxLinkRange
        message(paste0("Dataset '", object@name, "': Estimated molecule size: ", object@parameters@maxLinkRange, "."))
    }
    return(invisible(object))
}

estimateMetaBinsMoleculeSize <- function(object) {
    if (! is(object, "msscafClass")) {
        stop("Parameter should be a msscafClass.")
    }
    object@data <- purrr::map(object@data, .estimateMetaBinsMoleculeSize, sizes = object@sizes)
    return(invisible(object))
}


.estimateRowCount <- function(object, sizes) {
    if (! is(object, "msscafExp")) {
        stop("Parameter should be a msscafExp.")
    }
    if (object@parameters@metaSize > 1) {
        colSums <- object@interactionMatrix %>%
            computeSymmetricColSumMeta(sizes, object@parameters@metaSize) %>%
            tibble::as_tibble() %>%
            # do not use the last bin of the ref
            dplyr::mutate(size = sizes[ref]) %>%
            dplyr::filter(size - bin >= object@parameters@metaSize)
    }
    else {
        colSums <- object@interactionMatrix %>%
            computeSymmetricColSum(sizes) %>%
            tibble::as_tibble()
    }
    tmp <- colSums %>%
        dplyr::filter(sum > 0) %>%
        dplyr::pull(sum)
    quartiles <- quantile(tmp, prob = c(.25, .75))
    iqr       <- quartiles[[2]] - quartiles[[1]]
    firstOutlier <- quartiles[[1]] - 1.5 * iqr
    firstOutlier <- 0
    lastOutlier  <- quartiles[[2]] + 1.5 * iqr
    tmp <- tmp %>%
        purrr::keep(~ .x >= firstOutlier) %>%
        purrr::keep(~ .x <= lastOutlier)

    tryCatch({
            fitNB                          <- fitdistrplus::fitdist(tmp, "nbinom")
            object@parameters@rowCountSize <- fitNB$estimate[[1]]
            object@parameters@rowCountMu   <- fitNB$estimate[[2]]
            t <- tmp[pnbinom(tmp, size = fitNB$estimate[[1]], mu = fitNB$estimate[[2]]) < 0.001]
            if (length(t) > 0) {
                object@parameters@minRowCount <- max(t)
            }
            else {
                object@parameters@minRowCount <- object@parameters@minCount
            }
            t <- tmp[pnbinom(tmp, size = fitNB$estimate[[1]], mu = fitNB$estimate[[2]], lower.tail = FALSE) < 0.001]
            if (length(t) > 0) {
                object@parameters@maxRowCount <- min(t)
            }
            else {
                object@parameters@maxRowCount <- max(tmp) + 1
            }
            message(paste0("Dataset '", object@name, "': Estimated min./max. row counts: ", object@parameters@minRowCount, "/", object@parameters@maxRowCount))
        },
        error = function(e) {
            message(paste0("Dataset '", object@name, "': Cannot estimate min./max. row count, keeping default ", object@parameters@minRowCount, "/", object@parameters@maxRowCount))
            object@parameters@minRowCount <- object@parameters@minCount
            object@parameters@maxRowCount <- max(colSums$sum) + 1
        }
    )
    nRows <- nrow(colSums)
    object@outlierBins <- colSums %>%
        dplyr::filter(sum < object@parameters@minRowCount | sum > object@parameters@maxRowCount) %>%
        dplyr::select(- sum)
    message(paste0("\t", nrow(object@outlierBins), "/", nRows, " outlier bins."))
    return(invisible(object))
}

estimateRowCount <- function(object) {
    if (! is(object, "msscafClass")) {
        stop("Parameter should be a msscafClass.")
    }
    object@data <- purrr::map(object@data, .estimateRowCount, sizes = object@sizes)
    return(invisible(object))
}

.estimateMetaSizes <- function(object, sizes) {
    if (! is(object, "msscafExp")) {
        stop("Parameter should be a msscafExp.")
    }
    metaSize <- estimateMetaSizeCpp(object@interactionMatrix, object@outlierBins, sizes)
    if (metaSize == 0) {
        stop(paste0("Error!  Data '", object@name, "' is almost empty.  Cannot analyze it."))
    }
    object@parameters@metaSize <- metaSize
    message(paste0("Dataset '", object@name, "': Estimated metabin size: ", metaSize, "."))
    return(invisible(object))
}


estimateMetaSizes <- function(object) {
    if (! is(object, "msscafClass")) {
        stop("Parameter should be a msscafClass.")
    }
    object@data <- purrr::map(object@data, .estimateMetaSizes, sizes = object@sizes)
    return(invisible(object))
}

# Smoothen a (distance, count) distribution
smoothenDistribution <- function(corner) {
    # Window can be too small, use median instead
    tryCatch(
        corner %>%
            dplyr::select(distance, count) %>%
            dplyr::mutate(count = predict(loess(count ~ distance, data = .))) %>%
            dplyr::arrange(distance) %>%
            dplyr::distinct() %>%
            return(),
        error = function(e) {
            corner %>%
                dplyr::select(distance, count) %>%
                dplyr::group_by(distance) %>%
                dplyr::summarise(count = median(count), .groups = "drop") %>%
                dplyr::arrange(distance) %>%
                return()
        }
    )
}

estimateDistanceCount <- function(object, sizes) {
    if (! is(object, "msscafExp")) {
        stop("Parameter should be a msscafExp.")
    }
    message("\t\tEstimating distance/count distribution.")
    # Set genome to (ref1, bin1, distance), and transform missing values to 0
#   if (sum(sizes) <= object@parameters@sampleSize) {
#       distanceCount <- object@interactionMatrix %>%
#           dplyr::filter(ref1 == ref2) %>%
#           dplyr::mutate(distance = bin1 - bin2) %>%
#           dplyr::filter(distance <= object@parameters@maxLinkRange) %>%
#           dplyr::select(ref1, bin1, distance, count) %>%
#           tidyr::complete(nesting(ref1, bin1), distance, fill = list(count = 0)) %>%
#           dplyr::rename(ref = ref1, bin = bin1)
#   }
#   # Cannot do the loess on the whole genome (too big): sample a few positions
#   else {
#       sampleSize <- max(100, object@parameters@sampleSize / (object@parameters@maxLinkRange + 1))
#       lines <- sizes %>%
#           tibble::enframe(name = "ref", value = "size") %>%
#           dplyr::mutate(size = size - object@parameters@maxLinkRange) %>%
#           dplyr::filter(size > 0) %>%
#           dplyr::sample_n(sampleSize, replace = TRUE, weigth = size) %>%
#           dplyr::mutate(pos = runif(nrow(.))) %>%
#           dplyr::mutate(bin = as.integer(pos * size)) %>%
#           dplyr::select(ref, bin) %>%
#           dplyr::mutate(ref = factor(ref, levels = names(sizes)))
#       distanceCount <- extractLines(object@interactionMatrix, lines, object@parameters@maxLinkRange) %>%
#           tibble::as_tibble()
#   }
#   distanceCount <- object@interactionMatrix %>%
#       dplyr::filter(ref1 == ref2) %>%
#       dplyr::mutate(distance = bin1 - bin2) %>%
#       dplyr::filter(distance <= object@parameters@maxLinkRange) %>%
#       dplyr::select(ref1, bin1, distance, count) %>%
#       tidyr::complete(nesting(ref1, bin1), distance, fill = list(count = 0)) %>%
#       dplyr::select(distance, count) %>%
#       dplyr::sample_n(min(nrow(.), object@parameters@sampleSize * object@parameters@maxLinkRange)) %>%
#       dplyr::right_join(tibble(distance = seq.int(from = 0, to = object@parameters@maxLinkRange, by = 1)), by = "distance") %>%
#       tidyr::replace_na(list(count = 0))

#   distanceCount <- estimateDistanceCountCpp(object@interactionMatrix, object@outlierBins, sizes, object@parameters@maxLinkRange) %>%
#       dplyr::sample_n(min(nrow(.), object@parameters@sampleSize * object@parameters@maxLinkRange)) %>%
#       dplyr::right_join(tibble(distance = seq.int(from = 0, to = object@parameters@maxLinkRange, by = 1)), by = "distance") %>%
#       tidyr::replace_na(list(count = 0))
    distanceCount <- estimateDistanceCountCpp(object@interactionMatrix, object@outlierBins, sizes, object@parameters@maxLinkRange, object@parameters@metaSize, object@parameters@sampleSize) %>%
        tibble::as_tibble()
    object@parameters@distanceCount <- smoothenDistribution(distanceCount)
    return(object)
}

# Offset (distance, count) distribution
offsetDistribution <- function(distribution, offset) {
    distribution %>%
        dplyr::mutate(distance = distance - offset) %>%
        dplyr::filter(distance >= 0)
}

# Compute the difference between two (distance, count) distribution
cornerDifference <- function(observedDistribution, backgroundDistribution) {
#message(str(observedDistribution))
#message(str(backgroundDistribution))
#message(str(nrow(observedDistribution)))
#message(str(nrow(backgroundDistribution)))
    n <- nrow(observedDistribution)
    #if (n == 0) return(0)
    if (nrow(backgroundDistribution) != n) {
        stop(paste0("Cannot compute corner distance is the lengths differ: # observed: ", nrow(observedDistribution), ", # background: ", nrow(backgroundDistribution), "."))
    }
#   v <- observedDistribution %>%
#       dplyr::left_join(backgroundDistribution, by = "distance", suffix = c("_obs", "_back")) %>%
#       dplyr::mutate(d = count_obs - count_back) %>%
#       dplyr::pull(d) %>%
#       sum()
    v <- observedDistribution %>%
        dplyr::left_join(backgroundDistribution, by = "distance", suffix = c("_obs", "_back")) %>%
        dplyr::mutate(d = (count_obs - count_back) ^ 2) %>%
        dplyr::pull(d) %>%
        sum() %>%
        sqrt()
#message(str(v / n))
    return(v / n)
#   return(v)
}

# Apply an offset to the corner and background distribution, and compare them.
computeCornerDifferenceOffset <- function(offset, corner, background, maxDistance, bothOffsets) {
    background <- offsetDistribution(background, offset)
    if (bothOffsets) {
        corner <- offsetDistribution(corner, offset)
    }
    else {
        corner <- corner %>%
            dplyr::filter(distance <= maxDistance - offset)
    }
    cornerDifference(corner, background)
}

# Smoothen distribution, compare to background distribution with various offsets
computeCornerDifferenceOffsets <- function(corner, background, distance, bothOffsets, pb) {
    pb$tick()
    corner    <- smoothenDistribution(corner)
    offsets   <- seq.int(from = 0, to = distance - 1, by = 1)
    if (bothOffsets) {
        distances <- purrr::map_dbl(offsets, computeCornerDifferenceBothOffsetCpp, corner, background, distance)
    }
    else {
        distances <- purrr::map_dbl(offsets, computeCornerDifferenceOffsetCpp, corner, background, distance)
    }
    #distances <- purrr::map_dbl(offsets, computeCornerDifferenceOffset, corner, background, distance, bothOffsets)
    tibble::tibble(distance = offsets, score = distances)
}

extractCornerFromPoint <- function(parameters, object, pb) {
    pb$tick()
    corner <- object@interactionMatrix %>%
        dplyr::filter(ref1 == parameters$ref) %>%
        dplyr::filter(ref2 == parameters$ref) %>%
        dplyr::mutate(bin1 = bin1 - parameters$bin) %>%
        dplyr::mutate(bin2 = bin2 - parameters$bin) %>%
        dplyr::filter(bin1 >= 0, bin1 <= object@parameters@maxLinkRange) %>%
        dplyr::filter(bin2 >= 0, bin2 <= object@parameters@maxLinkRange) %>%
        dplyr::filter(bin1 >= bin2)
    fillCorner(corner, object@parameters@maxLinkRange, FALSE) %>%
        dplyr::select(distance, count)
#   corner <- object@interactionMatrix %>%
#       dplyr::filter(ref1 == parameters$ref1) %>%
#       dplyr::filter(ref2 == parameters$ref2) %>%
#       dplyr::mutate(bin1 = bin1 - parameters$bin1) %>%
#       dplyr::mutate(bin2 = bin2 - parameters$bin2) %>%
#       dplyr::filter(bin1 >= 0, bin1 <= object@parameters@maxLinkRange) %>%
#       dplyr::filter(bin2 >= 0, bin2 <= object@parameters@maxLinkRange) %>%
#       dplyr::filter(bin1 >= bin2)
#   fillCorner(corner, object@parameters@maxLinkRange, FALSE)
}

# Get random points on the diagonal
findRandomCornerPoints <- function(object, sizes, nSamples) {
    sizes %>%
        tibble::enframe(name = "ref", value = "size") %>%
        dplyr::slice_sample(n = nSamples, weight_by = size, replace = TRUE) %>%
        dplyr::mutate(size = size - object@parameters@maxLinkRange) %>%
        dplyr::mutate(bin = as.integer(round(runif(nrow(.)) * size)))
#   # Extract all observed pairs of refs
#   object@interactionMatrix %>%
#       dplyr::select(ref1, ref2) %>%
#       dplyr::filter(ref1 != ref2) %>%
#       dplyr::distinct() %>%
#   # Compute ref sizes
#       dplyr::mutate(size1 = sizes[ref1]) %>%
#       dplyr::mutate(size2 = sizes[ref2]) %>%
#       dplyr::mutate(size  = size1 * size2) %>%
#   # Sample, weighted by size
#       dplyr::slice_sample(n = nSamples, weight_by = size, replace = TRUE) %>%
#   # Get random point
#       dplyr::mutate(size1 = size1 - object@parameters@maxLinkRange) %>%
#       dplyr::mutate(size2 = size2 - object@parameters@maxLinkRange) %>%
#       dplyr::mutate(bin1 = as.integer(round(runif(nrow(.)) * size1))) %>%
#       dplyr::mutate(bin2 = as.integer(round(runif(nrow(.)) * size2))) %>%
#       dplyr::select(ref1, bin1, ref2, bin2)
#       dplyr::mutate(bin2 = as.integer(round(runif(nrow(.)) * size2))) %>%
#       dplyr::select(ref1, bin1, ref2, bin2)
}

estimateCornerVariance <- function(object, sizes, pvalueThreshold) {
    if (! is(object, "msscafExp")) {
        stop("Parameter should be a msscafExp.")
    }
    message("\t\tEstimating distance/count variance.")
    fitTriangleDifference <- function(distribution) {
        output <- tibble::tibble(shape = NA_real_, rate = NA_real_)
        try({
            f <- MASS::fitdistr(distribution$score, "gamma")
            output <- tibble::tibble(shape = f$estimate[["shape"]], rate = f$estimate[["rate"]])}, silent = TRUE)
        return(output)
    }
    nSamples <- 10 / pvalueThreshold
    pb <- progress_bar$new(total = nSamples)
    # Sample a few triangles
    object@parameters@cornerScores <- sampleTriangles(object@interactionMatrix, object@outlierBins, sizes, object@parameters@maxLinkRange, object@parameters@metaSize, nSamples) %>%
        tibble::as_tibble() %>%
        dplyr::group_by(index) %>%
        dplyr::group_split() %>%
        # Compute 
        purrr::map_dfr(computeCornerDifferenceOffsets, object@parameters@distanceCount, object@parameters@maxLinkRange, TRUE, pb) %>%
        dplyr::group_by(distance) %>%
        dplyr::group_split() %>%
        purrr::map_dfr(fitTriangleDifference, .id = "distance") %>%
        dplyr::mutate(distance = as.integer(distance)) %>%
        tidyr::drop_na()
#       # This parameter should be tuned!
#       dplyr::slice_min(score, prop = 0.5, with_ties = FALSE) %>%
#       dplyr::slice_max(score, n = 1, with_ties = FALSE) %>%
#       dplyr::ungroup()

#   points   <- findRandomCornerPoints(object, sizes, nSamples) %>%
#       purrr::transpose()
#   message("\t\tExtracting random matrices.")
#   pb       <- progress_bar$new(total = length(points))
#   corners  <- purrr::map(points, extractCornerFromPoint, object, pb) %>%
#   	purrr::map(~ dplyr::bind_cols(.x %>% dplyr::select(distance) %>% dplyr::sample_n(nrow(.)),
#   	                              .x %>% dplyr::select(count)))
#   message("\t\tComputing stats.")
#   pb       <- progress_bar$new(total = length(points))
#message("distance count")
#message(str(object@parameters@distanceCount))
#message("corner[1]")
#message(str(corners[[1]]))
#message("max link range")
#message(str(object@parameters@maxLinkRange))
#   object@parameters@cornerScores <- purrr::map_dfr(corners, computeCornerDifferenceOffsets, object@parameters@distanceCount, object@parameters@maxLinkRange, TRUE, pb) %>%
#       group_by(distance) %>%
#       # This parameter should be tuned!
#       slice_min(score, prop = 0.5, with_ties = FALSE) %>%
#       # slice_max(score, prop = pvalueThreshold, with_ties = FALSE) %>%
#       slice_max(score, n = 1, with_ties = FALSE) %>%
#       ungroup()
    return(invisible(object))
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
    if (! is(object, "msscafExp")) {
        stop("Parameter should be a msscafExp.")
    }
    message(paste0("\tDataset '", object@name, "':"))
    object   <- estimateDistanceCount(object, sizes)
    object   <- estimateCornerVariance(object, sizes, pvalueThreshold)
    object@parameters@cornerLimit <- estimateCornerLimits(object, minNBins)
    return(invisible(object))
}

estimateCorners <- function(object, pvalueThreshold) {
    if (! is(object, "msscafClass")) {
        stop("Parameter should be a msscafClass.")
    }
    message("Join shape estimations.")
    object@data <- purrr::map(object@data, .estimateCorners, object@sizes, pvalueThreshold, object@minNBins)
    return(invisible(object))
}

estimateDistributions <- function(object) {
    if (! is(object, "msscafClass")) {
        stop("Parameter should be a msscafClass.")
    }
    object <- estimateBackgroundCounts(object)
    object <- estimateMetaBinsMoleculeSize(object)
    object <- estimateRowCount(object)
    return(invisible(object))
}
