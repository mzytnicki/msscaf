# Perform test
# counts is a tibble with 2 columns:
#  - count: the counts
#  - type: the corner names, or "interior"
#computeTest <- function(counts, name) {
#    data <- counts %>%
#	dplyr::filter(type == name) %>%
#	dplyr::pull(count)
#    background <- counts %>%
#	dplyr::filter(type != name) %>%
#	dplyr::pull(count)
computePValue <- function(background, data) {
    if (length(data) == 0) {
        return(1.0)
    }
    if (max(data) <= 1) {
        return(1.0)
    }
    if ((length(background) != 0) & (length(data) == 0)) {
        return(1.0)
    }
    return(tryCatch(wilcox.test(background,
                                data,
                                paired      = FALSE,
                                exact       = FALSE,
                                alternative = "less")$p.value,
                    error = function(c) 1.0))
}

# # Check whether one corner has more counts than the rest
# # counts is a tibble with 3 columns:
# #  - chromosome1: the first reference
# #  - chromosome2: the second reference
# #  - count:       the counts
# #  - type:        the corner names, or "interior"
#     if (! is(object, "msscaf2RefExp")) {
# 	stop("Parameter should be a 'msscaf2RefExp'.")
#     }
#     counts     <- counts %>%
# 		      dplyr::filter(chromosome1 == object@chromosome1,
# 				    chromosome2 == object@chromosome2) %>%
# 		      dplyr::select(count, type)
#     testUL     <- computeTest(counts, "cornerUL")
#     testUR     <- computeTest(counts, "cornerUR")
#     testBL     <- computeTest(counts, "cornerBL")
#     testBR     <- computeTest(counts, "cornerBR")
#     testPlot <- ggplot(counts, aes(x = type, y = count)) + 
# 	geom_violin() +
# 	scale_y_log10()
#     object@interactionMatrix <- object@interactionMatrix %>%
# 				    dplyr::filter(count > 0)
#     mapPlot <- plot.10X2Ref(object, circles = TRUE)
#     progressBar$tick()
#     return(list(ref1     = object@chromosome1,
# 		ref2     = object@chromosome2,
# 		testPlot = testPlot,
# 		mapPlot  = mapPlot,
# 		tUL      = testUL,
# 		tUR      = testUR,
# 		tBL      = testBL,
# 		tBR      = testBR))
# }

# Extract the right corner, fill missing values with 0, and summarize as distance/count tibble.
extractCorner <- function(object, after1, after2) {
    interactionMatrix <- object@interactionMatrix
    if (after1) {
        interactionMatrix <- interactionMatrix %>%
            dplyr::mutate(bin1 = object@size1 - bin1)
    }
    if (after2) {
        interactionMatrix <- interactionMatrix %>%
            dplyr::mutate(bin2 = object@size2 - bin2)
    }
    fillCorner(interactionMatrix, object@parameters@minLinkRange, TRUE)
}


# Set missing values as 0
# Subsample the distribution, if needed
computeSample <- function(observed, nExpected, sampleSize = 10000) {
    # Matrix is sparse.  Fill missing values with zeros.
    # But only sample a part.
    nObserved <- length(observed)
    if (nObserved >= nExpected) {
        return(observed)
    }
    nMissing <- nExpected - nObserved
    if (nExpected < sampleSize) {
        return(c(observed, rep(0, nMissing)))
    }
    ratio <- sampleSize / nExpected
    return(c(sample(observed, nObserved * ratio), rep(0, nMissing * ratio)))
}

# Read a count tibble. Negative counts are non-corner.
# Set missing values as 0, for the corner, and the rest.
# Return a p-value.
computeTest <- function(counts, size1, size2, minDistance, maxDistance, sampleSize = 10000) {
    counts       <- counts %>% dplyr::pull(count)
    cornerCounts <- counts %>% purrr::keep(~ .x > 0)
    otherCounts  <- counts %>% purrr::keep(~ .x < 0)
    cornerSize   <- computeCornerSize(size1, size2, minDistance)
    otherSize    <- computeOtherSize(size1, size2, maxDistance)
    cornerCounts <- computeSample(cornerCounts, cornerSize)
    otherCounts  <- computeSample(- otherCounts, otherSize)
    return(computePValue(otherCounts, cornerCounts))
}

# # Extract a corner, given the corner point
# getCorner <- function(object, pointX, pointY) {
#     object@interactionMatrix %>%
# 	dplyr::filter((abs(bin1 - pointX) + abs(bin2 - pointY)) < object@parameters@minLinkRange,
# 		      abs(bin1 - pointX) < object@size1 / 2,
# 		      abs(bin2 - pointY) < object@size2 / 2) %>%
# 	dplyr::pull(count)
# }

# # Extract everything but the corners
# getInterior <- function(object) {
#     object@interactionMatrix %>%
# 	dplyr::filter((abs(bin1 - 0)            + abs(bin2 - 0))            >= object@parameters@minLinkRange,
# 		      (abs(bin1 - 0)            + abs(bin2 - object@size2)) >= object@parameters@minLinkRange,
# 		      (abs(bin1 - object@size1) + abs(bin2 - 0))            >= object@parameters@minLinkRange,
# 		      (abs(bin1 - object@size1) + abs(bin2 - object@size2)) >= object@parameters@minLinkRange) %>%
# 	dplyr::pull(count)
# }

# # Extract corners and interior of an interaction matrix
# extractRegions <- function(object) {
#     cornerSize   <- computeCornerSize(object)
#     interiorSize <- object@size1 * object@size2 - 4 * cornerSize
#     dataUL       <- getCorner(object, 0,            0           )
#     dataUR       <- getCorner(object, object@size1, 0           )
#     dataBL       <- getCorner(object, 0,            object@size2)
#     dataBR       <- getCorner(object, object@size1, object@size2)
#     interior     <- getInterior(object)
#     dataUL       <- computeSample(dataUL, cornerSize)
#     dataUR       <- computeSample(dataUR, cornerSize)
#     dataBL       <- computeSample(dataBL, cornerSize)
#     dataBR       <- computeSample(dataBR, cornerSize)
#     interior     <- computeSample(interior, interiorSize)
#     bind_rows(tibble(count = dataUL,   type = "cornerUL"),
# 	      tibble(count = dataUR,   type = "cornerUR"),
# 	      tibble(count = dataBL,   type = "cornerBL"),
# 	      tibble(count = dataBR,   type = "cornerBR"),
# 	      tibble(count = interior, type = "interior")) %>%
#       dplyr::mutate(chromosome1 = object@chromosome1) %>%
#       dplyr::mutate(chromosome2 = object@chromosome2)
# }

# # Keep region iff count is large enough
# filterRegions <- function(object, counts) {
#     tmp <- counts %>%
# 	dplyr::filter(type != "interior") %>%
# 	dplyr::summarize(maxCount = max(count), .groups = "drop") %>%
# 	dplyr::filter(maxCount > object@parameters@minCount)
#     if (nrow(tmp) > 0) {
# 	return(counts)
#     }
#     tmp <- counts %>%
# 	dplyr::filter(type != "interior") %>%
# 	dplyr::mutate(zeros = dplyr::if_else(count == 0, 1, 0)) %>%
# 	dplyr::group_by(type) %>%
# 	dplyr::summarize(nZeros    = sum(zeros),
# 			 nNonZeros = sum(1 - zeros),
# 			 .groups = "drop") %>%
# 	dplyr::filter(nNonZeros > nZeros / 2)
#     if (nrow(tmp) > 0) {
# 	return(counts)
#     }
#     dplyr::slice_head(counts, n = 0)
# }

# # Extract corners and interior of an interaction matrix
# extractAndFilterRegions <- function(object, progressBar) {
#     progressBar$tick()
#     count <- extractRegions(object)
#     filterRegions(object, count)
# }

# Extract the counts from the current join, and perform test.
testJoin <- function(counts, sizes, minDistance, maxDistance, pb) {
    firstRow <- counts[1, ]
    r1       <- firstRow$ref1
    r2       <- firstRow$ref2
    c        <- firstRow$corner
    size1    <- as.numeric(sizes[[r1]])
    size2    <- as.numeric(sizes[[r2]])
    p        <- computeTest(counts, size1, size2, minDistance, maxDistance, sampleSize = 10000)
    pb$tick()
    tibble::tibble(ref1 = r1, ref2 = r2, corner = c, pvalue = p)
}

# Add plots to experiment
getJoinInfo <- function(object, parameters) {
    testPlot <- classifyCornerPointsCpp(object@interactionMatrix, object@size1, object@size2, object@parameters@minLinkRange) %>%
        tibble::as_tibble() %>%
        computeAllSamples(object@size1, object@size2, object@parameters@minLinkRange) %>%
        tibble::enframe(name = "type", value = "count") %>%
        tidyr::unnest(count) %>%
        ggplot(aes(x = type, y = count)) + 
        geom_violin() +
        scale_y_log10()
    object@interactionMatrix <- object@interactionMatrix %>%
                                    dplyr::filter(count > 0)
    mapPlot <- plot.10X2Ref(object, circles = TRUE)
    return(list(testPlot = testPlot, mapPlot = mapPlot))
}

.checkJoins <- function(object, sizes, pvalueThreshold) {
    if (! is(object, "msscafExp")) {
        stop("Parameter should be a msscafExp.")
    }
    message(paste0("\tDataset '", object@name, "'."))
    minNBins <- 10
    selectedRefs <- sumCornerCpp(object@interactionMatrix, object@outlierBins, sizes, object@parameters@minLinkRange, object@parameters@metaSize) %>%
        tibble::as_tibble() %>%
        dplyr::filter(n > minNBins * (minNBins + 1) / 4) %>%
        dplyr::mutate(avg = count / n) %>%
        dplyr::filter(avg >= 1) %>%
        dplyr::select(ref1, ref2, corner)
#   cornerSums <- sumCornerCpp(object@interactionMatrix, object@outlierBins, sizes, object@parameters@minLinkRange, object@parameters@metaSize) %>%
#       tibble::as_tibble() %>%
#       # Distribution is inflated in 0, remove it
#       dplyr::filter(count > 0) %>%
#       # Remove strange even/odd pattern
#       dplyr::mutate(count = as.numeric(trunc((count - 1)/2)))
#   if (nrow(cornerSums) == 0) {
#       message("\t\tNo join found.")
#       joinsObject <- new("msscafJoins")
#       joinsObject@data <- tibble::tibble(ref1 = integer(), ref2 = integer(), after1 = logical(), after2 = logical(), pvalue = numeric(), pvalueCorner = numeric())
#       object@joins <- joinsObject
#       return(object)
#   }
#   # Select outlier corner sums
#   # Not meaningful if you have too few
#   if (nrow(cornerSums) >= 30) {
#       fitNB <- fitdistrplus::fitdist(cornerSums$count, "nbinom")
#       threshold <- qnbinom(0.99, size = fitNB$estimate[["size"]], mu = fitNB$estimate[["mu"]])
#       selectedRefs <- cornerSums %>%
#           dplyr::filter(count > threshold)
#   }
#   else {
#       selectedRefs <- cornerSums
#   }
    message(paste0("\t\t", nrow(selectedRefs), " selected joins."))
    pb <- progress_bar$new(total = nrow(selectedRefs))
    selectedRefs <- extractCornersCpp(object@interactionMatrix, selectedRefs, object@outlierBins, sizes, object@parameters@minLinkRange, object@parameters@maxLinkRange, object@parameters@metaSize) %>%
        tibble::as_tibble() %>%
        dplyr::group_by(ref1, ref2, corner) %>%
        dplyr::group_split() %>%
        purrr::map_dfr(testJoin, sizes = sizes, minDistance = object@parameters@minLinkRange, maxDistance = object@parameters@maxLinkRange, pb = pb) %>%
#       dplyr::mutate(size1 = sizes[ref1]) %>%
#       dplyr::mutate(size2 = sizes[ref2]) %>%
#       dplyr::mutate(minSize = pmin(size1, size2)) %>%
#       # Do not apply p-value filter if refs are too short: corners overlap
#       dplyr::filter((pvalue <= pvalueThreshold) | (minSize <= 4 * object@parameters@minLinkRange)) %>%
        dplyr::filter(pvalue <= pvalueThreshold) %>%
        dplyr::group_by(ref1, ref2) %>% # Possible 2 joins for the same pairs of refs?
        dplyr::slice_min(pvalue, n = 1, with_ties = FALSE) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(corner = as.character(corner)) %>%
        dplyr::mutate(after1 = (stringr::str_sub(corner, 1, 1) == "E")) %>%
        dplyr::mutate(after2 = (stringr::str_sub(corner, 2, 2) == "E")) %>%
        dplyr::select(ref1, ref2, after1, after2, pvalue)
    message(paste0("\t\tKeeping ", nrow(selectedRefs), " of them."))
#   selectedCounts <- keepScaffoldsPairsCpp(object@interactionMatrix, selectedRefs) %>% tibble::as_tibble()
#   objects        <- splitBy2RefFromMatrix(object, selectedCounts, sizes)
#   joinInfo       <- purrr::map2(objects, purrr::transpose(selectedRefs), getJoinInfo) %>% purrr::transpose()
    joinsObject <- new("msscafJoins")
    joinsObject@data <- selectedRefs
#   joinsObject@testPlots <- joinInfo$testPlot
#   joinsObject@mapPlots <- joinInfo$mapPlot
    object@joins <- joinsObject
    return(object)

#   pb <- progress_bar$new(total = length(objects))
#   counts <- purrr::map_dfr(objects, extractAndFilterRegions, progressBar = pb)
#   keptRefs <- counts %>%
#       dplyr::select(chromosome1, chromosome2) %>%
#       tidyr::unite("refs", c("chromosome1", "chromosome2"), sep = "_") %>%
#       dplyr::distinct() %>%
#       dplyr::pull(refs)
#   if (! all(keptRefs %in% names(objects))) stop("Problem while selecting references.") 
#   objects <- objects[keptRefs]
#   pb <- progress_bar$new(total = length(objects))
#   joins <- purrr::map(objects, checkJoin, counts, progressBar = pb)
#   if (length(joins) == 0) {
#       message("\tNo join found.")
#       joinsObject <- new("msscafJoins")
#       joinsObject@data <- tibble(ref1   = character(),
#                                  ref2   = character(),
#                                  vert   = character(),
#                                  hor    = character(),
#                                  pvalue = numeric())
#       joinsObject@testPlots <- c()
#       joinsObject@mapPlots <- c()
#       object@joins <- joinsObject
#       return(object)
#   }
#   joins <- tibble(ref1     = purrr::map_chr(joins, "ref1"),
#                   ref2     = purrr::map_chr(joins, "ref2"),
#                   tUL      = purrr::map(joins, "tUL"),
#                   tUR      = purrr::map(joins, "tUR"),
#                   tBL      = purrr::map(joins, "tBL"),
#                   tBR      = purrr::map(joins, "tBR"),
#                   testPlot = purrr::map(joins, "testPlot"),
#                   mapPlot  = purrr::map(joins, "mapPlot")) %>%
#       tidyr::gather(key = "key", value = "test", tUL, tUR, tBL, tBR) %>%
#       dplyr::filter(!unlist(purrr::map(test, rlang::is_null))) %>%
#       dplyr::mutate(pvalue = unlist(test)) %>%
#       tidyr::separate(key, c(1, 2), into = c(NA, "vert", "hor")) %>%
#       dplyr::mutate(vert = factor(vert)) %>%
#       dplyr::mutate(hor = factor(hor)) %>%
#       dplyr::filter(pvalue <= pvalueThreshold) %>%
#       dplyr::arrange(pvalue)
#   joinsObject <- new("msscafJoins")
#   joinsObject@data <- joins %>%
#       dplyr::select(-c(test, testPlot, mapPlot))
#   joinsObject@testPlots <- joins %>%
#       dplyr::select(ref1, ref2, testPlot) %>%
#       tidyr::unite("name", ref1, ref2) %>%
#       deframe()
#   joinsObject@mapPlots <- joins %>%
#       dplyr::select(ref1, ref2, mapPlot) %>%
#       tidyr::unite("name", ref1, ref2) %>%
#       deframe()
#   object@joins <- joinsObject
#   return(object)
}

checkJoins <- function(object, pvalueThreshold) {
    if (! is(object, "msscafClass")) {
        stop("Parameter should be a msscafClass.")
    }
    message("Finding putative joins.")
    object@data <- purrr::map(object@data, .checkJoins, sizes = object@sizes, pvalueThreshold)
    return(invisible(object))
}

transformRefAfterToRefRef <- function(joins) {
    dplyr::bind_rows(joins %>%
            dplyr::select(ref1, after1, pvalue) %>%
            dplyr::rename(ref = ref1, after = after1),
        joins %>%
            dplyr::select(ref2, after2, pvalue) %>%
            dplyr::rename(ref = ref2, after = after2))
}

discardJoinsFromRefRef <- function(joins, toDiscard) {
    # Transform to ref1/bin1
    refBin1 <- toDiscard %>%
        dplyr::rename(ref1 = ref) %>%
        dplyr::rename(after1 = after) %>%
        dplyr::select(ref1, after1)
    # Transform to ref2/bin2
    refBin2 <- toDiscard %>%
        dplyr::rename(ref2 = ref) %>%
        dplyr::rename(after2 = after) %>%
        dplyr::select(ref2, after2)
    joins %>%
        dplyr::anti_join(refBin1, by = c("ref1", "after1")) %>%
        dplyr::anti_join(refBin2, by = c("ref2", "after2"))
}

.removeAmbiguousJoins <- function(joins, name = NULL) {
    if (! is.null(name)) {
        message(paste0("\tData '", name, "'."))
    }
    # Transform joins to before/after
    nJoins <- nrow(joins)
    if (nJoins == 0) {
        return(joins)
    }
    ambiguousJoins <- joins %>%
        transformRefAfterToRefRef() %>%
    # Find duplicates
        dplyr::group_by(ref, after) %>%
        dplyr::slice_min(order_by = pvalue, n = 2, with_ties = FALSE) %>%
        dplyr::mutate(logPvalue = -log10(pvalue)) %>%
        dplyr::summarize(pvalueDiff = max(logPvalue) - min(logPvalue), n = dplyr::n(), .groups = "drop") %>%
    # Select duplicates
        dplyr::filter(n > 1, pvalueDiff <= 2) %>%
        dplyr::select(ref, after)
    joins <- discardJoinsFromRefRef(joins, ambiguousJoins)
    if (! is.null(name)) {
        message(paste0("\t\tKept ", nrow(joins), "/", nJoins, "."))
    }
    return(joins)
}

removeAmbiguousJoins <- function(object) {
    if (! is(object, "msscafClass")) {
        stop("Parameter should be a msscafClass.")
    }
    message("Removing ambiguous joins")
    f <- function(d) {
        d@joins@data <- .removeAmbiguousJoins(d@joins@data, d@name)
        return(d)
    }
    object@data <- purrr::map(object@data, f)
    return(invisible(object))
}

removeSuboptimalJoins <- function(joins) {
    nJoins <- nrow(joins)
    if (nJoins == 0) {
        return(joins)
    }
    maxDistance <- joins %>%
        dplyr::pull(distance) %>%
        max()
    minPvalue <- joins %>%
        dplyr::filter(pvalue > 0) %>%
        dplyr::pull(pvalue) %>%
        min()
    minPvalueCorner <- joins %>%
        dplyr::filter(pvalueCorner > 0) %>%
        dplyr::pull(pvalueCorner) %>%
        min()
    # Try to favor close interactions, and good p-values
    joins <- joins %>%
        dplyr::mutate(score = -log(pvalue + minPvalue / 2) - log(pvalueCorner + minPvalueCorner / 2) * (maxDistance - distance)) %>%
        dplyr::arrange(desc(score))
    selected <- joins %>% dplyr::slice_head(n = 0)
    # Iteratively select best join, and discard possible conflicting joins
    nJoins <- nrow(joins)
    message("\tRemoving sub-optimal joins.")
    pb <- progress_bar$new(total = nJoins)
    while (nJoins > 0) {
        best  <- joins %>% dplyr::slice(1)
        selected <- selected %>%
            dplyr::bind_rows(best)
        bestList <- best %>% purrr::transpose() %>% purrr::pluck(1)
        joins <- joins %>%
            dplyr::slice(-1) %>%
            dplyr::filter((ref1 != bestList$ref1) | (after1 != bestList$after1)) %>%
            dplyr::filter((ref2 != bestList$ref1) | (after2 != bestList$after1)) %>%
            dplyr::filter((ref1 != bestList$ref2) | (after1 != bestList$after2)) %>%
            dplyr::filter((ref2 != bestList$ref2) | (after2 != bestList$after2))
        pb$tick(nJoins - nrow(joins))
        nJoins <- nrow(joins)
    }
    return(selected)
}

..checkCornersOld <- function(parameters, object, sizes, pb) {
    objectRef <- extract2Ref(object, parameters$ref1, parameters$ref2, sizes[[parameters$ref1]], sizes[[parameters$ref2]])
    corner    <- extractCorner(objectRef, parameters$after1, parameters$after2)
    background <- object@parameters@distanceCount
    values    <- computeCornerDifferenceOffsets(corner, background, object@parameters@minLinkRange, FALSE, pb) %>%
        dplyr::left_join(object@parameters@cornerScores, by = "distance", suffix = c("_corner", "_background")) %>%
        # Score of -1 means that nothing matched (because of outlier bins)
        dplyr::filter(score_corner >= 0) %>%
        dplyr::filter(score_corner <= score_background)
    if (nrow(values) == 0) return(-1)
    values %>%
        dplyr::slice_min(distance, n = 1, with_ties = FALSE) %>%
        dplyr::pull(distance) %>%
        return()
}

.checkCornersOld <- function(object, sizes, pvalueThreshold) {
    message(paste0("\tDataset '", object@name, "'."))
    nJoins <- nrow(object@joins@data)
    pb     <- progress_bar$new(total = nrow(object@joins@data))
    values <- object@joins@data %>%
        dplyr::mutate(ref1 = as.character(ref1)) %>%
        dplyr::mutate(ref2 = as.character(ref2)) %>%
        purrr::transpose() %>%
        purrr::map_dbl(..checkCorners, object, sizes, pb)
    object@joins@data <- object@joins@data %>%
        dplyr::mutate(offset = values) %>%
        dplyr::filter(offset >= 0)
    message(paste0("\t\tKept ", nrow(object@joins@data), "/", nJoins, "."))
    return(object)
}

..checkCorners <- function(corner, object, sizes, pb) {
    maxDistance    <- object@parameters@minLinkRange - 8
    values     <- computeCornerDifferenceOffsets(corner, object@parameters@distanceCount, object@parameters@minLinkRange, FALSE, pb) %>%
        dplyr::filter(distance > 0) %>%
        dplyr::filter(distance <= maxDistance) %>%
        dplyr::left_join(object@parameters@cornerScores, by = "distance") %>%
#       dplyr::slice_min(score, n = 1, with_ties = FALSE) %>%
#       dplyr::mutate(pvalueCorner = pgamma(score, shape = shape, rate = rate)) %>%
        dplyr::mutate(pvalueCorner = pgamma(score, shape = shape, rate = rate)) %>%
        dplyr::slice_min(pvalueCorner, n = 1, with_ties = FALSE) %>%
        dplyr::select(distance, pvalueCorner) %>%
        dplyr::mutate(distance = as.integer(distance) * object@parameters@metaSize)
#   background <- object@parameters@distanceCount
#   values     <- computeCornerDifferenceOffsets(corner, background, object@parameters@minLinkRange, FALSE, pb) %>%
#       dplyr::left_join(object@parameters@cornerScores, by = "distance", suffix = c("_corner", "_background")) %>%
#       dplyr::filter(score_corner <= score_background)
#   if (nrow(values) == 0) return(-1)
#   values %>%
#       dplyr::slice_min(distance, n = 1, with_ties = FALSE) %>%
#       dplyr::pull(distance) %>%
#       return()
}

.checkCorners <- function(object, sizes, pvalueThreshold) {
    corners <- extractCornersFullCpp(object@interactionMatrix, object@joins@data, object@outlierBins, sizes, object@parameters@minLinkRange, object@parameters@metaSize) %>%
        tibble::as_tibble() %>%
        # Unused cells are set to -1
        dplyr::filter(count >= 0)
    nJoins <- nrow(object@joins@data)
    message(paste0("\tSelected ", nJoins, " joins."))
    if (nJoins == 0) {
        return(object)
    }
    pb     <- progress_bar$new(total = nJoins)
    values <- corners %>%
        dplyr::group_by(index) %>%
        dplyr::group_split() %>%
        purrr::map_dfr(..checkCorners, object, sizes, pb)
    object@joins@data <- object@joins@data %>%
        dplyr::bind_cols(values) %>%
        dplyr::filter(pvalueCorner < 1 - pvalueThreshold)
    message(paste0("\t\tKept ", nrow(object@joins@data), "."))
    return(object)
}

checkCorners <- function(object, pvalueThreshold) {
    if (! is(object, "msscafClass")) {
        stop("Parameter should be a msscafClass.")
    }
    message("Checking corner shape.")
    object@data <- purrr::map(object@data, .checkCorners, object@sizes, pvalueThreshold)
    return(invisible(object))
}

mergeJoins <- function(object, pvalueThreshold) {
    if (! is(object, "msscafClass")) {
        stop("Parameter should be a msscafClass.")
    }
    object@joins <- dplyr::bind_rows(purrr::map(purrr::map(object@data, "joins"), "data")) %>%
        dplyr::distinct() %>%
        dplyr::arrange(pvalue) %>%
        dplyr::mutate(pvalue = pmin(pvalue, pvalueCorner)) %>%
        dplyr::mutate(p_adj = p.adjust(pvalue, method = "BH")) %>%
        dplyr::filter(p_adj <= pvalueThreshold) %>%
        dplyr::group_by(ref1, ref2, after1, after2) %>%
        dplyr::arrange(pvalue) %>%
        dplyr::slice_head(n = 1) %>%
        dplyr::ungroup() %>%
        .removeAmbiguousJoins() %>%
        removeSuboptimalJoins()
    message(paste0("\t", nrow(object@joins), " joins found in total."))
    return(invisible(object))                                                                                                                                                                                      
}

..addJoinPlots <- function(object, ref1, ref2, size1, size2, left, up) {
    x <- ifelse(left, 0, size1)
    y <- ifelse(up, 0, size2)
    object2Ref <- extract2Ref(object, ref1, ref2, size1, size2)
    return(plot.10X2Ref(object2Ref, logColor = TRUE, circles = TRUE, center = list(x, y), radius = 100))
}

.addJoinPlots <- function(parameters, object, pb) {
    ref1             <- parameters$ref1
    ref2             <- parameters$ref2
    joinPlots        <- purrr::map(object@data, ..addJoinPlots, ref1 = ref1, ref2 = ref2, size1 = object@sizes[[ref1]], size2 = object@sizes[[ref2]], left = parameters$left, up = parameters$up)
    names(joinPlots) <- purrr::map(object@data, "name")
    pb$tick()
    return(joinPlots)
}

addJoinPlots <- function(object) {
    if (! is(object, "msscafClass")) {
	stop("Parameter should be a msscafClass.")
    }
    message(paste0("\tComputing plots."))
    pb                      <- progress_bar$new(total = nrow(object@joins))
    object@joinPlots        <- purrr::map(object@joins %>%
				       dplyr::select(ref1, ref2, after1, after2) %>%
				       purrr::transpose(),
				   .addJoinPlots, object = object, pb = pb)
    names(object@joinPlots) <- object@joins %>%
        dplyr::mutate(after1 = dplyr::if_else(after1, "E", "B")) %>%
        dplyr::mutate(after2 = dplyr::if_else(after2, "E", "B")) %>%
        tidyr::unite("name", c(after1, after2), sep = "") %>%
        tidyr::unite("name", c(ref1, ref2, name), sep = "_") %>%
        dplyr::pull(name)
    return(object)                                                                                                                                                                                      
}

findJoins <- function(object, pvalue = 0.05) {
    if (! is(object, "msscafClass")) {
        stop("Parameter should be a msscafClass.")
    }
    object <- estimateCorners(object, pvalue)
    object <- checkJoins(object, pvalue)
    object <- removeAmbiguousJoins(object)
    object <- checkCorners(object, pvalue)
    object <- mergeJoins(object, pvalue)
    #object <- addJoinPlots(object)
    gc(verbose = FALSE)
    return(invisible(object))                                                                                                                                                                                      
}
