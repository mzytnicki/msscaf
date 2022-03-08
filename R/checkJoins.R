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
computeTest <- function(background, data) {
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
# checkJoin <- function(object, counts, progressBar) {
#     if (! is(object, "tenxchecker2RefExp")) {
# 	stop("Parameter should be a 'tenxchecker2RefExp'.")
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
    fillCorner(interactionMatrix, object@parameters@maxLinkRange, TRUE)
}


# Set missing values as 0
# Subsample the distribution, if needed
computeSample <- function(observed, nExpected, sampleSize = 10000) {
    # Matrix is sparse.  Fill missing values with zeros.
    # But only sample a part.
    nObserved <- length(observed)
    if (nObserved >= nExpected) {
        return(nObserved)
    }
    nMissing <- nExpected - nObserved
    if (nExpected < sampleSize) {
        return(c(observed, rep(0, nMissing)))
    }
    ratio <- sampleSize / nExpected
    return(c(sample(observed, nObserved * ratio), rep(0, nMissing * ratio)))
}

# Read a (type, count) tibble.
# Set missing values as 0, for all the types.
# Return a list of distributions.
computeAllSamples <- function(counts, size1, size2, maxDistance, sampleSize = 10000) {
    cornerTypes   <- counts %>% dplyr::pull(type) %>% levels()
    cornerSize    <- computeCornerSize(size1, size2, maxDistance)
    interiorSize  <- size1 * size2 - 4 * cornerSize
    splitCounts   <- counts %>%
        dplyr::select(type, count)
    splitCounts <- split(splitCounts$count, splitCounts$type)
    expectedSizes <- c(rep.int(cornerSize, 4), interiorSize)
    samples <- purrr::map2(splitCounts, expectedSizes, computeSample)
    names(samples) <- cornerTypes
    return(samples)
}

# # Extract a corner, given the corner point
# getCorner <- function(object, pointX, pointY) {
#     object@interactionMatrix %>%
# 	dplyr::filter((abs(bin1 - pointX) + abs(bin2 - pointY)) < object@parameters@maxLinkRange,
# 		      abs(bin1 - pointX) < object@size1 / 2,
# 		      abs(bin2 - pointY) < object@size2 / 2) %>%
# 	dplyr::pull(count)
# }

# # Extract everything but the corners
# getInterior <- function(object) {
#     object@interactionMatrix %>%
# 	dplyr::filter((abs(bin1 - 0)            + abs(bin2 - 0))            >= object@parameters@maxLinkRange,
# 		      (abs(bin1 - 0)            + abs(bin2 - object@size2)) >= object@parameters@maxLinkRange,
# 		      (abs(bin1 - object@size1) + abs(bin2 - 0))            >= object@parameters@maxLinkRange,
# 		      (abs(bin1 - object@size1) + abs(bin2 - object@size2)) >= object@parameters@maxLinkRange) %>%
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
# 	dplyr::mutate(zeros = if_else(count == 0, 1, 0)) %>%
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
testJoin <- function(parameters, counts, sizes, maxDistance, pb) {
    counts <- counts %>%
        dplyr::mutate(ref1 = as.numeric(ref1)) %>%
        dplyr::mutate(ref2 = as.numeric(ref2)) %>%
        dplyr::filter(ref1 == parameters$ref1) %>%
        dplyr::filter(ref2 == parameters$ref2)
    size1         <- sizes[parameters$ref1]
    size2         <- sizes[parameters$ref2]
    samples       <- computeAllSamples(counts, size1, size2, maxDistance, sampleSize = 10000)
    testedCorner  <- samples[[parameters$corner]]
    samples[[parameters$corner]]  <- NULL
    pb$tick()
    min(purrr::map_dbl(samples, computeTest, testedCorner))
}

# Add plots to experiment
getJoinInfo <- function(object, parameters) {
    testPlot <- classifyCornerPointsCpp(object@interactionMatrix, object@size1, object@size2, object@parameters@maxLinkRange) %>%
        as_tibble() %>%
        computeAllSamples(object@size1, object@size2, object@parameters@maxLinkRange) %>%
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
    if (! is(object, "tenxcheckerExp")) {
        stop("Parameter should be a tenxcheckerExp.")
    }
    message(paste0("\tDataset '", object@name, "'."))
    selectedRefs   <- filterCornersCpp(object@interactionMatrix, sizes, object@parameters@maxLinkRange, object@parameters@metaSize) %>% as_tibble()
    message(paste0("\t\t", nrow(selectedRefs), " selected joins."))
    pb <- progress_bar$new(total = nrow(selectedRefs))
    selectedCounts <- extractCornersCpp(object@interactionMatrix, selectedRefs, sizes, object@parameters@maxLinkRange, object@parameters@metaSize) %>% as_tibble()
    minPValues     <- purrr::map_dbl(purrr::transpose(selectedRefs), testJoin, counts = selectedCounts, sizes = sizes, maxDistance = object@parameters@maxLinkRange, pb = pb)
    selectedRefs   <- selectedRefs %>%
        dplyr::mutate(pvalue = minPValues) %>%
        dplyr::filter(pvalue <= pvalueThreshold) %>%
        dplyr::group_by(ref1, ref2) %>% # Possible 2 joins for the same pairs of refs?
        dplyr::slice_min(pvalue, n = 1, with_ties = FALSE) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(corner = as.character(corner)) %>%
        dplyr::mutate(after1 = (stringr::str_sub(corner, 1, 1) == "E")) %>%
        dplyr::mutate(after2 = (stringr::str_sub(corner, 2, 2) == "E")) %>%
        dplyr::select(ref1, ref2, after1, after2, pvalue)
    message(paste0("\t\tKeeping ", nrow(selectedRefs), " of them."))
#   selectedCounts <- keepScaffoldsPairsCpp(object@interactionMatrix, selectedRefs) %>% as_tibble()
#   objects        <- splitBy2RefFromMatrix(object, selectedCounts, sizes)
#   joinInfo       <- purrr::map2(objects, purrr::transpose(selectedRefs), getJoinInfo) %>% purrr::transpose()
    joinsObject <- new("tenxcheckerJoins")
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
#       joinsObject <- new("tenxcheckerJoins")
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
#   joins <- tibble(ref1     = map_chr(joins, "ref1"),
#                   ref2     = map_chr(joins, "ref2"),
#                   tUL      = map(joins, "tUL"),
#                   tUR      = map(joins, "tUR"),
#                   tBL      = map(joins, "tBL"),
#                   tBR      = map(joins, "tBR"),
#                   testPlot = map(joins, "testPlot"),
#                   mapPlot  = map(joins, "mapPlot")) %>%
#       tidyr::gather(key = "key", value = "test", tUL, tUR, tBL, tBR) %>%
#       dplyr::filter(!unlist(map(test, is_null))) %>%
#       dplyr::mutate(pvalue = unlist(test)) %>%
#       tidyr::separate(key, c(1, 2), into = c(NA, "vert", "hor")) %>%
#       dplyr::mutate(vert = factor(vert)) %>%
#       dplyr::mutate(hor = factor(hor)) %>%
#       dplyr::filter(pvalue <= pvalueThreshold) %>%
#       dplyr::arrange(pvalue)
#   joinsObject <- new("tenxcheckerJoins")
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
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    message("Finding putative joins.")
    object@data <- map(object@data, .checkJoins, sizes = object@sizes, pvalueThreshold)
    return(invisible(object))
}

.removeDuplicateJoins <- function(object) {
    if (! is(object, "tenxcheckerExp")) {
        stop("Parameter should be a tenxcheckerExp.")
    }
    message(paste0("\tDataset '", object@name, "'."))
    # Transform joins to before/after
    nJoins <- nrow(object@joins@data)
#   tmp <- object@joins@data
#       dplyr::mutate(hor = as.character(hor)) %>%
#       dplyr::mutate(vert = as.character(vert)) %>%
    ambiguousJoins <- dplyr::bind_rows(object@joins@data %>%
                                           dplyr::select(ref1, after1, pvalue) %>%
                                           dplyr::rename(ref = ref1, after = after1),
                                       object@joins@data %>%
                                           dplyr::select(ref2, after2, pvalue) %>%
                                           dplyr::rename(ref = ref2, after = after2)) %>%
    # Find duplicates
        dplyr::group_by(ref, after) %>%
        dplyr::slice_min(order_by = pvalue, n = 2, with_ties = FALSE) %>%
        dplyr::mutate(logPvalue = -log10(pvalue)) %>%
        dplyr::summarize(pvalueDiff = max(logPvalue) - min(logPvalue), n = n(), .groups = "drop") %>%
    # Select duplicates
        dplyr::filter(n > 1, pvalueDiff <= 2) %>%
        dplyr::select(ref, after)
    # Transform to ref1/bin1
    refBin1 <- ambiguousJoins %>%
        dplyr::rename(ref1 = ref) %>%
        dplyr::rename(after1 = after) %>%
        dplyr::select(ref1, after1)
    # Transform to ref2/bin2
    refBin2 <- ambiguousJoins %>%
        dplyr::rename(ref2 = ref) %>%
        dplyr::rename(after2 = after) %>%
        dplyr::select(ref2, after2)
    object@joins@data <- object@joins@data %>%
        dplyr::anti_join(refBin1, by = c("ref1", "after1")) %>%
        dplyr::anti_join(refBin2, by = c("ref2", "after2"))
    message(paste0("\t\tKept ", nrow(object@joins@data), "/", nJoins, "."))
    return(object)
}

removeDuplicateJoins <- function(object) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    message("Removing ambiguous joins")
    object@data <- map(object@data, .removeDuplicateJoins)
    return(invisible(object))
}

..checkCorners <- function(parameters, object, sizes, pb) {
    objectRef <- extract2Ref(object, parameters$ref1, parameters$ref2, sizes[[parameters$ref1]], sizes[[parameters$ref2]])
    corner    <- extractCorner(objectRef, parameters$after1, parameters$after2)
    background <- object@parameters@distanceCount
    values    <- computeCornerDifferenceOffsets(corner, background, object@parameters@maxLinkRange, FALSE, pb) %>%
        dplyr::left_join(object@parameters@cornerScores, by = "distance", suffix = c("_corner", "_background")) %>%
        dplyr::filter(score_corner <= score_background)
    if (nrow(values) == 0) return(-1)
    values %>%
        dplyr::slice_min(distance, n = 1, with_ties = FALSE) %>%
        dplyr::pull(distance) %>%
        return()
}

.checkCorners <- function(object, sizes, pvalueThreshold) {
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

checkCorners <- function(object, pvalueThreshold) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    message("Checking corner shape.")
    object@data <- map(object@data, .checkCorners, object@sizes, pvalueThreshold)
    return(invisible(object))
}

mergeJoins <- function(object, pvalueThreshold) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    object@joins <- dplyr::bind_rows(map(map(object@data, "joins"), "data")) %>%
        dplyr::distinct() %>%
        dplyr::arrange(pvalue) %>%
        dplyr::mutate(p_adj = p.adjust(pvalue, method = "BH")) %>%
        dplyr::filter(p_adj <= pvalueThreshold) %>%
        dplyr::group_by(ref1, ref2, after1, after2) %>%
        dplyr::arrange(pvalue) %>%
        dplyr::slice_head(n = 1) %>%
        dplyr::ungroup()
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
    joinPlots        <- map(object@data, ..addJoinPlots, ref1 = ref1, ref2 = ref2, size1 = object@sizes[[ref1]], size2 = object@sizes[[ref2]], left = parameters$left, up = parameters$up)
    names(joinPlots) <- map(object@data, "name")
    pb$tick()
    return(joinPlots)
}

addJoinPlots <- function(object) {
    if (! is(object, "tenxcheckerClass")) {
	stop("Parameter should be a tenxcheckerClass.")
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
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    object <- estimateCorners(object, pvalue)
    object <- checkJoins(object, pvalue)
    object <- removeDuplicateJoins(object)
    object <- checkCorners(object, pvalue)
    object <- mergeJoins(object, pvalue)
    #object <- addJoinPlots(object)
    gc(verbose = FALSE)
    return(invisible(object))                                                                                                                                                                                      
}
