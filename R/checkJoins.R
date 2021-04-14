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

computeTest <- function(background, nBackground, data, nData) {
#   message(paste0("# back seen: ", length(background)))
#   message(paste0("# back exp: ", nBackground))
#   message(paste0("# data: ", length(data)))
#   message(paste0("# data > 0: ", length(data[data != 0])))
#   message(paste0("max back: ", max(background)))
#   message(paste0("max data: ", max(data)))
#   message(paste0("med back: ", median(background)))
#   message(paste0("med data: ", median(data)))
    if (length(data) == 0) {
        return(1.0)
    }
    if (max(data) <= 1) {
        return(1.0)
    }
#   message("check 1")
    if ((length(background) != 0) & (length(data) == 0)) {
        return(1.0)
    }
#   message("check 2")
#   if ((length(background) == 0) & (length(data) != 0)) {
#       if (length(data) >= 10) return(0.0)
#       return(1.0)
#   }
#   message("check 3")
    background <- computeSample(background, nBackground)
    data       <- computeSample(data,       nData)
#   if (length(data[data != 0]) < (length(data) / 2)) {
#       return(1.0)
#   }
#   message("check 4")
#   if ((max(background) < 3) | (max(data) < 3)) {
#       return(1.0)
#   }
#   if (median(background) > (median(data) / 2)) {
#       return(1.0)
#   }
#   message("back: ")
#   message(str(table(background)))
#   message("data: ")
#   message(str(table(data)))
#   message(str(wilcox.test(background,                                                                                                                                                                        
#                               data,                                                                                                                                                                              
#                               paired      = FALSE,                                                                                                                                                               
#                               exact       = FALSE,                                                                                                                                                               
#                               alternative = "less")))
#   message("check 6")
    return(tryCatch(wilcox.test(background,
                                data,
                                paired      = FALSE,
                                exact       = FALSE,
                                alternative = "less")$p.value,
                    error = function(c) 1.0))
}

compareRegions <- function(object, rest, pointX, pointY, cornerSize, name) {
    # message(str(name))
    f <- function (x, y) {
        #return(((sqrt((x - pointX)^2 + (y - pointY)^2)) < object@parameters@maxLinkRange) &
        return(((abs(x - pointX) + abs(y - pointY)) < object@parameters@maxLinkRange) &
               (abs(x - pointX) < object@size1 / 2) &
               (abs(y - pointY) < object@size2 / 2))
    }
    corner <- object@interactionMatrix %>%
        filter(f(bin1, bin2)) %>%
        select(count) %>%
        mutate(type = name)
    interior <- object@interactionMatrix %>%
        filter(! f(bin1, bin2))
    rest <- rest %>%
        filter(! f(bin1, bin2))
#   message(name)
#   message(paste0("size 1: ", object@size1))
#   message(paste0("size 2: ", object@size2))
    test <- computeTest(interior$count, as.numeric(object@size1) * as.numeric(object@size2) - cornerSize, corner$count, cornerSize)
#   message(str(test))
    return(list(data = corner, test = test, rest = rest))
}

computeCornerSize <- function(object) {
    cornerSize <- object@parameters@maxLinkRange * (object@parameters@maxLinkRange + 1) / 2
    size1      <- object@parameters@maxLinkRange - (object@size1 / 2)
    size2      <- object@parameters@maxLinkRange - (object@size2 / 2)
    if (size1 > 0) {
        cornerSize <- cornerSize - (size1 * (size1 + 1) / 2)
    }
    if (size2 > 0) {
        cornerSize <- cornerSize - (size2 * (size2 + 1) / 2)
    }
    return(cornerSize)
}

checkJoin <- function(object, progressBar) {
    if (! is(object, "tenxchecker2RefExp")) {
        stop("Parameter should be a 'tenxchecker2RefExp'.")
    }
#   message(paste0("chr1: ", object@chromosome1))
#   message(paste0("chr2: ", object@chromosome2))
    cornerSize <- computeCornerSize(object)
    rest       <- object@interactionMatrix
#   object     <- tryCatch(normalizeNotSquareICE(object),
#       error = function(cond) { 
#           object@interactionMatrix <- makeSymmetric(object@interactionMatrix)
#           return(object)
#       })
    dataUL     <- compareRegions(object, rest,        1,            1,            cornerSize, "cornerUL")
    dataUR     <- compareRegions(object, dataUL$rest, object@size1, 1,            cornerSize, "cornerUR")
    dataBL     <- compareRegions(object, dataUR$rest, 1,            object@size2, cornerSize, "cornerBL")
    dataBR     <- compareRegions(object, dataBL$rest, object@size1, object@size2, cornerSize, "cornerBR")
    rest       <- dataBR$rest %>%
        select(count) %>%
        mutate(type = "interior")
    dataUL$test <- max(dataUL$test,
                       computeTest(dataUR$data$count, cornerSize, dataUL$data$count, cornerSize),
                       computeTest(dataBL$data$count, cornerSize, dataUL$data$count, cornerSize),
                       computeTest(dataBR$data$count, cornerSize, dataUL$data$count, cornerSize),
                       na.rm = TRUE)
#   message(paste0("UR1: ", dataUR$test))
    dataUR$test <- max(dataUR$test,
                       computeTest(dataUL$data$count, cornerSize, dataUR$data$count, cornerSize),
                       computeTest(dataBL$data$count, cornerSize, dataUR$data$count, cornerSize),
                       computeTest(dataBR$data$count, cornerSize, dataUR$data$count, cornerSize),
                       na.rm = TRUE)
#   message(paste0("UR2: ", computeTest(dataUL$data$count, cornerSize, dataUR$data$count, cornerSize)))
#   message(paste0("UR3: ", computeTest(dataBL$data$count, cornerSize, dataUR$data$count, cornerSize)))
#   message(paste0("UR4: ", computeTest(dataBR$data$count, cornerSize, dataUR$data$count, cornerSize)))
#   message(paste0("UR5: ", dataUR$test))
    dataBL$test <- max(dataBL$test,
                       computeTest(dataUL$data$count, cornerSize, dataBL$data$count, cornerSize),
                       computeTest(dataUL$data$count, cornerSize, dataBL$data$count, cornerSize),
                       computeTest(dataBR$data$count, cornerSize, dataBL$data$count, cornerSize),
                       na.rm = TRUE)
    dataBR$test <- max(dataBR$test,
                       computeTest(dataUL$data$count, cornerSize, dataBR$data$count, cornerSize),
                       computeTest(dataUR$data$count, cornerSize, dataBR$data$count, cornerSize),
                       computeTest(dataBL$data$count, cornerSize, dataBR$data$count, cornerSize),
                       na.rm = TRUE)
    counts <- bind_rows(rest,
                        dataUL$data,
                        dataUR$data,
                        dataBL$data,
                        dataBR$data) %>%
        filter(count != 0)
#   message("check 7")
    testPlot <- ggplot(counts, aes(x = type, y = count)) + 
        geom_violin() +
        scale_y_log10()
#   message("check 8")
    object@interactionMatrix <- object@interactionMatrix %>% filter(count > 0)
#   message("check 9")
    mapPlot <- plot.10X2Ref(object, circles = TRUE)
#   message("check 10")
    progressBar$tick()
    return(list(ref1     = object@chromosome1,
                ref2     = object@chromosome2,
                testPlot = testPlot,
                mapPlot  = mapPlot,
                tUL      = dataUL$test,
                tUR      = dataUR$test,
                tBL      = dataBL$test,
                tBR      = dataBR$test))
}

cornerEuclideanDistance <- function(observedDistribution, backgroundDistribution) {
# message(str(observedDistribution))
# message(str(backgroundDistribution))
    n <- nrow(observedDistribution)
    #if (n == 0) return(0)
    v <- observedDistribution %>%
        dplyr::left_join(backgroundDistribution, by = "distance") %>%
        dplyr::mutate(d = (count - meanCount) ^ 2) %>%
        dplyr::pull(d) %>%
        sum() %>%
        sqrt()
    return(v / n)
}

computeCornerDistanceOffset <- function(offset, backgroundDistribution, corner, maxDistance) {
# message(str(offset))
    backgroundDistribution <- backgroundDistribution %>%
        dplyr::mutate(distance = distance - offset) %>%
        dplyr::filter(distance >= 0)
    corner <- corner %>%
        dplyr::filter(distance <= maxDistance - offset)
    cornerEuclideanDistance(corner, backgroundDistribution)
}

# Extract the rignt corner, fill missing values with 0, and summarize as distance/count tibble.
extractCorner <- function(object, vert, hor) {
# message(str(object@interactionMatrix))
# message(str(vert))
# message(str(hor))
    interactionMatrix <- object@interactionMatrix
    if (vert == "B") {
        interactionMatrix <- interactionMatrix %>%
            dplyr::mutate(bin1 = object@size1 - bin1)
    }
    if (hor == "R") {
        interactionMatrix <- interactionMatrix %>%
            dplyr::mutate(bin2 = object@size2 - bin2)
    }
# message(str(interactionMatrix))
    interactionMatrix <- interactionMatrix %>%
        dplyr::mutate(binMin = pmin(bin1, bin2)) %>%
        dplyr::mutate(binMax = pmax(bin1, bin2)) %>%
        dplyr::select(-c(bin1, bin2)) %>%
        dplyr::rename(bin1 = binMax, bin2 = binMin)
# message(str(interactionMatrix))
# message(str(computeDistanceCount(interactionMatrix, object@parameters@maxLinkRange)))
# message(str(computeDistanceCount(interactionMatrix, object@parameters@maxLinkRange) %>% filter(count > 0)))
    computeDistanceCount(interactionMatrix, object@parameters@maxLinkRange)
}

.checkJoins <- function(object, sizes, pvalueThreshold) {
    if (! is(object, "tenxcheckerExp")) {
        stop("Parameter should be a tenxcheckerExp.")
    }
    message(paste0("Checking joins for '", object@name, "'."))
    objects <- splitBy2Ref(object, sizes)
    pb <- progress_bar$new(total = length(objects))
    #joins <- bplapply(objects, checkJoin, progressBar = pb)
    joins <- lapply(objects, checkJoin, progressBar = pb)
    if (length(joins) == 0) {
        message("\tNo join found.")
        joinsObject <- new("tenxcheckerJoins")
        joinsObject@data <- tibble(ref1   = character(),
                                   ref2   = character(),
                                   vert   = character(),
                                   hor    = character(),
                                   pvalue = numeric())
        joinsObject@testPlots <- c()
        joinsObject@mapPlots <- c()
        object@joins <- joinsObject
        return(object)
    }
    joins <- tibble(ref1     = map_chr(joins, "ref1"),
                    ref2     = map_chr(joins, "ref2"),
                    tUL      = map(joins, "tUL"),
                    tUR      = map(joins, "tUR"),
                    tBL      = map(joins, "tBL"),
                    tBR      = map(joins, "tBR"),
                    testPlot = map(joins, "testPlot"),
                    mapPlot  = map(joins, "mapPlot")) %>%
        tidyr::gather(key = "key", value = "test", tUL, tUR, tBL, tBR) %>%
        dplyr::filter(!unlist(map(test, is_null))) %>%
        dplyr::mutate(pvalue = unlist(test)) %>%
        tidyr::separate(key, c(1, 2), into = c(NA, "vert", "hor")) %>%
        dplyr::mutate(vert = factor(vert)) %>%
        dplyr::mutate(hor = factor(hor)) %>%
        dplyr::filter(pvalue <= pvalueThreshold) %>%
        dplyr::arrange(pvalue)
    joinsObject <- new("tenxcheckerJoins")
    joinsObject@data <- joins %>%
        dplyr::select(-c(test, testPlot, mapPlot))
    joinsObject@testPlots <- joins %>%
        dplyr::select(ref1, ref2, testPlot) %>%
        tidyr::unite("name", ref1, ref2) %>%
        deframe()
    joinsObject@mapPlots <- joins %>%
        dplyr::select(ref1, ref2, mapPlot) %>%
        tidyr::unite("name", ref1, ref2) %>%
        deframe()
    object@joins <- joinsObject
    return(object)
}

checkJoins <- function(object, pvalueThreshold) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    object@data <- map(object@data, .checkJoins, sizes = object@sizes, pvalueThreshold)
    return(invisible(object))
}

.removeDuplicateJoins <- function(object) {
    if (! is(object, "tenxcheckerExp")) {
        stop("Parameter should be a tenxcheckerExp.")
    }
    message(paste0("\t\tDataset '", object@name, "'."))
    # Transform joins to before/after
    tmp <- object@joins@data %>%
        dplyr::mutate(hor = as.character(hor)) %>%
        dplyr::mutate(vert = as.character(vert)) %>%
        dplyr::mutate(after1 = (hor == "R")) %>%
        dplyr::mutate(after2 = (vert == "B"))
    ambiguousJoins <- dplyr::bind_rows(tmp %>%
                                           dplyr::select(ref1, after1, pvalue) %>%
                                           dplyr::rename(ref = ref1, after = after1),
                                       tmp %>%
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
    # Transfor to ref1/bin1
    refBin1 <- ambiguousJoins %>%
        dplyr::rename(ref1 = ref) %>%
        dplyr::mutate(hor = dplyr::if_else(after, "R", "L")) %>%
        dplyr::select(ref1, hor)
    # Transfor to ref2/bin2
    refBin2 <- ambiguousJoins %>%
        dplyr::rename(ref2 = ref) %>%
        dplyr::mutate(vert = dplyr::if_else(after, "B", "U")) %>%
        dplyr::select(ref2, vert)
    object@joins@data <- object@joins@data %>%
        dplyr::anti_join(refBin1, by = c("ref1", "hor")) %>%
        dplyr::anti_join(refBin2, by = c("ref2", "vert"))
    return(object)
}

removeDuplicateJoins <- function(object) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    message("\tRemoving ambiguous joins")
    object@data <- map(object@data, .removeDuplicateJoins)
    return(invisible(object))
}

..checkCorners <- function(parameters, object, sizes, pb) {
    objectRef <- extract2Ref(object, parameters$ref1, parameters$ref2, sizes[[parameters$ref1]], sizes[[parameters$ref2]])
    corner    <- extractCorner(objectRef, parameters$vert, parameters$hor)
    values    <- purrr::map(seq.int(from = 0, to = objectRef@parameters@maxLinkRange, by = 1), computeCornerDistanceOffset, object@parameters@distanceCount, corner, object@parameters@maxLinkRange) %>% unlist()
    distance  <- object@parameters@cornerScores %>%
        dplyr::mutate(observed = values) %>%
        dplyr::mutate(difference = observed - score) %>%
        dplyr::filter(difference < 0.0) %>%
        dplyr::filter(distance <= object@parameters@cornerLimit) %>%
        dplyr::slice_min(difference) %>%
        dplyr::pull(distance)
    pb$tick()
    if (length(distance) == 0) return(-1)
    return(distance)
#message(str(values))
}

.checkCorners <- function(object, sizes, pvalueThreshold) {
    message(paste0("\t\tDataset '", object@name, "'."))
    nJoins <- nrow(object@joins@data)
    pb     <- progress_bar$new(total = nrow(object@joins@data))
    values <- object@joins@data %>%
        dplyr::mutate(hor = as.character(hor)) %>%
        dplyr::mutate(vert = as.character(vert)) %>%
        purrr::transpose() %>%
        purrr::map(..checkCorners, object, sizes, pb) %>%
        unlist()
    object@joins@data <- object@joins@data %>%
        dplyr::mutate(offset = values) %>%
        dplyr::filter(offset >= 0)
    message(paste0("\t\t\tKept ", nJoins, "/", nrow(object@joins@data), "."))
    return(object)
}

checkCorners <- function(object, pvalueThreshold) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    message("\tChecking corner shape.")
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
        dplyr::filter(p_adj <= pvalueThreshold)
    message(paste0(nrow(object@joins), " joins found."))
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
                                       dplyr::select(ref1, ref2, hor, vert) %>%
                                       dplyr::mutate(left = hor  == "L") %>%
                                       dplyr::mutate(up   = vert == "U") %>%
                                       purrr::transpose(),
                                   .addJoinPlots, object = object, pb = pb)
    names(object@joinPlots) <- object@joins %>% unite("name", c(ref1, ref2, vert), sep = "_") %>% unite("name", c(name, hor), sep = "") %>% pull(name)
    return(object)                                                                                                                                                                                      
}

findJoins <- function(object, pvalue = 0.05) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    object <- estimateCornerDistanceThreshold(object, pvalue)
    object <- checkJoins(object, pvalue)
    object <- removeDuplicateJoins(object)
    object <- checkCorners(object, pvalue)
    object <- mergeJoins(object, pvalue)
    object <- addJoinPlots(object)
    return(invisible(object))                                                                                                                                                                                      
}
