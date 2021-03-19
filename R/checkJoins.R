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
    if ((length(background) == 0) & (length(data) == 0)) {
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
#   message(paste0("back: ", background))
#   message(paste0("data: ", data))
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

.checkJoins <- function(object, sizes) {
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
        gather(key = "key", value = "test", tUL, tUR, tBL, tBR) %>%
        filter(!unlist(map(test, is_null))) %>%
        mutate(pvalue = unlist(test)) %>%
        separate(key, c(1, 2), into = c(NA, "vert", "hor")) %>%
        mutate(vert = factor(vert)) %>%
        mutate(hor = factor(hor)) %>%
        arrange(pvalue)
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

checkJoins <- function(object) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    object@data <- map(object@data, .checkJoins, sizes = object@sizes)
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

..addJoinPlots <- function(object, ref1, ref2, size1, size2) {
    object2Ref <- extract2Ref(object, ref1, ref2, size1, size2)
    return(plot.10X2Ref(object2Ref))
}


.addJoinPlots <- function(parameters, object, pb) {
    ref1             <- parameters$ref1
    ref2             <- parameters$ref2
    joinPlots        <- map(object@data, ..addJoinPlots, ref1 = ref1, ref2 = ref2, size1 = object@sizes[[ref1]], size2 = object@sizes[[ref2]])
    names(joinPlots) <- map(object@data, "name")
    pb$tick()
    return(joinPlots)
}

addJoinPlots <- function(object) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    message(paste0("Computing plots."))
    pb                      <- progress_bar$new(total = nrow(object@joins))
    object@joinPlots        <- map(object@joins %>% select(ref1, ref2) %>% transpose(), .addJoinPlots, object = object, pb = pb)
    names(object@joinPlots) <- object@joins %>% unite("name", c(ref1, ref2, vert), sep = "_") %>% unite("name", c(name, hor), sep = "") %>% pull(name)
    return(object)                                                                                                                                                                                      
}

findJoins <- function(object, pvalue = 0.05) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    object <- checkJoins(object)
    object <- mergeJoins(object, pvalue)
    object <- addJoinPlots(object)
    return(invisible(object))                                                                                                                                                                                      
}
