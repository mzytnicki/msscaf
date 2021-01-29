computeTest <- function(object, xmin, xmax, ymin, ymax, name) {
    corner <- object@interactionMatrix %>%
        filter(bin1 >= xmin,
               bin1 <= xmax,
               bin2 >= ymin,
               bin2 <= ymax) %>%
        select(count) %>%
        mutate(type = name)
    rest <- object@interactionMatrix %>%
        filter(bin1 < xmin | bin1 > xmax | bin2 < ymin | bin2 > ymax) %>%
        select(count) %>%
        mutate(type = "rest")
    tmp <- corner %>% pull(count)
    if (length(tmp) == 0) {
        test <- NULL
    }
    else if (max(tmp) < 10) {
        test <- NULL
    }
    else {
        test <- tryCatch(t.test(rest   %>% pull(count),
                                corner %>% pull(count),
                                alternative = "less"),
                         error = function(c) NULL)
    }
    return(list(data = corner, test = test))
}

checkJoin <- function(object, progressBar) {
    if (! is(object, "tenxchecker2RefExp")) {
        stop("Parameter should be a 'tenxchecker2RefExp'.")
    }
    # message(paste0("  Matrix ",
    #                object@chromosome1,
    #                "/",
    #                object@chromosome2))
    if (object@size1 <= 2 * object@parameters@maxLinkRange) {
        firstLim1 <- object@size1 / 2
        lastLim1  <- object@size1 / 2
    } else {
        firstLim1 <- object@parameters@maxLinkRange
        lastLim1  <- object@size1 - object@parameters@maxLinkRange
    }
    if (object@size2 <= 2 * object@parameters@maxLinkRange) {
        firstLim2 <- object@size2 / 2
        lastLim2  <- object@size2 / 2
    } else {
        firstLim2 <- object@parameters@maxLinkRange
        lastLim2  <- object@size2 - object@parameters@maxLinkRange
    }
    #object <- normalizeKR(object)
    interior <- object@interactionMatrix %>%
        filter((bin1 >= firstLim1 | bin2 >= firstLim2) &
               (bin1 <= lastLim1  | bin2 >= firstLim2) &
               (bin1 >= firstLim1 | bin2 <= lastLim2) &
               (bin1 <= lastLim1  | bin2 <= lastLim2)) %>%
        select(count) %>%
        mutate(type = "interior")
    dataUL <- computeTest(object,
                          0, firstLim1,
                          0, firstLim2,
                          "cornerUL")
    dataUR <- computeTest(object,
                          lastLim1, object@size1,
                          0, firstLim2,
                          "cornerUR")
    dataBL <- computeTest(object,
                          0, firstLim1,
                          lastLim2, object@size2,
                          "cornerBL")
    dataBR <- computeTest(object,
                          lastLim1, object@size1,
                          lastLim2, object@size2,
                          "cornerBR")
    counts <- bind_rows(interior,
                        dataUL$data,
                        dataUR$data,
                        dataBL$data,
                        dataBR$data)
    testPlot <- ggplot(counts, aes(x = type, y = count)) + 
        geom_boxplot() +
        geom_jitter() +
        scale_y_log10()
    mapPlot <- plot.10X2Ref(object,
                          x1 = firstLim1,
                          x2 = lastLim1,
                          y1 = firstLim2,
                          y2 = lastLim2)
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

.checkJoins <- function(object, sizes, pvalueThreshold) {
    if (! is(object, "tenxcheckerExp")) {
        stop("Parameter should be a tenxcheckerExp.")
    }
    message(paste0("Checking joins for '", object@name, "'."))
    objects <- splitBy2Ref(object, sizes)
    pb <- progress_bar$new(total = length(objects))
    joins <- lapply(objects, checkJoin, progressBar = pb)
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
        mutate(pvalue = map_dbl(test, "p.value")) %>%
        separate(key, c(1, 2), into = c(NA, "vert", "hor")) %>%
        mutate(vert = factor(vert)) %>%
        mutate(hor = factor(hor)) %>%
        arrange(pvalue) %>%
        filter(pvalue <= pvalueThreshold)
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

checkJoins <- function(object, pvalue) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    object@data <- map(object@data, .checkJoins, sizes = object@sizes, pvalue)
    return(invisible(object))
}

mergeJoins <- function(object) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    joins <- dplyr::bind_rows(map(map(object@data, "joins"), "data")) %>%
        dplyr::distinct() %>%
        dplyr::arrange(ref1, ref2)
    object@joins <- joins                                                                                                                                                                                        
    return(invisible(object))                                                                                                                                                                                      
}

findJoins <- function(object, pvalue = 0.05) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    object <- checkJoins(object, pvalue)
    object <- mergeJoins(object)
}
