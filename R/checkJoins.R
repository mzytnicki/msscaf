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

checkJoin <- function(object) {
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
    plot1 <- ggplot(counts, aes(x = type, y = count)) + 
        geom_boxplot() +
        geom_jitter() +
        scale_y_log10()
    plot2 <- plot.10X2Ref(object,
                          x1 = firstLim1,
                          x2 = lastLim1,
                          y1 = firstLim2,
                          y2 = lastLim2)
    return(list(ref1  = object@chromosome1,
                ref2  = object@chromosome2,
                plot1 = plot1,
                plot2 = plot2,
                tUL   = dataUL$test,
                tUR   = dataUR$test,
                tBL   = dataBL$test,
                tBR   = dataBR$test))
}

checkJoins <- function(object) {
    message("Splitting matrix.")
    objects <- splitBy2Ref(object)
    joins <- lapply(objects, checkJoin)
    tibble(ref1  = map_chr(joins, "ref1"),
           ref2  = map_chr(joins, "ref2"),
           tUL   = map(joins, "tUL"),
           tUR   = map(joins, "tUR"),
           tBL   = map(joins, "tBL"),
           tBR   = map(joins, "tBR"),
           plot1 = map(joins, "plot1"),
           plot2 = map(joins, "plot2")) %>%
        gather(key = "key", value = "test", tUL, tUR, tBL, tBR) %>%
        filter(!unlist(map(test, is_null))) %>%
        mutate(pvalue = map_dbl(test, "p.value")) %>%
        separate(key, c(1, 2), into = c(NA, "vert", "hor")) %>%
        mutate(vert = factor(vert)) %>%
        mutate(hor = factor(hor))
}

filterJoins <- function(object, joins) {
    joins %>%
        filter(pvalue <= object@parameters@pvalueThreshold)
}