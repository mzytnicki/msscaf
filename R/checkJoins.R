checkJoin <- function(object) {
    message(paste0("  Matrix ",
                   object@chromosome1,
                   "/",
                   object@chromosome2))
    if (object@size1 <= 2 * object@parameters@maxLinkRange) {
        firstLim1 <- object@size1 / 2
        lastLim1  <- object@size1 / 2
    }
    else {
        firstLim1 <- object@parameters@maxLinkRange
        lastLim1  <- object@size1 - object@parameters@maxLinkRange
    }
    if (object@size2 <= 2 * object@parameters@maxLinkRange) {
        firstLim2 <- object@size2 / 2
        lastLim2  <- object@size2 / 2
    }
    else {
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
    cornerUL <- object@interactionMatrix %>%
        filter(bin1 <= firstLim1,
               bin2 <= firstLim2) %>%
        select(count) %>%
        mutate(type = "cornerUL")
    cornerUR <- object@interactionMatrix %>%
        filter(bin1 >= lastLim1,
               bin2 <= firstLim2) %>%
        select(count) %>%
        mutate(type = "cornerUR")
    cornerBL <- object@interactionMatrix %>%
        filter(bin1 <= firstLim1,
               bin2 >= lastLim2) %>%
        select(count) %>%
        mutate(type = "cornerBL")
    cornerBR <- object@interactionMatrix %>%
        filter(bin1 >= lastLim1,
               bin2 >= lastLim2) %>%
        select(count) %>%
        mutate(type = "cornerBR")
    counts <- bind_rows(interior,
                        cornerUL,
                        cornerUR,
                        cornerBL,
                        cornerBR)
    plot1 <- ggplot(counts, aes(x = type, y = count)) + 
        geom_boxplot() +
        geom_jitter() +
        scale_y_log10()
    plot2 <- plot.10X2Ref(object)
    tUL <- t.test(interior %>% pull(count), cornerUL %>% pull(count), alternative = "less")
    tUR <- t.test(interior %>% pull(count), cornerUR %>% pull(count), alternative = "less")
    tBL <- t.test(interior %>% pull(count), cornerBL %>% pull(count), alternative = "less")
    tBR <- t.test(interior %>% pull(count), cornerBR %>% pull(count), alternative = "less")
    return(list(ref1  = object@chromosome1,
                ref2  = object@chromosome2,
                plot1 = plot1,
                plot2 = plot2,
                tUL   = tUL,
                tUR   = tUR,
                tBL   = tBL,
                tBR   = tBR))
}

checkJoins <- function(object) {
    message("Splitting matrix.")
    objects <- splitBy2Ref(object)
    joins <- bplapply(objects, checkJoin)
    tibble(ref1  = map_chr(joins, "ref1"),
           ref2  = map_chr(joins, "ref2"),
           tUL   = map(joins, "tUL"),
           tUR   = map(joins, "tUR"),
           tBL   = map(joins, "tBL"),
           tBR   = map(joins, "tBR"),
           plot1 = map(joins, "plot1"),
           plot2 = map(joins, "plot2"))
}