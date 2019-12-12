randomizeData <- function(data) {
    data %>% mutate(count = sample(count))
}

computeMeanTriangle <- function(value, data = data) {
    data %>%
        filter(bin1 <= value) %>%
        filter(bin2 >= value) %>%
        summarise(nCells = n(), meanCount = mean(count))
}

computeMeanTriangles <- function(data) {
    n <- max(data$bin1, data$bin2)
    map_df(seq(n), computeMeanTriangle, data = data) %>%
        rowid_to_column(var = "bin")
}

computeMeanTrianglesRandom <- function(i, data = data) {
    dataRandom <- randomizeData(data)
    computeMeanTriangles(data) %>%
        add_column(randMeanCount = computeMeanTriangles(dataRandom)$meanCount) %>%
        mutate(fcMeanCount = meanCount - randMeanCount)
}

checkBreak <- function(object) {
    message("    Checking breaks.")
    data <- object@interactionMatrix
    trianglesList <- bplapply(seq.int(object@parameters@nRandomizations),
                          computeMeanTrianglesRandom,
                          data = data)
    triangles <- bind_rows(trianglesList) %>%
        select(bin, nCells, fcMeanCount) %>%
        group_by(bin) %>%
        summarise(nCells = mean(nCells), fcMeanCount = mean(fcMeanCount)) %>%
        ungroup()
    breakPointBin <- NA
    breakPoint <- triangles %>%
        arrange(fcMeanCount, desc(nCells)) %>%
        slice(1) %>%
        select(bin, fcMeanCount, nCells) %>%
        as.list()
    if (breakPoint$fcMeanCount <= object@parameters@breakThreshold &
        breakPoint$nCells >= object@parameters@breakNCells) {
        breakPointBin <- breakPoint$bin
    }
    plot1 <- triangles %>%
        gather(key = "type", value = "value", fcMeanCount, nCells) %>%
        ggplot(aes(x = bin, y = value)) +
            geom_line() +
            facet_grid(rows = vars(type), scales = "free_y")
    if (!is.na(breakPointBin)) {
        plot1 <- plot1 + geom_vline(xintercept = breakPointBin,
                                    colour = "red",
                                    linetype = "longdash")
    }
    plot2 <- triangles %>%
        ggplot(aes(x = fcMeanCount)) +
            geom_histogram()
    plot3 <- plot.10XRef(object, FALSE, breakPointBin)
    return(list(plot1 = plot1,
                plot2 = plot2,
                plot3 = plot3,
                ref   = object@chromosome,
                bin   = breakPointBin))
}

normalizeAndBreak <- function(object) {
    message(paste0("  Working on ", object@chromosome, "."))
    object <- normalizeKR(object)
    object <- normalizeMD(object)
    object <- removeFarFromDiagonal(object)
    checkBreak(object)
}

checkBreaks <- function(object) {
    message("Splitting matrix.")
    objects <- splitByRef(object)
    message("Done.")
    breaks <- lapply(objects, normalizeAndBreak)
    tibble(ref   = map_chr(breaks, "ref"),
           bin   = map_dbl(breaks, "bin"), 
           plot1 = map(breaks, "plot1"),
           plot2 = map(breaks, "plot2"),
           plot3 = map(breaks, "plot3"))
}