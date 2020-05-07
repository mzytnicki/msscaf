plotTriangles <- function(triangles, bins = NULL) {
    p <- triangles %>%
        gather(key = "type", value = "value", fcMeanCount, nCells) %>%
        ggplot(aes(x = bin, y = value)) +
            geom_line() +
            facet_grid(rows = vars(type), scales = "free_y")
    #if ((length(bins) >= 1) & (!is.null(bins)) & (!is.na(bins))) {
    if (!is.null(bins)) {
        for (bin in bins) {
          p <- p + geom_vline(xintercept = bin,
                              colour = "red",
                              linetype = "longdash")
        }
    }
    return(p)
}

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
    # breakPointBin <- NA
    # breakPoint <- triangles %>%
    #     arrange(fcMeanCount, desc(nCells)) %>%
    #     slice(1) %>%
    #     select(bin, fcMeanCount, nCells) %>%
    #     as.list()
    # if (breakPoint$fcMeanCount <= object@parameters@breakThreshold &
    #     breakPoint$nCells >= object@parameters@breakNCells) {
    #     breakPointBin <- breakPoint$bin
    # }
    plot1 <- plotTriangles(triangles)
    # if (!is.na(breakPointBin)) {
    #     plot1 <- plot1 + geom_vline(xintercept = breakPointBin,
    #                                 colour = "red",
    #                                 linetype = "longdash")
    # }
    plot2 <- triangles %>%
        ggplot(aes(x = fcMeanCount)) +
            geom_histogram()
    plot3 <- plot.10XRef(object, FALSE)
    return(list(data = triangles,
                plot1 = plot1,
                plot2 = plot2,
                plot3 = plot3))
    # return(list(plot1 = plot1,
    #             plot2 = plot2,
    #             plot3 = plot3,
    #             ref   = object@chromosome,
    #             bin   = breakPointBin))
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
    return(list(data = map_dfr(breaks, "data", .id = "ref") %>%
                    mutate(ref = factor(ref)), 
                plot1 = map(breaks, "plot1"),
                plot2 = map(breaks, "plot2"),
                plot3 = map(breaks, "plot3")))
}

filterBreak <- function(reference, object = object, breaks = breaks) {
    bins <- c()
    plot1 <- c()
    objectRef <- extractRef(object, reference)
    breaksRef <- breaks %>%
        filter(ref == reference)
    originalBreakRef <- breaksRef
    repeat {
        breakPointBin <- breaksRef %>%
            filter(fcMeanCount <= object@parameters@breakThreshold) %>%
            filter(bin         >= object@parameters@maxLinkRange) %>%
            filter(bin         <= objectRef@size - object@parameters@maxLinkRange) %>%
            filter(nCells      >= object@parameters@breakNCells) %>%
            arrange(fcMeanCount, desc(nCells)) %>%
            slice(1) %>%
            pull(bin)
        if (length(breakPointBin) == 0) {
            message(bins)
            return(list(ref   = as.character(reference),
                        bins  = bins,
                        plot1 = plotTriangles(originalBreakRef, bins),
                        plot2 = plot.10XRef(objectRef, TRUE, bins)))
        }
        bins <- c(bins, breakPointBin)
        breakPointBin <- breakPointBin[[1]]
        breaksRef %<>%
            filter(abs(breakPointBin - bin) >= object@parameters@maxLinkRange)
    }
}

filterBreaks <- function(object, breaks) {
    selectedRefs <- breaks %>%
        filter(fcMeanCount <= object@parameters@breakThreshold) %>%
        # TODO: find a better filter here?
        #filter(nCells      >= object@parameters@breakNCells) %>%
        select(ref) %>%
        distinct() %>%
        pull()
    if (length(selectedRefs) == 0) {
        message("No break found.")
        return(list(breaks = c(), plots = c()))
    }
    #selectedBreaks <- lapply(selectedRefs, filterBreak, object = object, breaks = breaks)
    selectedBreaks <- bplapply(selectedRefs, filterBreak, object = object, breaks = breaks)
    selectedBreaks <- as_tibble(transpose(selectedBreaks))
    selectedPlots <- selectedBreaks %>%
        select(ref, plot1, plot2) %>%
        mutate(ref = flatten_chr(ref))
    selectedSplits <- selectedBreaks %>%
        select(ref, bins) %>%
        mutate(ref = flatten_chr(ref)) %>%
        unnest(bins) %>%
        rename(bin = bins)
    return(list(breaks = selectedSplits, plots = selectedPlots))
}