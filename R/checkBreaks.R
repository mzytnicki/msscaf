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
    #message("    Checking breaks.")
    data <- object@interactionMatrix
    trianglesList <- bplapply(seq.int(object@parameters@nRandomizations),
    #trianglesList <- lapply(seq.int(object@parameters@nRandomizations),
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

normalizeAndBreak <- function(object, progressBar) {
    progressBar$tick()
    if (isMatrixEmpty(object@interactionMatrix)) {
        return(list(data = NULL, plot1 = NULL, plot2 = NULL, plot3 = NULL))
    }
    #message(paste0("  Working on ", object@chromosome, "."))
    object <- normalizeKR(object)
    object <- normalizeMD(object)
    object <- removeFarFromDiagonal(object)
    checkBreak(object)
}

checkBreaks <- function(object) {
    message("Splitting matrix.")
    objects <- splitByRef(object)
    message("Done.")
    pb <- progress_bar$new(total = length(objects))
    breaks <- lapply(objects, normalizeAndBreak, progressBar = pb)
    if (length(breaks) == 0) {
        return(list())
    }
    return(list(data = map_dfr(breaks, "data", .id = "ref") %>%
                    mutate(ref = factor(ref)), 
                plot1 = map(breaks, "plot1"),
                plot2 = map(breaks, "plot2"),
                plot3 = map(breaks, "plot3")))
}

filterBreak <- function(parameters) {
    bins <- c()
    plot1 <- c()
    reference <- parameters$ref
    objectRef <- parameters$object
    breaksRef <- parameters$breaks
    originalBreakRef <- breaksRef
    repeat {
        breakPointBin <- breaksRef %>%
            filter(fcMeanCount <= objectRef@parameters@breakThreshold) %>%
            filter(bin         >= objectRef@parameters@maxLinkRange) %>%
            filter(bin         <= objectRef@size - objectRef@parameters@maxLinkRange) %>%
            filter(nCells      >= objectRef@parameters@breakNCells) %>%
            arrange(fcMeanCount, desc(nCells)) %>%
            slice(1) %>%
            pull(bin)
        if (length(breakPointBin) == 0) {
            return(list(ref   = as.character(reference),
                        bins  = bins,
                        plot1 = plotTriangles(originalBreakRef, bins),
                        plot2 = plot.10XRef(objectRef, TRUE, bins)))
        }
        bins <- c(bins, breakPointBin)
        breakPointBin <- breakPointBin[[1]]
        breaksRef %<>%
            filter(abs(breakPointBin - bin) >= objectRef@parameters@maxLinkRange)
    }
}

filterBreaks <- function(object, breaks) {
    if (length(breaks) == 0) {
        return(list())
    }
    selectedBreaks <- breaks %>%
        filter(fcMeanCount <= object@parameters@breakThreshold)
        # TODO: find a better filter here?
        #filter(nCells      >= object@parameters@breakNCells) %>%
    selectedRefs <- selectedBreaks %>%
        select(ref) %>%
        distinct() %>%
        pull() %>%
        as.character()
    if (length(selectedRefs) == 0) {
        message("No break found.")
        return(list(breaks = c(), plots = c()))
    }
    splitObject <- splitByRef(object)[selectedRefs]
    splitBreaks <- selectedBreaks %>%
        group_by(ref) %>%
        group_split()
    functionParameters <- transpose(list(ref = selectedRefs, object = splitObject, breaks = splitBreaks))
    #selectedBreaks <- lapply(selectedRefs, filterBreak, object = object, breaks = breaks)
    # parallel seems to need too much RAM
    selectedBreaks <- bplapply(functionParameters, filterBreak)
    selectedBreaks <- as_tibble(transpose(selectedBreaks))
    if (is.null(unlist(selectedBreaks$bins))) {
        message("No break passed the filter.")
        return(list(breaks = c(), plots = c()))
    }
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
