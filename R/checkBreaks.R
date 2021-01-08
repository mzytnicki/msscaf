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
        dplyr::filter(bin1 >= value) %>%
        dplyr::filter(bin2 <= value) %>%
        dplyr::summarise(nCells = n(), meanCount = mean(count)) %>%
        tidyr::replace_na(list(meanCount = 0))
}

computeMeanTriangles <- function(data) {
    n <- max(data$bin1, data$bin2)
    if (n < 0) {
        stop("Error while computing triangles: matrix is empty.")
    }
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
    if (nrow(data) == 0) {
       # Matrix is empty, skip
       return(list(data =  NULL,
                   plot1 = NULL,
                   plot2 = NULL,
                   plot3 = NULL))
    }
    triangles <- computeMeanTriangles(data) %>%
        rename(fcMeanCount = meanCount) %>%
        mutate(fcMeanCount = fcMeanCount - median(fcMeanCount))
    # trianglesList <- bplapply(seq.int(object@parameters@nRandomizations),
    # trianglesList <- lapply(seq.int(object@parameters@nRandomizations),
    #                       computeMeanTrianglesRandom,
    #                       data = data)
    # triangles <- bind_rows(trianglesList) %>%
    #     dplyr::select(bin, nCells, fcMeanCount) %>%
    #     group_by(bin) %>%
    #     summarise(nCells = mean(nCells), fcMeanCount = mean(fcMeanCount)) %>%
    #     ungroup()
    # breakPointBin <- NA
    # breakPoint <- triangles %>%
    #     arrange(fcMeanCount, desc(nCells)) %>%
    #     slice(1) %>%
    #     dplyr::select(bin, fcMeanCount, nCells) %>%
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

normalizeAndBreak <- function(object, progressBar, kr, md, diag) {
    progressBar$tick()
    if (isMatrixEmpty(object@interactionMatrix)) {
        return(list(data = NULL, plot1 = NULL, plot2 = NULL, plot3 = NULL))
    }
    #message(paste0("\n  Working on ", object@chromosome, ".\n"))
    if (kr) {
        object <- normalizeKR(object)
    }
    if (md) {
        object <- normalizeMD(object)
    }
    if (diag) {
        object <- removeFarFromDiagonal(object)
    }
    checkBreak(object)
}

checkBreaks <- function(object, kr = TRUE, md = TRUE, diag = TRUE) {
    message("Splitting matrix.")
    objects <- splitByRef(object)
    message("Done.")
    pb <- progress_bar$new(total = length(objects))
    #breaks <- bplapply(objects, normalizeAndBreak, progressBar = pb)
    breaks <- lapply(objects, normalizeAndBreak, progressBar = pb, kr = kr, md = md, diag = diag)
    data   <- map_dfr(breaks, "data", .id = "ref")
    if (nrow(data) == 0) {
        return(list())
    }
    return(list(data = data %>% mutate(ref = factor(ref)), 
                plot1 = map(breaks, "plot1"),
                plot2 = map(breaks, "plot2"),
                plot3 = map(breaks, "plot3")))
}

estimateBreakThreshold <- function(object, breaks, pvalue) {
    object@parameters@breakThreshold <- breaks %>%
        filter(nCells == object@parameters@breakNCells) %>%
        mutate(absFcMeanCount = abs(fcMeanCount)) %>%
        arrange(desc(absFcMeanCount)) %>%
        mutate(class = if_else(fcMeanCount >= 0, 1, 0)) %>%
        mutate(cumSumClass = cumsum(class)) %>%
        mutate(value = cumSumClass / max(cumSumClass)) %>%
        filter(class == 0) %>%
        filter(value <= pvalue) %>%
        tail(n = 1) %>%
        pull(fcMeanCount)
    return(object)
}

filterBreak <- function(parameters) {
    maxNBreaks <- 1000
    bins <- c()
    plot1 <- c()
    reference <- parameters$ref
    objectRef <- parameters$object
    breaksRef <- parameters$breaks
    originalBreakRef <- breaksRef
    breaksRef <- breaksRef %>%
        filter(fcMeanCount <= objectRef@parameters@breakThreshold) %>%
        filter(nCells      >= objectRef@parameters@breakNCells) %>%
        filter(bin         >= objectRef@parameters@maxLinkRange) %>%
        filter(bin         <= objectRef@size - objectRef@parameters@maxLinkRange)
    repeat {
        breakPointBin <- breaksRef %>%
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
        if (length(bins) >= maxNBreaks) {
            stop(paste0("Found more than ",
                        maxNBreaks,
                        " filtered breaks in ref '",
                        as.character(reference),
                        "'.  This is probably too much, please decrease ",
                        "'object@parameters@breakThreshold'.")) 
        }
        breakPointBin <- breakPointBin[[1]]
        breaksRef %<>%
            filter(abs(breakPointBin - bin) >= objectRef@parameters@maxLinkRange)
    }
}

filterBreaks <- function(object, breaks) {
    if (length(breaks) == 0) {
        return(list())
    }
    selectedRefs <- breaks %>%
        dplyr::select(ref) %>%
        distinct() %>%
        pull() %>%
        as.character()
    if (length(selectedRefs) == 0) {
        message("No break found.")
        return(list(breaks = c(), plots = c()))
    }
    splitObject <- splitByRef(object)[selectedRefs]
    splitBreaks <- breaks %>%
        group_by(ref) %>%
        group_split()
    functionParameters <- transpose(list(ref = selectedRefs, object = splitObject, breaks = splitBreaks))
    # parallel seems to need too much RAM
    selectedBreaks <- bplapply(functionParameters, filterBreak)
    #selectedBreaks <- lapply(functionParameters, filterBreak)
    selectedBreaks <- as_tibble(transpose(selectedBreaks))
    if (is.null(unlist(selectedBreaks$bins))) {
        message("No break passed the filter.")
        return(list(breaks = c(), plots = c()))
    }
    selectedPlots <- selectedBreaks %>%
        dplyr::select(ref, plot1, plot2) %>%
        mutate(ref = flatten_chr(ref))
    selectedSplits <- selectedBreaks %>%
        dplyr::select(ref, bins) %>%
        mutate(ref = flatten_chr(ref)) %>%
        unnest(bins) %>%
        rename(bin = bins)
    return(list(breaks = selectedSplits, plots = selectedPlots))
}
