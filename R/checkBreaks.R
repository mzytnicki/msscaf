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
        mutate(randMeanCount = computeMeanTriangles(dataRandom)$meanCount) %>%
        mutate(fcMeanCount = meanCount - randMeanCount)
}

checkBreak <- function(object) {
    #message("    Checking breaks.")
    data <- object@interactionMatrix
    if (nrow(data) == 0) {
       # Matrix is empty, skip
       return(list(data           = NULL,
                   changePlot     = NULL,
                   changeDistPlot = NULL,
                   mapPlot        = NULL))
    }
    triangles <- computeMeanTriangles(data) %>%
        rename(fcMeanCount = meanCount) %>%
        mutate(fcMeanCount = fcMeanCount - median(fcMeanCount))
    changePlot     <- plotTriangles(triangles)
    changeDistPlot <- triangles %>%
        ggplot(aes(x = fcMeanCount)) +
            geom_histogram()
    mapPlot <- plot.10XRef(object, FALSE)
    return(list(data           = triangles,
                changePlot     = changePlot,
                changeDistPlot = changeDistPlot,
                mapPlot        = mapPlot))
}

normalizeAndBreak <- function(object, progressBar, kr, md, diag) {
    progressBar$tick()
    if (isMatrixEmpty(object@interactionMatrix)) {
        return(list(data = NULL, changePlot = NULL, changeDistPlot = NULL, mapPlot = NULL))
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

.checkBreaks <- function(object, chromosomes, sizes) {
    message(paste0("\tDataset '", object@name , "'.\n\t\tSplitting matrix."))
    md   = FALSE
    diag = FALSE
    kr   = FALSE
    if (object@name == "HiC") {
        kr = TRUE
    }
    objects <- splitByRef(object, chromosomes, sizes)
    pb <- progress_bar$new(total = length(objects))
    #breaks <- bplapply(objects, normalizeAndBreak, progressBar = pb)
    message("\t\tComputing stats.")
    breaks <- lapply(objects, normalizeAndBreak, progressBar = pb, kr = kr, md = md, diag = diag)
    data   <- map_dfr(breaks, "data", .id = "ref") %>% mutate(ref = factor(ref))
    if (nrow(data) == 0) {
        return(list())
    }
    breaksObject                 <- new("tenxcheckerBreaks")
    breaksObject@data            <- data
    breaksObject@changePlots     <- map(breaks, "changePlot")
    breaksObject@changeDistPlots <- map(breaks, "changeDistPlot")
    breaksObject@mapPlots        <- map(breaks, "mapPlot")
    object@breaks <- breaksObject
    return(object)
}

checkBreaks <- function(object) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    message("Finding statistics.")
    object@data <- map(object@data, .checkBreaks, chromosomes = object@chromosomes, sizes = object@sizes)
    return(invisible(object))
}

.estimateBreakThreshold <- function(object, pvalue) {
    if (! is(object, "tenxcheckerExp")) {
        stop(paste0("Parameter should be a tenxcheckerExp, it is a ", is(object), " ."))
    }
    message(paste0("\tDataset '", object@name , "'."))
    # Take half of the 9th decile as the minimum number of cells.
    object@parameters@breakNCells <- object@breaks@data %>%
        dplyr::slice_max(nCells, prop = 0.9) %>%
        dplyr::arrange(nCells) %>%
        dplyr::slice_head(n = 1) %>%
        dplyr::pull(nCells)
    object@parameters@breakNCells <- object@parameters@breakNCells / 2
    object@parameters@breakThreshold <- object@breaks@data %>%
        dplyr::filter(nCells >= object@parameters@breakNCells) %>%
        dplyr::mutate(absFcMeanCount = abs(fcMeanCount)) %>%
        dplyr::arrange(desc(absFcMeanCount)) %>%
        dplyr::mutate(class = if_else(fcMeanCount >= 0, 1, 0)) %>%
        dplyr::mutate(cumSumClass = cumsum(class)) %>%
        dplyr::mutate(value = cumSumClass / max(cumSumClass)) %>%
        dplyr::filter(class == 0) %>%
        dplyr::filter(value <= pvalue) %>%
        dplyr::slice_tail(n = 1) %>%
        dplyr::pull(fcMeanCount)
    return(object)
}

estimateBreakThreshold <- function(object, pvalue) {
    if (! is(object, "tenxcheckerClass")) {                                                                                                                                                                        
        stop("Parameter should be a tenxcheckerClass.")                                                                                                                                                            
    }                                                                                                                                                                                                              
    message("Estimating break thresholds")
    object@data <- map(object@data, .estimateBreakThreshold, pvalue = pvalue)
    return(invisible(object))
}

filterBreak <- function(parameters) {
    maxNBreaks       <- 10000
    bins             <- c()
    reference        <- parameters$ref
    objectRef        <- parameters$object
    breaksRef        <- parameters$breaks
    originalBreakRef <- breaksRef
    breaksRef <- breaksRef %>%
        dplyr::filter(fcMeanCount <= objectRef@parameters@breakThreshold) %>%
        dplyr::filter(nCells      >= objectRef@parameters@breakNCells) %>%
        dplyr::filter(bin         >= objectRef@parameters@maxLinkRange) %>%
        dplyr::filter(bin         <= objectRef@size - objectRef@parameters@maxLinkRange)
    repeat {
        breakPointBin <- breaksRef %>%
            arrange(fcMeanCount, desc(nCells)) %>%
            slice(1) %>%
            pull(bin)
        if (length(breakPointBin) == 0) {
            return(list(ref        = as.character(reference),
                        bins       = bins,
                        changePlot = plotTriangles(originalBreakRef, bins),
                        mapPlot    = plot.10XRef(objectRef, TRUE, bins)))
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
            filter(abs(breakPointBin - bin) > objectRef@parameters@maxLinkRange)
    }
}

.filterBreaks <- function(object, chromosomes, sizes) {
    if (! is(object, "tenxcheckerExp")) {
        stop(paste0("Parameter should be a tenxcheckerExp, it is a ", is(object), " ."))
    }
    object@breaks@filteredData <- tibble(ref = factor(c(), levels = object@breaks@data), bin = integer())
    object@breaks@changePlots  <- c()
    object@breaks@mapPlots     <- c()
    message(paste0("\tDataset '", object@name , "': Filtering breaks."))
    if (nrow(object@breaks@data) == 0) {
        return(object)
    }
    selectedRefs <- object@breaks@data %>%
        dplyr::select(ref) %>%
        distinct() %>%
        pull() %>%
        as.character()
    if (length(selectedRefs) == 0) {
        message("\t\tNo break found.")
        return(object)
    }
    splitObject <- splitByRef(object, chromosomes, sizes)[selectedRefs]
    splitBreaks <- object@breaks@data %>%
        group_by(ref) %>%
        group_split()
    functionParameters <- transpose(list(ref = selectedRefs, object = splitObject, breaks = splitBreaks))
    # parallel seems to need too much RAM
    selectedBreaks <- bplapply(functionParameters, filterBreak)
    #selectedBreaks <- lapply(functionParameters, filterBreak)
    selectedBreaks <- as_tibble(transpose(selectedBreaks))
    if (is.null(unlist(selectedBreaks$bins))) {
        message("\t\tNo break passed the filter.")
        return(object)
    }
    selectedBreaks <- selectedBreaks %>%
        mutate(ref = flatten_chr(ref))
    changePlots <- selectedBreaks %>%
        dplyr::select(ref, changePlot) %>%
        deframe()
    mapPlots <- selectedBreaks %>%
        dplyr::select(ref, mapPlot) %>%
        deframe()
    selectedSplits <- selectedBreaks %>%
        dplyr::select(ref, bins) %>%
        unnest(bins) %>%
        rename(bin = bins)
    object@breaks@filteredData <- selectedSplits
    object@breaks@changePlots  <- changePlots
    object@breaks@mapPlots     <- mapPlots
    return(object)
}

filterBreaks <- function(object) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    message("Filtering results.")
    object@data <- map(object@data, .filterBreaks, chromosomes = object@chromosomes, sizes = object@sizes)
    return(invisible(object))
}

..compareBreaks <- function(object1, object2) {
    if (! is(object1, "tenxcheckerExp")) {
        stop("Parameter 1 should be a tenxcheckerExp.")
    }
    if (! is(object2, "tenxcheckerExp")) {
        stop("Parameter 2 should be a tenxcheckerExp.")
    }
    if ((nrow(object1@breaks@filteredData) == 0) |
        (nrow(object2@breaks@filteredData) == 0)) {
        return(object1)
    } 
    newData <- object1@breaks@filteredData %>%
        dplyr::left_join(object2@breaks@data, by = c("ref", "bin")) %>%
        dplyr::filter(nCells < object2@parameters@breakNCells | fcMeanCount <= 0) %>%
        dplyr::select("ref", "bin")
    object1@breaks@filteredData <- newData
    return(object1)
}

.compareBreaks <- function(object1, objects) {
    if (! is(object1, "tenxcheckerExp")) {
        stop("Parameter 1 should be a tenxcheckerExp.")
    }
    object1 <- purrr::reduce(objects, ..compareBreaks, .init = object1)
    return(object1)
}

compareBreaks <- function(object) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    object@data <- map(object@data, .compareBreaks, objects = object@data)
    return(invisible(object))
}

mergeBreaks <- function(object) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    breaks <- dplyr::bind_rows(map(map(object@data, "breaks"), "filteredData")) %>%
        dplyr::distinct() %>%
        dplyr::arrange(ref, desc(bin))
    object@breaks <- breaks
    return(invisible(object))
}

findBreaks <- function(object, pvalue = 0.05) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    object <- checkBreaks(object)
    object <- estimateBreakThreshold(object, pvalue)
    object <- filterBreaks(object)
    object <- compareBreaks(object)
    object <- mergeBreaks(object)
}
