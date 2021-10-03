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

computeMeanTriangles <- function(data, distance) {
    n <- max(data$bin1, data$bin2)
    if (n < 0) {
        stop("Error while computing triangles: matrix is empty.")
    }
    computeMeanTrianglesCpp(data, distance) %>% as_tibble()
}


checkBreak <- function(object) {
    #message("    Checking breaks.")
    data <- object@interactionMatrix
    if (nrow(data) == 0) {
       # Matrix is empty, skip
       return(list(data           = tibble(ref = character(), fcMeanCount = numeric(), nCells = integer()),
                   changePlot     = NULL,
                   changeDistPlot = NULL,
                   mapPlot        = NULL))
    }
    triangles      <- computeMeanTriangles(data, object@parameters@maxLinkRange)
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
        return(list(data           = tibble(ref = character(), fcMeanCount = numeric(), nCells = integer()),
                    changePlot     = NULL,
                    changeDistPlot = NULL,
                    mapPlot        = NULL))
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
#   if (object@name == "HiC") {
#       kr = TRUE
#   }
    objects <- splitByRef(object, chromosomes, sizes)
    pb <- progress_bar$new(total = length(objects))
    message("\t\tComputing stats.")
    breaks <- lapply(objects, normalizeAndBreak, progressBar = pb, kr = kr, md = md, diag = diag)
    #message(str(breaks))
    data   <- map_dfr(breaks, "data", .id = "ref") %>% dplyr::mutate(ref = factor(ref, levels = chromosomes))
    breaksObject                 <- new("tenxcheckerBreaks")
    breaksObject@data            <- data
    breaksObject@changePlots     <- map(breaks, "changePlot")
    breaksObject@changeDistPlots <- map(breaks, "changeDistPlot")
    breaksObject@mapPlots        <- map(breaks, "mapPlot")
    object@breaks                <- breaksObject
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

.computeBreakPvalue <- function(object, pvalue) {
    if (! is(object, "tenxcheckerExp")) {
        stop(paste0("Parameter should be a tenxcheckerExp, it is a ", is(object), " ."))
    }
    message(paste0("\tDataset '", object@name , "'."))
    # Take half of the 9th decile as the minimum number of cells.
    object@parameters@breakNCells <- object@breaks@data %>%
        dplyr::slice_min(nCells, prop = 0.9) %>%
        dplyr::slice_max(nCells, n = 1, with_ties = FALSE) %>%
        dplyr::pull(nCells)
    object@parameters@breakNCells <- object@parameters@breakNCells / 2
    object@breaks@data <- object@breaks@data %>%
        dplyr::filter(nCells >= object@parameters@breakNCells) %>%
        dplyr::mutate(absFcMeanCount = abs(fcMeanCount)) %>%
        dplyr::arrange(desc(absFcMeanCount)) %>%
        dplyr::mutate(class = if_else(fcMeanCount >= 0, 1, 0)) %>%
        dplyr::mutate(cumSumClass = cumsum(class)) %>%
        dplyr::mutate(pvalue = cumSumClass / max(cumSumClass)) %>%
        dplyr::select(ref, bin, nCells, fcMeanCount, pvalue)
    return(object)
}

computeBreakPvalue <- function(object) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    message("Estimating break thresholds")
    object@data <- map(object@data, .computeBreakPvalue)
    return(invisible(object))
}

removeNearEqualBreaks <- function(reference, breaks, distance, size, pvalueThreshold) {
    maxNBreaks <- 10000
    breaks     <- breaks %>%
        dplyr::filter(fcMeanCount <  0) %>%
        dplyr::filter(pvalue      <= pvalueThreshold) %>%
        #dplyr::filter(nCells      >= object@parameters@breakNCells) %>%
        dplyr::filter(bin         >= distance) %>%
        dplyr::filter(bin         <= size - distance)
    outputBreaks <- breaks %>% slice_head(n = 0)
    repeat {
        if (nrow(breaks) == 0) {
            return(outputBreaks)
        }
        breakPoint <- breaks %>%
            dplyr::slice_min(pvalue, n = 1, with_ties = FALSE)
        if (nrow(outputBreaks) >= maxNBreaks) {
            stop(paste0("Found more than ",
                        maxNBreaks,
                        " filtered breaks in ref '",
                        as.character(reference),
                        "'.  This is probably too much, please decrease ",
                        "'object@parameters@breakThreshold'.")) 
        }
        outputBreaks <- outputBreaks %>% add_row(breakPoint)
        breaks %<>%
            dplyr::filter(abs(breakPoint$bin[[1]] - bin) > distance)
    }
}

filterBreak <- function(parameters, pvalueThreshold, pb) {
    reference        <- parameters$ref
    objectRef        <- parameters$object
    breaksRef        <- parameters$breaks
    originalBreakRef <- breaksRef
    breaks           <- removeNearEqualBreaks(reference, breaksRef, objectRef@parameters@maxLinkRange, objectRef@size, pvalueThreshold)
    pb$tick()
    return(list(data       = breaks,
                changePlot = plotTriangles(originalBreakRef),
                mapPlot    = plot.10XRef(objectRef, TRUE, breaks$bin)))
}

.filterBreaks <- function(object, pvalue, chromosomes, sizes) {
    if (! is(object, "tenxcheckerExp")) {
        stop(paste0("Parameter should be a tenxcheckerExp, it is a ", is(object), " ."))
    }
    message(paste0("\tDataset '", object@name, "'..."))
    object@breaks@filteredData <- object@breaks@data %>% slice_head(n = 0)
    object@breaks@changePlots  <- c()
    object@breaks@mapPlots     <- c()
    if (nrow(object@breaks@data) == 0) {
        message("\t\tNo break found.")
        return(object)
    }
    splitBreaks <- object@breaks@data %>% dplyr::group_by(ref)
    splitNames  <- splitBreaks %>% dplyr::group_keys() %>% dplyr::pull(ref) %>% as.character()
    splitBreaks <- splitBreaks %>% dplyr::group_split()
    splitObject <- splitByRef(object, chromosomes, sizes)[splitNames]
    functionParameters <- transpose(list(ref = splitNames, object = splitObject, breaks = splitBreaks))
    # parallel seems to need too much RAM
    #selectedBreaks <- bplapply(functionParameters, filterBreak)
    pb             <- progress_bar$new(total = length(splitNames))
    selectedBreaks <- lapply(functionParameters, filterBreak, pvalueThreshold = pvalue, pb = pb)
    selectedBreaks <- transpose(selectedBreaks)
    if (length(selectedBreaks$data) == 0) {
        message("\t\tNo break passed the filter.")
        return(object)
    }
    object@breaks@filteredData <- bind_rows(selectedBreaks$data)
    object@breaks@changePlots  <- setNames(selectedBreaks$changePlot, splitNames)
    object@breaks@mapPlots     <- setNames(selectedBreaks$mapPlot, splitNames)
    message(paste0("\t\t... filtered to ", nrow(object@breaks@filteredData), " breaks."))
    return(object)
}

filterBreaks <- function(object, pvalue) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    message("Filtering results.")
    object@data <- map(object@data, .filterBreaks, pvalue = pvalue, chromosomes = object@chromosomes, sizes = object@sizes)
    return(invisible(object))
}

..compareBreaks <- function(object1, object2) {
    if (! is(object1, "tenxcheckerExp")) {
        stop("Parameter 1 should be a tenxcheckerExp.")
    }
    if (! is(object2, "tenxcheckerExp")) {
        stop("Parameter 2 should be a tenxcheckerExp.")
    }
    if (object1@name == object2@name) {
        return(object1)
    }
    if ((nrow(object1@breaks@filteredData) == 0) |
        (nrow(object2@breaks@filteredData) == 0)) {
        return(object1)
    } 
    newData <- object1@breaks@filteredData %>%
        dplyr::left_join(object2@breaks@data, by = c("ref", "bin"), suffix = c("", "_other")) %>%
        dplyr::filter(fcMeanCount_other <= 0) %>%
        dplyr::filter(pvalue_other < 0.5) %>%
        #dplyr::filter(nCells < object2@parameters@breakNCells | fcMeanCount <= 0) %>%
        dplyr::select(ref, bin, nCells, fcMeanCount, pvalue)
    object1@breaks@filteredData <- newData
    message(paste0("\tDataset '", object1@name, "', compared to ", object2@name, ", is filtered down to ", nrow(object1@breaks@filteredData), " breaks."))
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
    message("Comparing breaks")
    object@data <- purrr::map(object@data, .compareBreaks, objects = object@data)
    return(invisible(object))
}

mergeBreaks <- function(object, pvalueThreshold) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    minMaxLinkRange <- min(unlist(map(map(object@data, "parameters"), "maxLinkRange")))
    breaks          <- dplyr::bind_rows(map(map(object@data, "breaks"), "filteredData")) %>% group_by(ref)
    refs            <- dplyr::group_keys(breaks) %>% tibble::deframe()
    breaks          <- dplyr::group_split(breaks)
    sizes           <- object@sizes[refs]
    nGroups         <- length(refs)
    minMaxLinkRange <- rep(minMaxLinkRange, nGroups)
    pvalueThreshold <- rep(pvalueThreshold, nGroups)
    breaks          <- pmap(list(refs, breaks, minMaxLinkRange, sizes, pvalueThreshold), removeNearEqualBreaks)
    object@breaks   <- breaks %>% bind_rows()
    if (nrow(object@breaks) == 0) {
        object@breaks <- dplyr::slice_head(object@data[[1]]@breaks@data, n = 0)
    }
    object@breaks   <- object@breaks %>% dplyr::arrange(ref, bin)
    return(invisible(object))
}

..addBreakPlots <- function(object, ref, bin, size) {
    objectRef <- extractRef(object, ref, size)
    zoomSize <- 100
    return(plot.10XRef(objectRef, bins = c(bin), lim = list(bin - zoomSize, bin + zoomSize)))
}

.addBreakPlots <- function(parameters, object, pb) {
    ref               <- parameters$ref
    bin               <- parameters$bin
    breakPlots        <- map(object@data, ..addBreakPlots, ref = ref, bin = bin, size = object@sizes[[ref]])
    names(breakPlots) <- map(object@data, "name")
    pb$tick()
    return(breakPlots)
}

addBreakPlots <- function(object) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    message(paste0("Computing plots."))
    pb                       <- progress_bar$new(total = nrow(object@breaks))
    object@breakPlots        <- map(object@breaks %>% dplyr::mutate(ref = as.character(ref)) %>% transpose(), .addBreakPlots, object = object, pb = pb)
    names(object@breakPlots) <- object@breaks %>% unite("name", c(ref, bin), sep = "_") %>% pull(name)
    return(object)                                                                                                                                                                                      
}

findBreaks <- function(object, pvalue = 0.05) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    object <- checkBreaks(object)
    object <- computeBreakPvalue(object)
    object <- filterBreaks(object, pvalue)
    object <- compareBreaks(object)
    object <- mergeBreaks(object, pvalue)
    if (nrow(object@breaks) <= 1000) object <- addBreakPlots(object)
    gc(verbose = FALSE)
    return(invisible(object))
}
