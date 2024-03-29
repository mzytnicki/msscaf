.checkBreaks <- function(object, chromosomes, sizes) {
    message(paste0("\tDataset '", object@name , "'.\n\t\tComputing stats."))
    breaksObject                 <- new("msscafBreaks")
    breaksObject@data            <- computeMeanTrianglesCpp(object@interactionMatrix, object@parameters@minLinkRange, object@parameters@metaSize, sizes, object@outlierBins) %>%
                                        tibble::as_tibble()
    # Possibly rescale by ref
    factors <- breaksObject@data %>%
        dplyr::filter(nCells >= object@parameters@breakNCells) %>%
        dplyr::select(ref, fcMeanCount) %>%
        dplyr::group_by(ref) %>%
        dplyr::summarize(factor = mean(fcMeanCount), .groups = "drop") %>%
        dplyr::mutate(factor = dplyr::if_else(factor > 0, 0, factor))
    breaksObject@data <- breaksObject@data %>%
        dplyr::left_join(factors, by = "ref") %>%
        dplyr::mutate(fcMeanCount = fcMeanCount - factor) %>%
        dplyr::select(- factor)
    breaksObject@changePlots     <- NULL
    breaksObject@changeDistPlots <- NULL
    breaksObject@mapPlots        <- NULL
    object@breaks                <- breaksObject
    return(object)
}

checkBreaks <- function(object) {
    if (! is(object, "msscafClass")) {
        stop("Parameter should be a msscafClass.")
    }
    message("Finding statistics.")
    object@data <- purrr::map(object@data, .checkBreaks, chromosomes = object@chromosomes, sizes = object@sizes)
    return(invisible(object))
}

.computeNCells <- function(object) {
    if (object@parameters@metaSize > 1) {
        object@parameters@breakNCells <- object@parameters@minLinkRange * (object@parameters@minLinkRange - 1) / 4
    }
    else {
        object@parameters@breakNCells <- (object@parameters@minLinkRange - 1) * (object@parameters@minLinkRange - 2) / 4
    }
    return(object)
}

computeNCells <- function(object) {
    if (! is(object, "msscafClass")) {
        stop("Parameter should be a msscafClass.")
    }
    message("Estimating break thresholds")
    object@data <- purrr::map(object@data, .computeNCells)
    return(invisible(object))
}

.computeBreakPvalue <- function(object, pvalue, sizes) {
    if (! is(object, "msscafExp")) {
        stop(paste0("Parameter should be a msscafExp, it is a ", is(object), " ."))
    }
    message(paste0("\tDataset '", object@name , "'."))
    tmp <- object@breaks@data %>%
        dplyr::filter(nCells >= object@parameters@breakNCells) %>%
        dplyr::filter(fcMeanCount >= 0) %>%
        dplyr::pull(fcMeanCount)
    standardDev <- sd(c(tmp, -tmp))
    nTestedPvalues <- object@breaks@data %>%
        dplyr::filter(fcMeanCount < 0) %>%
        dplyr::filter(nCells >= object@parameters@breakNCells) %>%
        nrow()
    object@breaks@data <- object@breaks@data %>%
        dplyr::mutate(pvalue = pnorm(fcMeanCount,
            mean = 0.0, sd = standardDev)) %>%
        dplyr::mutate(pvalue = dplyr::if_else(fcMeanCount >= 0, 1, pvalue)) %>%
        dplyr::mutate(pvalue = dplyr::if_else(nCells < object@parameters@breakNCells, 1, pvalue)) %>%
        dplyr::mutate(testedPvalue = dplyr::if_else(fcMeanCount < 0 &
                                                    nCells >= object@parameters@breakNCells,
                                         pvalue, NA_real_)) %>%
        dplyr::mutate(padj = p.adjust(testedPvalue, method = "BH", n = nTestedPvalues)) %>%
        dplyr::select(-testedPvalue)
    # Fill holes if use meta bins
    if (object@parameters@metaSize > 1) {
        object@breaks@data <- tibble::enframe(sizes, name = "ref", value = "bin") %>%
            dplyr::mutate(ref = factor(ref, levels = names(sizes))) %>%
            dplyr::mutate(bin = purrr::map(bin, ~ seq.int(from = 1, to = .x))) %>%
            tidyr::unchop(bin) %>%
            dplyr::left_join(object@breaks@data %>%
                    dplyr::mutate(padj = dplyr::if_else(is.na(padj), -1, padj)),
                by = c("ref", "bin")) %>%
            tidyr::fill(fcMeanCount, nCells, pvalue, padj) %>%
             dplyr::mutate(padj = dplyr::if_else(padj == -1, NA_real_, padj))
    }
    return(object)
}

computeBreakPvalue <- function(object) {
    if (! is(object, "msscafClass")) {
        stop("Parameter should be a msscafClass.")
    }
    message("Estimating break thresholds")
    object@data <- purrr::map(object@data, .computeBreakPvalue, sizes = object@sizes)
    return(invisible(object))
}

removeNearEqualBreaks <- function(reference, breaks, distance, size, pvalueThreshold) {
    maxNBreaks <- 10000
    breaks     <- breaks %>%
        dplyr::filter(fcMeanCount <  0) %>%
        dplyr::filter(padj        <= pvalueThreshold) %>%
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
    breaks           <- removeNearEqualBreaks(reference, breaksRef, objectRef@parameters@minLinkRange, objectRef@size, pvalueThreshold)
    pb$tick()
    return(list(data       = breaks,
                changePlot = plotTriangles(originalBreakRef),
                mapPlot    = plot.10XRef(objectRef, TRUE, breaks$bin)))
}

.filterBreaks <- function(object, pvalueThreshold, chromosomes, sizes) {
    if (! is(object, "msscafExp")) {
        stop(paste0("Parameter should be a msscafExp, it is a ", is(object), " ."))
    }
    message(paste0("\tDataset '", object@name, "'..."))
    object@breaks@filteredData <- object@breaks@data %>%
                                          dplyr::filter(padj <= pvalueThreshold)
#   object@breaks@filteredData <- removeNearEqualBreaksCpp(object@breaks@data %>%
#                                         dplyr::filter(padj <= pvalueThreshold) %>%
#                                         dplyr::arrange(padj),
#                                     object@parameters@minLinkRange * object@parameters@metaSize) %>%
#       tibble::as_tibble()
    object@breaks@changePlots  <- c()
    object@breaks@mapPlots     <- c()
#   object@breaks@filteredData <- object@breaks@data %>% slice_head(n = 0)
#   object@breaks@changePlots  <- c()
#   object@breaks@mapPlots     <- c()
#   if (nrow(object@breaks@data) == 0) {
#       message("\t\tNo break found.")
#       return(object)
#   }
#   splitBreaks <- object@breaks@data %>% dplyr::group_by(ref)
#   splitNames  <- splitBreaks %>% dplyr::group_keys() %>% dplyr::pull(ref) %>% as.character()
#   splitBreaks <- splitBreaks %>% dplyr::group_split()
#   splitObject <- splitByRef(object, chromosomes, sizes)[splitNames]
#   functionParameters <- transpose(list(ref = splitNames, object = splitObject, breaks = splitBreaks))
#   # parallel seems to need too much RAM
#   #selectedBreaks <- bplapply(functionParameters, filterBreak)
#   pb             <- progress_bar$new(total = length(splitNames))
#   selectedBreaks <- lapply(functionParameters, filterBreak, pvalueThreshold = pvalue, pb = pb)
#   selectedBreaks <- transpose(selectedBreaks)
#   if (length(selectedBreaks$data) == 0) {
#       message("\t\tNo break passed the filter.")
#       return(object)
#   }
#   object@breaks@filteredData <- bind_rows(selectedBreaks$data)
#   object@breaks@changePlots  <- setNames(selectedBreaks$changePlot, splitNames)
#   object@breaks@mapPlots     <- setNames(selectedBreaks$mapPlot, splitNames)
    message(paste0("\t\t... filtered to ", nrow(object@breaks@filteredData), " breaks."))
    return(object)
}

filterBreaks <- function(object, pvalue) {
    if (! is(object, "msscafClass")) {
        stop("Parameter should be a msscafClass.")
    }
    message("Filtering results.")
    object@data <- purrr::map(object@data, .filterBreaks, pvalue = pvalue, chromosomes = object@chromosomes, sizes = object@sizes)
    return(invisible(object))
}

..compareBreaks <- function(object1, object2) {
    if (! is(object1, "msscafExp")) {
        stop("Parameter 1 should be a msscafExp.")
    }
    if (! is(object2, "msscafExp")) {
        stop("Parameter 2 should be a msscafExp.")
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
        dplyr::filter(is.na(fcMeanCount_other) | (fcMeanCount_other <= 0)) %>%
        dplyr::filter(is.na(pvalue_other) | (pvalue_other < 0.5) | (nCells_other <= object2@parameters@breakNCells)) %>%
        #dplyr::filter(nCells < object2@parameters@breakNCells | fcMeanCount <= 0) %>%
        dplyr::select(ref, bin, nCells, fcMeanCount, pvalue, padj)
    object1@breaks@filteredData <- newData
    message(paste0("\tDataset '", object1@name, "', compared to ", object2@name, ", is filtered down to ", nrow(object1@breaks@filteredData), " breaks."))
    return(object1)
}

.compareBreaks <- function(object1, objects) {
    if (! is(object1, "msscafExp")) {
        stop("Parameter 1 should be a msscafExp.")
    }
    object1 <- purrr::reduce(objects, ..compareBreaks, .init = object1)
    return(object1)
}

compareBreaks <- function(object) {
    if (! is(object, "msscafClass")) {
        stop("Parameter should be a msscafClass.")
    }
    message("Comparing breaks")
    object@data <- purrr::map(object@data, .compareBreaks, objects = object@data)
    return(invisible(object))
}

.mergeBreaks <- function(object, pvalueThreshold) {
    object@breaks@filteredData <- object@breaks@filteredData %>%
        dplyr::arrange(padj) %>%
        removeNearEqualBreaksCpp(object@parameters@minLinkRange * object@parameters@metaSize) %>%
        tibble::as_tibble()
}

mergeBreaks <- function(object, pvalueThreshold) {
    if (! is(object, "msscafClass")) {
        stop("Parameter should be a msscafClass.")
    }
    linkRanges      <- unlist(purrr::map(purrr::map(object@data, "parameters"), "minLinkRange"))
    metaSizes       <- unlist(purrr::map(purrr::map(object@data, "parameters"), "metaSize"))
    maxRange        <- max(linkRanges * metaSizes)
    # object@breaks   <- purrr::map_dfr(object@data, .mergeBreaks) %>%
    object@breaks   <- dplyr::bind_rows(purrr::map(purrr::map(object@data, "breaks"), "filteredData")) %>%
        dplyr::select(ref, bin, pvalue)
    object@breaks <- object@breaks %>%
        dplyr::group_by(ref, bin) %>%
        dplyr::summarize(pvalue = exp(sum(log(pvalue))), .groups = "drop") %>%
        dplyr::arrange(pvalue) %>%
        removeNearEqualBreaksCpp(maxRange) %>%
        tibble::as_tibble() %>%
        dplyr::mutate(padj = p.adjust(pvalue, method = "BH", n = nrow(.))) %>%
        dplyr::arrange(ref, bin)
#   breaks          <- dplyr::bind_rows(purrr::map(purrr::map(object@data, "breaks"), "filteredData")) %>% group_by(ref)
#   refs            <- dplyr::group_keys(breaks) %>% tibble::deframe()
#   breaks          <- dplyr::group_split(breaks)
#   sizes           <- object@sizes[refs]
#   nGroups         <- length(refs)
#   minMaxLinkRange <- rep(minMaxLinkRange, nGroups)
#   pvalueThreshold <- rep(pvalueThreshold, nGroups)
#   breaks          <- pmap(list(refs, breaks, minMaxLinkRange, sizes, pvalueThreshold), removeNearEqualBreaks)
#   breaks          <- dplyr::bind_rows(purrr::map(purrr::map(object@data, "breaks"), "filteredData"))
#   object@breaks   <- removeNearEqualBreaksCpp(breaks %>% dplyr::arrange(padj), minMaxLinkRange) %>% tibble::as_tibble()
    # object@breaks   <- breaks %>% bind_rows()
    if (nrow(object@breaks) == 0) {
        object@breaks <- dplyr::slice_head(object@data[[1]]@breaks@data, n = 0)
        message("No break found.")
    }
    else {
        message(paste0("\t", nrow(object@breaks), " breaks found."))
    }
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
    breakPlots        <- purrr::map(object@data, ..addBreakPlots, ref = ref, bin = bin, size = object@sizes[[ref]])
    names(breakPlots) <- purrr::map(object@data, "name")
    pb$tick()
    return(breakPlots)
}

addBreakPlots <- function(object) {
    if (! is(object, "msscafClass")) {
        stop("Parameter should be a msscafClass.")
    }
    message(paste0("Computing plots."))
    pb                       <- progress_bar$new(total = nrow(object@breaks))
    object@breakPlots        <- purrr::map(object@breaks %>% dplyr::mutate(ref = as.character(ref)) %>% transpose(), .addBreakPlots, object = object, pb = pb)
    names(object@breakPlots) <- object@breaks %>% unite("name", c(ref, bin), sep = "_") %>% pull(name)
    return(object)                                                                                                                                                                                      
}

findBreaks <- function(object, pvalue = 0.05) {
    if (! is(object, "msscafClass")) {
        stop("Parameter should be a msscafClass.")
    }
    object <- computeNCells(object)
    object <- checkBreaks(object)
    object <- computeBreakPvalue(object)
    object <- filterBreaks(object, pvalue)
    object <- compareBreaks(object)
    object <- mergeBreaks(object, pvalue)
    # if (nrow(object@breaks) <= 1000) object <- addBreakPlots(object)
    gc(verbose = FALSE)
    return(invisible(object))
}
