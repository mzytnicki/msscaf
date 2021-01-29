.removeSmallScaffolds <- function(object, keptRefs) {
    if (! is(object, "tenxcheckerExp")) {
        stop("Parameter should be a tenxcheckerExp.")
    }
    object@interactionMatrix <- as_tibble(removeSmallScaffoldsCpp(object@interactionMatrix, keptRefs))
    return(object)
}

removeSmallScaffolds <- function(object) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    message(paste0("Removing small scaffolds (currently: ", length(object@chromosomes), ")."))
    object@chromosomes       <- names(object@sizes[object@sizes >= object@minNBins])
    object@sizes             <- object@sizes[object@chromosomes]
    object@data              <- map(object@data, .removeSmallScaffolds, keptRefs = object@chromosomes)
    message(paste0("\tKeeping ", length(object@chromosomes), " scaffolds with at least ", object@minNBins, " bins."))
    return(object)
}

.removeLowCountRows <- function(object, sizes) {
    if (! is(object, "tenxcheckerExp")) {
        stop("Parameter should be a tenxcheckerExp.")
    }
    message(paste0("Dataset '", object@name , "': Removing low counts rows."))
    object@interactionMatrix <- as_tibble(removeLowCountRowsCpp(object@interactionMatrix, sizes, object@parameters@minRowCount))
    return(object)
}

removeLowCountRows <- function(object) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    object@data <- map(object@data, .removeLowCountRows, sizes = object@sizes)
    return(object)
}

.removeLowCount <- function(object) {
    if (! is(object, "tenxcheckerExp")) {
        stop("Parameter should be a tenxcheckerExp.")
    }
    message(paste0("Dataset '", object@name , "': Removing low counts."))
    nCells <- object@interactionMatrix %>%
        nrow()
    nLowCounts <- object@interactionMatrix %>%
        filter(count < object@parameters@minCount) %>%
        nrow()
    object@interactionMatrix %<>%
        filter(count >= object@parameters@minCount)
    message(paste0("\tRemoved ", nLowCounts, " / ", nCells, " cells."))
    return(object)
}

removeLowCount <- function(object) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    object@data <- map(object@data, .removeLowCount)
    return(object)
}

removeFarFromDiagonal <- function(object) {
    #message("    Removing values far from diagonal.")
    object@interactionMatrix %<>%
        filter(abs(bin1 - bin2) <= object@parameters@maxLinkRange)
    return(object)
}

cleanData <- function(object) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    object <- removeLowCount(object)
    object <- removeLowCountRows(object)
    object <- removeSmallScaffolds(object)
    return(invisible(object))
}
