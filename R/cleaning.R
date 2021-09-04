removeSmallScaffolds <- function(object) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    message(paste0("Removing small scaffolds (currently: ", length(object@chromosomes), ")."))
    chromosomes <- names(object@sizes[object@sizes >= object@minNBins])
    object      <- keepScaffolds(object, chromosomes)
    seenRefs    <- mixedsort(unique(unlist(map(object@data, ~ unique(c(as.character(unique(.x@interactionMatrix$ref1)),
                                                                       as.character(unique(.x@interactionMatrix$ref2))))))))
    object      <- keepScaffolds(object, seenRefs)
    message(paste0("\tKeeping ", length(seenRefs), " scaffolds with at least ", object@minNBins, " bins."))
    return(object)
}

.removeLowCountRows <- function(object, sizes) {
    if (! is(object, "tenxcheckerExp")) {
        stop("Parameter should be a tenxcheckerExp.")
    }
    message(paste0("\tDataset '", object@name , "'."))
    object@interactionMatrix <- as_tibble(removeLowCountRowsCpp(object@interactionMatrix, sizes, object@parameters@minRowCount))
    return(object)
}

removeLowCountRows <- function(object) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    message(paste0("Removing low counts rows."))
    object@data <- map(object@data, .removeLowCountRows, sizes = object@sizes)
    return(object)
}

.removeLowCount <- function(object) {
    if (! is(object, "tenxcheckerExp")) {
        stop("Parameter should be a tenxcheckerExp.")
    }
    nCells <- object@interactionMatrix %>%
        nrow()
    nLowCounts <- object@interactionMatrix %>%
        filter(count < object@parameters@minCount) %>%
        nrow()
    object@interactionMatrix %<>%
        filter(count >= object@parameters@minCount)
    message(paste0("\tDataset '", object@name , "': removed ", nLowCounts, " / ", nCells, " cells."))
    return(object)
}

removeLowCount <- function(object) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    message("Removing low counts.")
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
    object <- removeSmallScaffolds(object)
    normalizeHighCountRows(object)
    return(invisible(object))
}
