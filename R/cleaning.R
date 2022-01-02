.removeSmallScaffoldsOutlierBins <- function(object, refs) {
    if (! is(object, "tenxcheckerExp")) {
        stop("Parameter should be a tenxcheckerExp.")
    }
    object@outlierBins <- object@outlierBins %>%
        dplyr::mutate(ref = as.character(ref)) %>%
        dplyr::filter(ref %in% refs) %>%
        dplyr::mutate(ref = factor(ref, levels = refs))
    return(object)
}

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
    object@data <- purrr::map(object@data, .removeSmallScaffoldsOutlierBins, refs = seenRefs)
    message(paste0("\tKeeping ", length(seenRefs), " scaffolds with at least ", object@minNBins, " bins."))
    return(object)
}

.removeLowCountRows <- function(object, sizes) {
    if (! is(object, "tenxcheckerExp")) {
        stop("Parameter should be a tenxcheckerExp.")
    }
    message(paste0("\tDataset '", object@name , "'."))
    output <- removeLowCountRowsCpp(object@interactionMatrix, sizes, object@parameters@minRowCount)
    object@interactionMatrix <- as_tibble(output$matrix)
    object@lowCountRows      <- as_tibble(output$rows)
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
    #normalizeHighCountRows(object)
    return(invisible(object))
}
