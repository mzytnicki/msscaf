.splitChromosome <- function(object, refs, prevRef, newRef, shiftedRef, splitPoint, comparator) {
    if (! is(object, "tenxcheckerExp")) {
        stop("First parameter should be a tenxcheckerExp.")
    }
    #message(str(prevRef))
    #message(str(newRef))
    #message(str(object@interactionMatrix))
    newRef <- factor(newRef, levels = refs)
    levels(object@interactionMatrix$ref1) <- refs
    levels(object@interactionMatrix$ref2) <- refs
    object@interactionMatrix <- object@interactionMatrix %>% 
        dplyr::mutate(ref1 = if_else((ref1 == prevRef) & (!comparator(bin1)), newRef, ref1)) %>%
        dplyr::mutate(ref2 = if_else((ref2 == prevRef) & (!comparator(bin2)), newRef, ref2)) %>%
        dplyr::mutate(bin1 = if_else(ref1 == shiftedRef, as.integer(bin1 - splitPoint + 1), bin1)) %>%
        dplyr::mutate(bin2 = if_else(ref2 == shiftedRef, as.integer(bin2 - splitPoint + 1), bin2))
    return(object)
}

splitChromosome <- function(object, parameters) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    prevRef    <- parameters$ref
    splitPoint <- parameters$bin
    newRef     <- parameters$newRef
    #message(str(parameters))
    #message(str(prevRef))
    prevSize   <- object@sizes[[prevRef]]
    if (splitPoint >= prevSize) {
        stop(paste0("Error while splitting chromosome '", prevRef, "' of size ", prevSize, " at point ", splitPoint, "."))
    }
    newSize    <- prevSize - splitPoint
    firstPart  <- (splitPoint > prevSize / 2)
    shiftedRef <- if (firstPart) newRef else prevRef
    if (firstPart) {
        comparator <- function(x) { (x < splitPoint) }
    } else {
        comparator <- function(x) { (x > splitPoint) }
    }
    object@data <- map(object@data, .splitChromosome, object@chromosomes, prevRef = prevRef, newRef = newRef, shiftedRef = shiftedRef, splitPoint = splitPoint, comparator = comparator)
    if (firstPart) {
        object@sizes[[prevRef]] <- splitPoint
        object@sizes[[newRef]]  <- newSize
    }
    else {
        object@sizes[[newRef]]  <- splitPoint
        object@sizes[[prevRef]] <- newSize
    }
    return(object)
}

splitChromosomes <- function(object) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    nSplits            <- nrow(object@breaks)
    newRefs            <- paste0("new_ref_", seq.int(nrow(object@breaks)))
    object@chromosomes <- c(object@chromosomes, newRefs)
    parameters <- object@breaks %>%
        arrange(ref, desc(bin)) %>%
        mutate(newRef = newRefs) %>%
        transpose()
    pb <- progress_bar$new(total = nSplits)
    for (param in parameters) {
        object <- splitChromosome(object, param)
        pb$tick()
    }
    return(invisible(object))
}
