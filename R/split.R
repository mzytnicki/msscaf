splitChromosome <- function(interactionMatrix, sizes, splitPoint, oldChromosome, newChromosome) {
    previousSize <- sizes[[oldChromosome]]
    newSize      <- previousSize - splitPoint
    firstPart    <- (splitPoint > previousSize / 2)
    shiftedChr   <- if (firstPart) newChromosome else oldChromosome
    if (firstPart) {
        comparator <- function(x) { (x < splitPoint) }
    } else {
        comparator <- function(x) { (x > splitPoint) }
    }
    interactionMatrix <- interactionMatrix %>% 
        mutate(ref1 = as.character(ref1)) %>%
        mutate(ref2 = as.character(ref2)) %>%
        mutate(ref1 = ifelse((ref1 == oldChromosome) & (!comparator(bin1)), newChromosome, ref1)) %>%
        mutate(ref2 = ifelse((ref2 == oldChromosome) & (!comparator(bin2)), newChromosome, ref2)) %>%
        mutate(bin1 = ifelse(ref1 == shiftedChr, bin1 - splitPoint + 1, bin1)) %>%
        mutate(bin2 = ifelse(ref2 == shiftedChr, bin2 - splitPoint + 1, bin2))
    if (firstPart) {
        sizes[[oldChromosome]] <- splitPoint
        sizes[[newChromosome]] <- newSize
    }
    else {
        sizes[[newChromosome]] <- splitPoint
        sizes[[oldChromosome]] <- newSize
    }
    return(list(interactionMatrix = interactionMatrix, sizes = sizes))
}


splitChromosomes <- function(object, splitPoints) {
    interactionMatrix <- object@interactionMatrix
    sizes             <- object@sizes
    newChromosomes    <- paste0("new_chr_", seq.int(nrow(splitPoints)))
    splitPoints       <- splitPoints %>% mutate(newChromosome = newChromosomes)
    nSplits           <- nrow(splitPoints)
    for (nSplit in seq.int(to = nSplits)) {
        message(paste0("  Splitting point ", nSplit, "/", nSplits, ": ", splitPoints$ref[[nSplit]]))
        newChromosome <- paste0("new_chr_", nSplit)
        out <- splitChromosome(interactionMatrix,
                               sizes,
                               splitPoints$bin[[nSplit]],
                               splitPoints$ref[[nSplit]],
                               splitPoints$newChromosome[[nSplit]])
        interactionMatrix <- out$interactionMatrix
        sizes             <- out$sizes
    }
    object@sizes                 <- sizes
    object@interactionMatrix     <- interactionMatrix
    object@chromosomes           <- c(object@chromosomes, newChromosomes)
    object@newChromosomes        <- splitPoints$ref
    names(object@newChromosomes) <- newChromosomes
    return(object)
}