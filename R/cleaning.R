removeSmallScaffolds <- function(object) {
    message(paste0("Removing small scaffolds (currently: ", object@interactionMatrix %>% select(ref1) %>% distinct() %>% nrow() , ")."))
    bigRefs <- object@interactionMatrix %>%
        filter(ref1 == ref2) %>%
        makeSymmetric() %>%
        group_by(ref1) %>%
        summarise(refSizes = max(bin1)) %>%
        ungroup() %>%
        filter(refSizes >= object@parameters@minNBins) %>%
        select(ref1) %>%
        pull()
    object@chromosomes <- as.vector(unique(sort(bigRefs)))
    message(paste0("Keeping ", length(object@chromosomes), " big references."))
    object@interactionMatrix %<>%
        filter(ref1 %in% object@chromosomes, ref2 %in% object@chromosomes) %>%
        mutate(ref1 = droplevels(ref1)) %>%
        mutate(ref1 = factor(ref1, levels = object@chromosomes)) %>%
        mutate(ref2 = droplevels(ref2)) %>%
        mutate(ref2 = factor(ref2, levels = object@chromosomes))
    object@sizes <- object@sizes[object@chromosomes]
    return(object)
}

removeLowCountRows <- function(object) {
    message("Removing low counts rows.")
    object@lowCounts <- object@interactionMatrix %>%
        makeSymmetric() %>%
        group_by(ref1, bin1) %>%
        summarise(countSum = sum(count)) %>%
        ungroup() %>%
        filter(countSum < object@parameters@minRowCount)
    message(paste0("Removing ", nrow(object@lowCounts), " rows."))
    object@interactionMatrix <- object@interactionMatrix %>%
        anti_join(object@lowCounts, by = c("ref1" = "ref1", "bin1" = "bin1"))
    return(object)
}

removeLowCount <- function(object) {
    message("Removing low counts.")
    nLowCounts <- object@interactionMatrix %>%
        filter(count < object@parameters@minCount) %>%
        nrow()
    object@interactionMatrix %<>%
        filter(count >= object@parameters@minCount)
    message(paste0("Removing ", nLowCounts, " cells."))
    return(object)
}

removeFarFromDiagonal <- function(object) {
    message("    Removing values far from diagonal.")
    object@interactionMatrix %<>%
        filter(abs(bin1 - bin2) <= object@parameters@maxLinkRange)
    return(object)
}

cleanData <- function(object) {
    object <- removeLowCountRows(object)
    object <- removeSmallScaffolds(object)
    object <- removeLowCount(object)
    return(object)
}