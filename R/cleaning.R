removeSmallScaffolds <- function(object) {
    message("Removing small scaffolds.")
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
    return(object)
}

removeLowCountRows <- function(object) {
    message("Removing low counts.")
    tmp <- object@interactionMatrix %>%
        makeSymmetric() %>%
        group_by(ref1, bin1) %>%
        summarise(countSum = sum(count)) %>%
        ungroup()
    object@lowCounts <- tmp %>%
        filter(countSum < object@parameters@minRowCount)
    message(paste0("Removing ", nrow(object@lowCounts), " rows."))
    object@interactionMatrix <- object@interactionMatrix %>%
        anti_join(object@lowCounts, by = c("ref1" = "ref1", "bin1" = "bin1"))
    return(object)
}

removeFarFromDiagonal <- function(object) {
    message("    Removing values far from diagonal.")
    object@interactionMatrix %<>%
        filter(abs(bin1 - bin2) <= object@parameters@maxLinkRange)
    return(object)
}


cleanData <- function(object) {
    object <- removeSmallScaffolds(object)
    object <- removeLowCountRows(object)
    return(object)
}