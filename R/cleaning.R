removeSmallScaffolds <- function (object) {
    message("Removing small scaffolds.")
    bigRefs <- object@interactionMatrix %>%
        makeSymmetric() %>%
        group_by(ref1) %>%
        summarise(refSizes = max(bin1)) %>%
        ungroup() %>%
        filter(refSizes >= object@minNBins) %>%
        select(ref1) %>%
        pull()
    object@chromosomes <- as.vector(unique(sort(bigRefs)))
    message(paste0("Keeping ", length(object@chromosomes), " big references."))
    object@interactionMatrix <- object@interactionMatrix %>%
        filter(ref1 %in% object@chromosomes, ref2 %in% object@chromosomes)
    return(object)
}

removeLowCountRows <- function (object) {
    message("Removing low counts.")
    tmp <- object@interactionMatrix %>%
        makeSymmetric() %>%
        group_by(ref1, bin1) %>%
        summarise(countSum = sum(count)) %>%
        ungroup()
    object@lowCounts <- tmp %>% filter(countSum < object@minRowCount)
    message(paste0("Removing ", nrow(object@lowCounts), " rows."))
    object@interactionMatrix <- object@interactionMatrix %>%
        anti_join(object@lowCounts, by = c("ref1" = "ref1", "bin1" = "bin1"))
    return(object)
}

cleanData <- function (object) {
    object <- removeSmallScaffolds(object)
    object <- removeLowCountRows(object)
    return(object)
}