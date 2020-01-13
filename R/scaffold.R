selectJoins <- function(joins) {
    joins1 <- joins %>%
        select(ref1, ref2, vert, hor, pvalue) %>%
        mutate(after     = (hor == "R")) %>%
        mutate(collinear = ((hor == "R") == (vert == "U"))) %>%
        select(-c(vert, hor))
    joins2 <- joins1 %>%
        rename(tmp = ref1) %>%
        rename(ref1 = ref2) %>%
        rename(ref2 = tmp) %>%
        mutate(after = ((!collinear) == after))
    filteredJoins <- bind_rows(joins1, joins2) %>%
        group_by(ref1, after) %>%
        filter(pvalue == min(pvalue)) %>%
        ungroup() %>%
        filter(ref1 < ref2)
    merges <- slice(filteredJoins, 0)
    while (nrow(filteredJoins) > 0) {
        bestLine <- filteredJoins %>% arrange(pvalue) %>% head(1) %>% as.list()
        merges <- bind_rows(merges, as_tibble(bestLine))
        filteredJoins <- filteredJoins %>%
            filter((ref1  != bestLine$ref1) |
                   (ref2  != bestLine$ref2) |
                   (after != bestLine$after))
    }
    merges %<>% arrange(pvalue)
    return(merges)
}

orderJoins <- function(object, joins) {
    joins %>%
        mutate(size1 = object@sizes[ref1]) %>%
        mutate(size2 = object@sizes[ref2]) %>%
        mutate(reference = ifelse(size1 >= size2, ref1, ref2)) %>%
        mutate(other = ifelse(size1 >= size2, ref2, ref1)) %>%
        mutate(after = ((!collinear) == after)) %>%
        select(reference, other, after, collinear) %>%
        filter(reference != other)
}

stitchChromosomePair <- function(object, reference, other, after, collinear) {
    referenceSize <- object@sizes[[reference]]
    otherSize     <- object@sizes[[other]]
    if (!collinear) {
        object@interactionMatrix %<>%
            mutate(bin1 = ifelse(ref1 == other, otherSize - bin1 + 1, bin1)) %>%
            mutate(bin2 = ifelse(ref2 == other, otherSize - bin2 + 1, bin2))
    }
    if (after) {
        object@interactionMatrix %<>%
            mutate(bin1 = ifelse(ref1 == other, bin1 + referenceSize, bin1)) %>%
            mutate(bin2 = ifelse(ref2 == other, bin2 + referenceSize, bin2))
    } else {
        object@interactionMatrix %<>%
            mutate(bin1 = ifelse(ref1 == reference, bin1 + otherSize, bin1)) %>%
            mutate(bin2 = ifelse(ref2 == reference, bin2 + otherSize, bin2))
    }
    object@interactionMatrix %<>%
        mutate(ref1 = ifelse(ref1 == other, reference, as.character(ref1))) %>%
        mutate(ref2 = ifelse(ref2 == other, reference, as.character(ref2)))
    object@sizes <- object@sizes[names(object@sizes) != other]
    object@chromosomes <- object@chromosomes[object@chromosomes != other]
    object@mergedChromosomes[other] <- reference
    return(object)
}

repairJoins <- function(joins, referenceRef, otherRef, afterJoin, collinearJoin) {
    joins %<>%
        mutate(after = ifelse((reference == otherRef),
                              (collinearJoin) == (after),
                              after)) %>%
        mutate(collinear = ifelse((reference == otherRef) | (other == otherRef),
                                  (collinearJoin) == (collinear),
                                  collinear)) %>%
        mutate(reference = ifelse(reference == otherRef,
                                  referenceRef,
                                  reference)) %>%
        mutate(other = ifelse(other == otherRef,
                              referenceRef,
                              other))
}

scaffold <- function(object, joins) {
    selectedJoins <- selectJoins(joins)
    orderedJoins <- orderJoins(object, selectedJoins)
    while (nrow(orderedJoins) != 0) {
        firstRow     <- orderedJoins %>% slice(1) %>% as.list()
        orderedJoins <- orderedJoins %>% slice(-1)
        reference    <- firstRow$reference
        other        <- firstRow$other
        after        <- firstRow$after
        collinear    <- firstRow$collinear
        message(paste0("  Stitching ", reference, " with ", other, ", ", nrow(orderedJoins), " remaining."))
        object       <- stitchChromosomePair(object, reference, other, after, collinear)
        if (nrow(orderedJoins) != 0) {
            orderedJoins <- repairJoins(orderedJoins, reference, other, after, collinear) %>%
                rename(ref1 = reference, ref2 = other)
            orderedJoins <- orderJoins(object, orderedJoins)
        }
    }
    return(object)
}
