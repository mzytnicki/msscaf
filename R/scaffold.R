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

.stitchChromosomePair <- function(object, sizes, reference, other, after, collinear) {
    referenceSize <- sizes[[reference]]
    otherSize     <- sizes[[other]]
    if (!collinear) {
        object@interactionMatrix %<>%
            mutate(bin1 = if_else(ref1 == other, as.integer(otherSize - bin1 + 1), as.integer(bin1))) %>%
            mutate(bin2 = if_else(ref2 == other, as.integer(otherSize - bin2 + 1), as.integer(bin2)))
    }
    if (after) {
        object@interactionMatrix %<>%
            mutate(bin1 = if_else(ref1 == other, as.integer(bin1 + referenceSize), as.integer(bin1))) %>%
            mutate(bin2 = if_else(ref2 == other, as.integer(bin2 + referenceSize), as.integer(bin2)))
    } else {
        object@interactionMatrix %<>%
            mutate(bin1 = if_else(ref1 == reference, as.integer(bin1 + otherSize), as.integer(bin1))) %>%
            mutate(bin2 = if_else(ref2 == reference, as.integer(bin2 + otherSize), as.integer(bin2)))
    }
    referenceFactor <- factor(rep(reference, nrow(object@interactionMatrix)), levels = levels(object@interactionMatrix$ref1))
    object@interactionMatrix %<>%
        mutate(ref1 = if_else(ref1 == other, referenceFactor, ref1)) %>%
        mutate(ref2 = if_else(ref2 == other, referenceFactor, ref2))
    return(object)
}

stitchChromosomePair <- function(object, reference, other, after, collinear) {
    object@data                     <- map(object@data, .stitchChromosomePair, sizes = object@sizes, reference = reference, other = other, after = after, collinear = collinear)
    object@sizes[[reference]]       <- object@sizes[[reference]] + object@sizes[[other]]
    object@sizes                    <- object@sizes[names(object@sizes) != other]
    object@chromosomes              <- object@chromosomes[object@chromosomes != other]
    object@mergedChromosomes[other] <- reference
    return(object)
}

repairJoins <- function(joins, referenceRef, otherRef, afterJoin, collinearJoin) {
    joins %<>%
        mutate(after = if_else((reference == otherRef),
                              (collinearJoin) == (after),
                              after)) %>%
        mutate(collinear = if_else((reference == otherRef) | (other == otherRef),
                                  (collinearJoin) == (collinear),
                                  collinear)) %>%
        mutate(reference = if_else(reference == otherRef,
                                  referenceRef,
                                  reference)) %>%
        mutate(other = if_else(other == otherRef,
                              referenceRef,
                              other))
}

scaffold <- function(object) {
    selectedJoins <- selectJoins(object@joins)
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
