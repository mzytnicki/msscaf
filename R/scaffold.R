# Change the UB-LR system to after/collinear
# When ref can be merged to several refs, choose the best one.
selectJoins <- function(joins) {
    # Change to the after/collinear system
    joins1 <- joins %>%
        dplyr::select(ref1, ref2, vert, hor, pvalue) %>%
        dplyr::mutate(after     = (hor == "R")) %>%
        dplyr::mutate(collinear = ((hor == "R") == (vert == "U"))) %>%
        dplyr::select(-c(vert, hor))
    # Create symmetric case (ref1 and 2 inverted)
    joins2 <- joins1 %>%
        dplyr::rename(tmp = ref1) %>%
        dplyr::rename(ref1 = ref2) %>%
        dplyr::rename(ref2 = tmp) %>%
        dplyr::mutate(after = ((!collinear) == after))
    # Choose best "after" ref
    filteredJoins <- dplyr::bind_rows(joins1, joins2) %>%
        dplyr::group_by(ref1, after) %>%
        dplyr::filter(pvalue == min(pvalue)) %>%
        dplyr::ungroup() %>%
        dplyr::filter(ref1 < ref2)
    merges <- dplyr::slice(filteredJoins, 0)
    # Remove possible duplicates
    while (nrow(filteredJoins) > 0) {
        bestLine <- filteredJoins %>% dplyr::arrange(pvalue) %>% dplyr::slice(1) %>% as.list()
        merges <- dplyr::bind_rows(merges, as_tibble(bestLine))
        filteredJoins <- filteredJoins %>%
            dplyr::filter((ref1  != bestLine$ref1) |
                          (ref2  != bestLine$ref2) |
                          (after != bestLine$after))
    }
    merges %<>% dplyr::arrange(pvalue)
    return(merges)
}

# Set the largest ref as reference
orderJoins <- function(object, joins) {
    joins %>%
        dplyr::mutate(size1     = object@sizes[ref1]) %>%
        dplyr::mutate(size2     = object@sizes[ref2]) %>%
        dplyr::mutate(reference = dplyr::if_else(size1 >= size2, ref1, ref2)) %>%
        dplyr::mutate(other     = dplyr::if_else(size1 >= size2, ref2, ref1)) %>%
        dplyr::mutate(after     = ((!collinear) == after)) %>%
        dplyr::select(reference, other, after, collinear) %>%
        dplyr::filter(reference != other)
}

.stitchChromosomePair <- function(object, sizes, reference, other, after, collinear) {
    referenceSize <- sizes[[reference]]
    otherSize     <- sizes[[other]]
    if (!collinear) {
        object@interactionMatrix %<>%
            dplyr::mutate(bin1 = dplyr::if_else(ref1 == other, as.integer(otherSize - bin1 + 1), as.integer(bin1))) %>%
            dplyr::mutate(bin2 = dplyr::if_else(ref2 == other, as.integer(otherSize - bin2 + 1), as.integer(bin2)))
    }
    if (after) {
        object@interactionMatrix %<>%
            dplyr::mutate(bin1 = dplyr::if_else(ref1 == other, as.integer(bin1 + referenceSize), as.integer(bin1))) %>%
            dplyr::mutate(bin2 = dplyr::if_else(ref2 == other, as.integer(bin2 + referenceSize), as.integer(bin2)))
    } else {
        object@interactionMatrix %<>%
            dplyr::mutate(bin1 = dplyr::if_else(ref1 == reference, as.integer(bin1 + otherSize), as.integer(bin1))) %>%
            dplyr::mutate(bin2 = dplyr::if_else(ref2 == reference, as.integer(bin2 + otherSize), as.integer(bin2)))
    }
    referenceFactor <- factor(rep(reference, nrow(object@interactionMatrix)), levels = levels(object@interactionMatrix$ref1))
    object@interactionMatrix %<>%
        dplyr::mutate(ref1 = dplyr::if_else(ref1 == other, referenceFactor, ref1)) %>%
        dplyr::mutate(ref2 = dplyr::if_else(ref2 == other, referenceFactor, ref2))
    return(object)
}

stitchChromosomePair <- function(object, reference, other, after, collinear) {
    object@data <- purrr::map(object@data, .stitchChromosomePair, sizes = object@sizes, reference = reference, other = other, after = after, collinear = collinear)
    otherSeq    <- object@sequences[[other]]
    if (! collinear) {
        otherSeq <- reverseComplement(otherSeq)
    }
    firstSeq    <- if (after) object@sequences[[reference]] else otherSeq
    secondSeq   <- if (after) otherSeq                      else object@sequences[[reference]]
    firstSize   <- if (after) object@sizes[[reference]]     else object@sizes[[other]]
    nMissingNts <- firstSize * object@binSize - length(firstSeq)
    if (nMissingNts > 0) {
        filler   <- DNAString(str_c(rep("N", nMissingNts), collapse = ""))
        firstSeq <- xscat(firstSeq, filler)
    }
    object@sequences[[reference]]   <- xscat(firstSeq, secondSeq)
    object@sizes[[reference]]       <- object@sizes[[reference]] + object@sizes[[other]]
    object@mergedChromosomes[other] <- reference
    object@sizes                    <- object@sizes[names(object@sizes) != other]
    object@chromosomes              <- object@chromosomes[object@chromosomes != other]
    object@sequences                <- object@sequences[names(object@sequences) != other]
    return(object)
}

repairJoins <- function(joins, referenceRef, otherRef, afterJoin, collinearJoin) {
    joins %>%
        dplyr::mutate(after = dplyr::if_else((reference == otherRef),
                              (collinearJoin) == (after),
                              after)) %>%
        dplyr::mutate(collinear = dplyr::if_else((reference == otherRef) | (other == otherRef),
                                  (collinearJoin) == (collinear),
                                  collinear)) %>%
        dplyr::mutate(reference = dplyr::if_else(reference == otherRef,
                                  referenceRef,
                                  reference)) %>%
        dplyr::mutate(other = dplyr::if_else(other == otherRef,
                              referenceRef,
                              other))
}

scaffold <- function(object) {
    selectedJoins <- selectJoins(object@joins)
    orderedJoins  <- orderJoins(object, selectedJoins)
    while (nrow(orderedJoins) != 0) {
        firstRow     <- orderedJoins %>% dplyr::slice(1) %>% as.list()
        orderedJoins <- orderedJoins %>% dplyr::slice(-1)
        reference    <- firstRow$reference
        other        <- firstRow$other
        after        <- firstRow$after
        collinear    <- firstRow$collinear
        message(paste0("  Stitching ", reference, " with ", other, ", ", nrow(orderedJoins), " remaining."))
        object       <- stitchChromosomePair(object, reference, other, after, collinear)
        if (nrow(orderedJoins) != 0) {
            orderedJoins <- repairJoins(orderedJoins, reference, other, after, collinear) %>%
                dplyr::rename(ref1 = reference, ref2 = other)
            orderedJoins <- orderJoins(object, orderedJoins)
        }
    }
    return(object)
}
