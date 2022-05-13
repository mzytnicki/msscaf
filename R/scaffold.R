# Change the UB-LR system to after/collinear
# When ref can be merged to several refs, choose the best one.
selectJoins <- function(object, joins) {
    # Change to the after/collinear system
    #joins1 <- joins %>%
    joins %>%
        dplyr::select(ref1, ref2, after1, after2, pvalue) %>%
        dplyr::group_by(ref1, ref2, after1, after2) %>%
        dplyr::slice_min(order_by = pvalue, n = 1, with_ties = FALSE) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(pvalue) %>%
        dplyr::mutate(ref1 = factor(ref1, levels = object@chromosomes)) %>%
        dplyr::mutate(ref2 = factor(ref2, levels = object@chromosomes))
#   # Create symmetric case (ref1 and 2 inverted)
#   joins2 <- joins1 %>%
#       dplyr::rename(tmp = ref1) %>%
#       dplyr::rename(ref1 = ref2) %>%
#       dplyr::rename(ref2 = tmp) %>%
#       dplyr::mutate(after = ((!collinear) == after))
#   # Choose best "after" ref
#   filteredJoins <- dplyr::bind_rows(joins1, joins2) %>%
#       dplyr::group_by(ref1, after) %>%
#       dplyr::filter(pvalue == min(pvalue)) %>%
#       dplyr::ungroup() %>%
#       dplyr::filter(ref1 < ref2)
#   merges <- dplyr::slice(filteredJoins, 0)
#   # Remove possible duplicates
#   while (nrow(filteredJoins) > 0) {
#       bestLine <- filteredJoins %>% dplyr::arrange(pvalue) %>% dplyr::slice(1) %>% as.list()
#       merges <- dplyr::bind_rows(merges, tibble::as_tibble(bestLine))
#       filteredJoins <- filteredJoins %>%
#           dplyr::filter((ref1  != bestLine$ref1) |
#                         (ref2  != bestLine$ref2) |
#                         (after != bestLine$after))
#   }
#   merges %<>% dplyr::arrange(pvalue)
#   return(merges)
}

# Set the largest ref as reference
orderJoins <- function(object, joins) {
    joins %>%
        dplyr::mutate(size1     = object@sizes[ref1]) %>%
        dplyr::mutate(size2     = object@sizes[ref2]) %>%
        dplyr::mutate(reference = dplyr::if_else(size1 >= size2, ref1, ref2)) %>%
        dplyr::mutate(other     = dplyr::if_else(size1 >= size2, ref2, ref1)) %>%
        dplyr::mutate(after     = ((!collinear) == after)) %>%
        dplyr::mutate(refSize   = object@sizes[reference]) %>%
        dplyr::mutate(otherSize = object@sizes[other]) %>%
        dplyr::select(reference, other, after, collinear, refSize, otherSize) %>%
        dplyr::filter(reference != other)
}

# Update scaffold places by adding one scaffold join
.addScaffold <- function(scaffoldPlaces, join) {
# message(str(join))
# scaffoldPlaces %>% head() %>% str() %>% message()
    forwardDirection <- scaffoldPlaces %>%
        dplyr::filter(ref == join$reference) %>%
        dplyr::slice(1) %>%
        dplyr::pull(forward)
    if (join$after & join$collinear & forwardDirection) {
        scaffoldPlaces <- scaffoldPlaces %>%
            dplyr::mutate(offset  = dplyr::if_else(newRef == join$other, as.integer(join$refSize + offset), offset)) %>%
            dplyr::mutate(newSize = dplyr::if_else(newRef == join$reference, newSize + join$otherSize, newSize)) %>%
            dplyr::mutate(forward = dplyr::if_else(newRef == join$other, forward == join$collinear, forward)) %>%
            dplyr::mutate(newRef  = dplyr::if_else(newRef == join$other, join$reference, newRef))
        return(scaffoldPlaces)
    }
    if (join$after & join$collinear) {
        scaffoldPlaces <- scaffoldPlaces %>%
            dplyr::mutate(offset  = dplyr::if_else(newRef == join$other, as.integer(join$refSize + offset), offset)) %>%
            dplyr::mutate(newSize = dplyr::if_else(newRef == join$reference, newSize + join$otherSize, newSize)) %>%
            dplyr::mutate(forward = dplyr::if_else(newRef == join$other, forward == join$collinear, forward)) %>%
            dplyr::mutate(newRef  = dplyr::if_else(newRef == join$other, join$reference, newRef))
        return(scaffoldPlaces)
    }
    if (join$after) {
        scaffoldPlaces <- scaffoldPlaces %>%
            dplyr::mutate(offset  = dplyr::if_else(newRef == join$other, as.integer(join$refSize + newSize - offset - size + 2), offset)) %>%
            dplyr::mutate(forward = dplyr::if_else(newRef == join$other, ! forward, forward)) %>%
            dplyr::mutate(newSize = dplyr::if_else(newRef == join$reference, newSize + join$otherSize, newSize)) %>%
            dplyr::mutate(newRef  = dplyr::if_else(newRef == join$other, join$reference, newRef))
        return(scaffoldPlaces)
    }
    if ((!join$after) & join$collinear) {
        scaffoldPlaces <- scaffoldPlaces %>%
            dplyr::mutate(offset  = dplyr::if_else(newRef == join$reference, as.integer(join$otherSize + offset), offset)) %>%
            dplyr::mutate(newSize = dplyr::if_else(newRef == join$reference, newSize + join$otherSize, newSize)) %>%
            dplyr::mutate(newRef  = dplyr::if_else(newRef == join$other, join$reference, newRef))
        return(scaffoldPlaces)
    }
    scaffoldPlaces %>%
        dplyr::mutate(offset  = dplyr::if_else(newRef == join$other, as.integer(newSize - offset - size + 2), offset)) %>%
        dplyr::mutate(offset  = dplyr::if_else(newRef == join$reference, as.integer(join$otherSize + offset), offset)) %>%
        dplyr::mutate(newSize = dplyr::if_else(newRef == join$reference, newSize + join$otherSize, newSize)) %>%
        dplyr::mutate(newRef  = dplyr::if_else(newRef == join$other, join$reference, newRef))
}

# Create a table, with the list of new places for the contigs
createScaffoldPlaces <- function(object, joins) {
    scaffoldPlaces <- object@sizes %>%
        tibble::enframe(name = "ref", value = "size") %>%
        dplyr::mutate(newRef = ref) %>%
        dplyr::mutate(offset = as.integer(0)) %>%
        dplyr::mutate(newSize = size) %>%
        dplyr::mutate(forward = TRUE)
    #scaffoldPlaces <- purrr::reduce2(purrr::transpose(joins), .addScaffold, .init = scaffoldPlaces)
    for (join in purrr::transpose(joins)) scaffoldPlaces <- .addScaffold(scaffoldPlaces, join)
    for (join in purrr::transpose(joins)) {
        scaffoldPlaces <- .addScaffold(scaffoldPlaces, join)
    }
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
        filler   <- DNAString(stringr::str_c(rep("N", nMissingNts), collapse = ""))
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

.scaffold <- function(object, orders, groupNames, sizes) {
    message(paste0("\tDataset '", object@name, "'."))
    l <- scaffoldCounts(object@interactionMatrix, object@outlierBins, orders, groupNames, sizes, object@parameters@metaSize)
    object@interactionMatrix <- l$interactionMatrix %>% tibble::as_tibble()
    object@outlierBins       <- l$outlierBins %>% tibble::as_tibble()
    return(object)
}

scaffold <- function(object) {
    if (nrow(object@joins) == 0) {
        message("No split found.")
        return(object)
    }
    selectedJoins    <- selectJoins(object, object@joins)
    #selectedJoins   <- selectJoins(object@joins)
    # orderedJoins   <- orderJoins(object, selectedJoins)
    groups           <- getRefOrders(selectedJoins, object@sizes)
    groupNames       <- as.numeric(factor(names(groups), levels = object@chromosomes))
    orderGroups      <- order(groupNames)
    groups           <- groups[orderGroups]
    groupNames       <- groupNames[orderGroups]
    message("Scaffolding sequences.")
    object@sequences   <- scaffoldContigs(object@sequences, groups, object@sizes, object@binSize)
    message("Scaffolding counts.")
    object@data        <- purrr::map(object@data, .scaffold, groups, groupNames, object@sizes)
    object@sizes       <- scaffoldSizes(groups, groupNames, object@sizes)
    object@chromosomes <- names(object@sizes)
    if (length(object@chromosomes) != length(object@sizes)) {
        stop(paste0("Size of objects differ after scaffolding: ", length(object@chromosomes), " vs ", length(object@sizes), "."))
    }
    checkSizeDifference(object)

#   pb <- progress_bar$new(total = nrow(orderedJoins))
#   while (nrow(orderedJoins) != 0) {
#       firstRow     <- orderedJoins %>% dplyr::slice(1) %>% as.list()
#       orderedJoins <- orderedJoins %>% dplyr::slice(-1)
#       reference    <- firstRow$reference
#       other        <- firstRow$other
#       after        <- firstRow$after
#       collinear    <- firstRow$collinear
#       # message(paste0("  Stitching ", reference, " with ", other, ", ", nrow(orderedJoins), " remaining."))
#       object       <- stitchChromosomePair(object, reference, other, after, collinear)
#       if (nrow(orderedJoins) != 0) {
#           orderedJoins <- repairJoins(orderedJoins, reference, other, after, collinear) %>%
#               dplyr::rename(ref1 = reference, ref2 = other)
#           orderedJoins <- orderJoins(object, orderedJoins)
#       }
#       pb$tick()
#   }
    gc(verbose = FALSE)
    return(object)
}
