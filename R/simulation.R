cleanMatrices <- function(interactionMatrices, refs) {
    interactionMatrices %>%
        dplyr::filter(count > 0) %>%
        dplyr::filter(ref1 != ref2 | bin1 >= bin2) %>%
        dplyr::mutate(ref1 = factor(ref1, levels = refs)) %>%
        dplyr::mutate(ref2 = factor(ref2, levels = refs))
}

generateCorner <- function(offset, r1, r2, refSize, maxNReads, diagSize, distance) {
    lambda <- as.integer(round(maxNReads * dgeom(distance + offset, prob = 1 / diagSize)))
    tibble::tibble(bin1 = seq.int(from = refSize - offset - 1, to = refSize - 1),
           bin2 = seq.int(from = 0, to = offset),
           count = as.integer(round(rpois(offset + 1, lambda)))) %>%
        dplyr::select(bin1, bin2, count) %>%
        dplyr::mutate(ref1 = r1) %>%
        dplyr::mutate(ref2 = r2)
}

addJoin <- function(interactionMatrices, refSize, maxNReads, ref1, ref2, distance, diagSize) {
    corner <- purrr::map_dfr(seq.int(refSize) - 1, generateCorner, ref1, ref2, refSize, maxNReads, diagSize, distance)
    dplyr::bind_rows(interactionMatrices, corner) %>%
        dplyr::group_by(ref1, bin1, ref2, bin2) %>%
        dplyr::summarise(count = sum(count), .groups = "drop")
}

addJoins <- function(interactionMatrices, refs, refSize, maxNReads, diagSize, ref1, ref2, distances) {
    for (i in seq_along(ref1)) {
        interactionMatrices <- addJoin(interactionMatrices, refSize, maxNReads, refs[[ref1[[i]]]], refs[[ref2[[i]]]], distances[[i]], diagSize)
    }
    return(interactionMatrices)
}

addSplit <- function(interactionMatrices, ref, position) {
    interactionMatrices %>%
        dplyr::mutate(count = dplyr::if_else(ref1 == ref & ref2 == ref & bin1 > position & bin2 < position, 0L, count))
}

addSplits <- function(interactionMatrices, refs, refIds, positions) {
    for (i in seq_along(refIds)) {
        interactionMatrices <- addSplit(interactionMatrices, refs[[refIds[[i]]]], positions[[i]])
    }
    return(interactionMatrices)
}

addNoise <- function(interactionMatrices, gaussianNoise) {
    emptyMatrices <- interactionMatrices %>%
        dplyr::select(ref1, bin1, ref2, bin2) %>%
        dplyr::mutate(count = as.integer(round(abs(rnorm(nrow(.), mean = 0, sd = gaussianNoise))))) %>%
        dplyr::filter(count > 0)
    dplyr::bind_rows(interactionMatrices, emptyMatrices) %>%
        dplyr::group_by(ref1, bin1, ref2, bin2) %>%
        dplyr::summarise(count = sum(count), .groups = "drop")
}

generateDiag <- function(distance, refs, refSize, maxNReads, diagSize) {
    nRefs  <- length(refs)
    lambda <- as.integer(round(maxNReads * dgeom(distance, prob = 1 / diagSize)))
    tibble::tibble(bin1 = as.integer(rep(seq.int(refSize) - 1, nRefs)),
           count = as.integer(round(rpois(refSize * nRefs, lambda)))) %>%
        dplyr::mutate(ref1 = rep(refs, each = refSize)) %>%
        dplyr::filter(count > 0) %>%
        dplyr::mutate(bin2 = as.integer(bin1 - distance)) %>%
        dplyr::filter(bin2 >= 0) %>%
        dplyr::select(ref1, bin1, bin2, count) %>%
        dplyr::mutate(ref2 = ref1, .after = 1)
}

addDiagonalMatrices <- function(interactionMatrices, refs, refSize, maxNReads, diagSize, gaussianNoise) {
    diags <- purrr::map_dfr(seq.int(refSize) - 1, generateDiag, refs, refSize, maxNReads, diagSize)
    dplyr::bind_rows(interactionMatrices, diags) %>%
        dplyr::group_by(ref1, bin1, ref2, bin2) %>%
        dplyr::summarise(count = as.integer(sum(count)), .groups = "drop")
}

makeEmptyInteractionMatrices <- function(refs, refSize) {
    bins    <- as.integer(seq.int(refSize) - 1)
    refBins <- purrr::cross2(refs, bins)
    tibble::tibble(ref1 = refBins %>% purrr::map(1) %>% purrr::flatten_chr(),
           bin1 = refBins %>% purrr::map(2) %>% purrr::flatten_int()) %>%
         dplyr::mutate(ref2 = ref1) %>%
         dplyr::mutate(bin2 = bin1) %>%
         tidyr::expand(ref1, bin1, ref2, bin2) %>%
         dplyr::mutate(count = 0)
}

makeSimData <- function(refs, refSize, maxNReads, diagSize, gaussianNoise, splitRefs, splitPositions, joinRef1, joinRef2, joinDistances) {
    interactionMatrices <- makeEmptyInteractionMatrices(refs, refSize)
    interactionMatrices <- addDiagonalMatrices(interactionMatrices, refs, refSize, maxNReads, diagSize, gaussianNoise)
    interactionMatrices <- addSplits(interactionMatrices, refs, splitRefs, splitPositions)
    interactionMatrices <- addJoins(interactionMatrices, refs, refSize, maxNReads, diagSize, joinRef1, joinRef2, joinDistances)
    interactionMatrices <- addNoise(interactionMatrices, gaussianNoise)
    interactionMatrices <- cleanMatrices(interactionMatrices, refs)
    msscafData(inputMatrix  = interactionMatrices,
                    maxLinkRange = NULL)
}

makeSimObject <- function(nData = 1, nRefs = 1, refSize = 100, maxNReads = 20, diagSize = 10, gaussianNoise = 0.2, splitRefs = c(), splitPositions = c(), joinRef1 = c(), joinRef2 = c(), joinDistances = c()) {
    if (length(splitRefs) != length(splitPositions)) {
        stop("# split references and # split positions should be equal.")
    }
    if (length(joinRef1) != length(joinRef2)) {
        stop("# join references #1 and # join references #2 should be equal.")
    }
    if (length(joinRef1) != length(joinDistances)) {
        stop("# join distances and # join references should be equal.")
    }
    if (length(joinDistances) > 0) {
        if (max(joinDistances) >= diagSize) {
            stop("Join distance should be less than diagonal size.")
        }
    }
    refs                <- paste0("ref_", seq.int(nRefs))
    sizes               <- rep(refSize, nRefs)
    names(sizes)        <- refs
    sequences           <- DNAStringSet(rep(paste0(rep("A", refSize), collapse = ""), nRefs))
    names(sequences)    <- refs
    tmpFileName         <- tempfile()
    writeXStringSet(sequences, tmpFileName)
    object              <- msscaf(tmpFileName, 1, 1)
    unlink(tmpFileName)
    for (i in seq.int(nData)) {
        data   <- makeSimData(refs, refSize, maxNReads, diagSize, gaussianNoise, splitRefs, splitPositions, joinRef1, joinRef2, joinDistances)
        object <- addExp(object, data, paste0("sim_", i))
    }
    return(invisible(object))
}
