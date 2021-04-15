generateRow <- function(distance, maxNReads, size, diagSize) {
    lambda <- as.integer(round(maxNReads * dgeom(distance, prob = 1 / diagSize)))
    tibble(bin1 = seq.int(size), count = as.integer(round(rpois(size, lambda)))) %>%
        dplyr::mutate(bin2 = as.integer(bin1 - distance)) %>%
        dplyr::filter(bin2 >= 0) %>%
        dplyr::select(bin1, bin2, count) %>%
        dplyr::filter(count > 0)
}

makePlainMatrix <- function(size, maxNReads, diagSize) {
    purrr::map_dfr(seq.int(size) - 1, generateRow, maxNReads, size, diagSize) %>%
        dplyr::mutate(ref1 = factor(c("dummyRef")), .before = 1) %>%
        dplyr::mutate(ref2 = factor(c("dummyRef")), .after = bin1)
}

makeNoiseMatrix <- function(data, size, gaussianSd) {
    noise <- tibble(bin1 = seq.int(size) - 1, bin2 = seq.int(size) - 1) %>%
        tidyr::complete(bin1, bin2) %>%
        dplyr::mutate(count = as.integer(round(abs(rnorm(nrow(.), mean = 0, sd = gaussianSd))))) %>%
        dplyr::filter(count > 0) %>%
        dplyr::mutate(ref1 = factor(c("dummyRef")), .before = 1) %>%
        dplyr::mutate(ref2 = factor(c("dummyRef")), .after = bin1)
    data %>% bind_rows(noise) %>%
        dplyr::group_by(ref1, bin1, ref2, bin2) %>%
        dplyr::summarise(count = sum(count)) %>%
        dplyr::ungroup()
}

makePlainNoisedData <- function(size, maxNReads, diagSize, gaussianNoise) {
    data <- makePlainMatrix(size, maxNReads, diagSize)
    data <- makeNoiseMatrix(data, size, gaussianNoise)
    object <- tenxcheckerData(inputMatrix  = data,
                              binSize      = 1,
                              maxLinkRange = NULL,
                              sizes        = c("dummyRef" = size))
    return(invisible(object))
}

cleanMatrices <- function(interactionMatrices, refs) {
    interactionMatrices %>%
        dplyr::filter(count > 0) %>%
        dplyr::filter(ref1 != ref2 | bin1 >= bin2) %>%
        dplyr::mutate(ref1 = factor(ref1, levels = refs)) %>%
        dplyr::mutate(ref2 = factor(ref2, levels = refs))
}

generateCorner <- function(distance, r1, r2, refSize, maxNReads, diagSize) {
    nRefs  <- length(refs)
    lambda <- as.integer(round(maxNReads * dgeom(distance, prob = 1 / diagSize)))
    tibble(bin1 = seq.int(from = refSize - distance, to = refSize - 1),
           count = as.integer(round(rpois(distance, lambda)))) %>%
        dplyr::mutate(bin2 = as.integer(bin1 + distance)) %>%
        dplyr::filter(bin2 <= refSize) %>%
        dplyr::select(bin1, bin2, count) %>%
        dplyr::mutate(ref1 = r1) %>%
        dplyr::mutate(ref2 = r2)
}

addJoin <- function(interactionMatrices, refSize, maxNReads, ref1, ref2, distance, diagSize) {
    corner <- purrr::map_dfr(seq.int(from = distance, to = diagSize), generateCorner, ref1, ref2, refSize, maxNReads, diagSize)
    dplyr::bind_rows(interactionMatrices, diags) %>%
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
        dplyr::mutate(count = if_else(ref1 == ref & ref2 == ref & bin1 > position & bin2 < position, 0, count))
}

addSplits <- function(interactionMatrices, refs, positions) {
    for (i in seq_along(refs)) {
        interactionMatrices <- addSplit(interactionMatrices, refs[[i]], positions[[i]])
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
    tibble(bin1 = rep(seq.int(refSize) - 1, nRefs),
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
        dplyr::summarise(count = sum(count), .groups = "drop")
}

makeEmptyInteractionMatrices <- function(refs, refSize) {
    bins    <- as.integer(seq.int(refSize) - 1)
    refBins <- purrr::cross2(refs, bins)
    tibble(ref1 = refBins %>% purrr::map(1) %>% purrr::flatten_chr(),
           bin1 = refBins %>% purrr::map(2) %>% purrr::flatten_int()) %>%
         dplyr::mutate(ref2 = ref1) %>%
         dplyr::mutate(bin2 = bin1) %>%
         tidyr::expand(ref1, bin1, ref2, bin2) %>%
         dplyr::mutate(count = 0)
}

makeSimData <- function(refs, refSize, maxNReads, diagSize, gaussianNoise, splitRefs, splitPositions, joinRef1, joinRef2, joinDistances) {
    interactionMatrices <- makeEmptyInteractionMatrices(refs, refSize)
    interactionMatrices <- addDiagonalMatrices(interactionMatrices, refs, refSize, maxNReads, diagSize, gaussianNoise)
    interactionMatrices <- addSplits(interactionMatrices, splitRefs, splitPositions)
    interactionMatrices <- addJoins(interactionMatrices, refs, refSize, maxNReads, diagSize, joinRef1, joinRef2, joinDistances)
    interactionMatrices <- addNoise(interactionMatrices, gaussianNoise)
    interactionMatrices <- cleanMatrices(interactionMatrices, refs)
    tenxcheckerData(inputMatrix  = interactionMatrices,
                    maxLinkRange = NULL)
}

makeSimObject <- function(nData, nRefs, refSize, maxNReads, diagSize, gaussianNoise, splitRefs, splitPositions, joinRef1, joinRef2, joinDistances) {
    refs                <- paste0("ref_", seq.int(nRefs))
    sizes               <- rep(refSize, nRefs)
    names(sizes)        <- refs
    sequences           <- DNAStringSet(rep(paste0(rep("A", refSize), collapse = ""), nRefs))
    names(sequences)    <- refs
    tmpFileName         <- tempfile()
    writeXStringSet(sequences, tmpFileName)
    object              <- tenxchecker(tmpFileName, 1, 1)
    unlink(tmpFileName)
    for (i in seq.int(nData)) {
        data   <- makeSimData(refs, refSize, maxNReads, diagSize, gaussianNoise, splitRefs, splitPositions, joinRef1, joinRef2, joinDistances)
        object <- addExp(object, data, paste0("sim_", i))
    }
    return(invisible(object))
}
