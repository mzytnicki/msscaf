computeInsertionBinDepth <- function(bin    = bin,
                                     object = object,
                                     depth  = depth) {
    start <- bin - depth
    end   <- bin + depth
    interactionMatrix <- object@interactionMatrix %>%
        makeSymmetric() %>%
        filter(bin1 >= bin2)
    triangle1 <- interactionMatrix %>%
        filter(bin1 <= end) %>%
        filter(bin2 >= start) %>%
        pull(count) %>%
        mean()
    triangle2 <- interactionMatrix %>%
        filter(bin2 < start) %>%
        filter(bin1 - bin2 <= depth) %>%
        filter(bin1 + bin2 >= 2 * start) %>%
        pull(count) %>%
        mean()
    triangle3 <- interactionMatrix %>%
        filter(bin1 > end) %>%
        filter(bin1 - bin2 <= depth) %>%
        filter(bin1 + bin2 <= 2 * end) %>%
        pull(count) %>%
        mean()
    otherTriangles <- mean(c(triangle2, triangle3))
    return(list(bin        = bin,
                distance   = depth,
                triangle   = triangle1,
                rest       = otherTriangles,
                difference = triangle1 - otherTriangles))
}

computeInsertionDepth <- function(depth, object = object) {
    purrr::map_df(seq(from = max(1, depth), to = object@size - depth),
           computeInsertionBinDepth,
           object = object,
           depth = depth)
}

computeInsertion <- function(object = object) {
    map_df(seq(from = 0, to = object@parameters@minLinkRange),
           computeInsertionDepth,
           object = object)
}


# checkInsertion <- function(object) {
#     rowSums <- object@interactionMatrix %>%
#         makeSymmetric() %>%
#         group_by(bin1) %>%
#         summarise(rowSum = sum(count)) %>%
#         ungroup()
#     diagonalValue <- object@interactionMatrix %>%
#         filter(bin1 == bin2) %>%
#         select(bin1, count) %>%
#         rename(diagCount = count)
#     inner_join(rowSums, diagonalValue, by = "bin1") %>%
#         mutate(ratio = diagonalValue / rowSums) %>%
#         filter(ratio >= 0.9)
# }

normalizeAndInsertion <- function(object) {
    message(paste0("  Working on ", object@chromosome, "."))
    originalObject <- object
    object <- normalizeKR(object)
    object <- normalizeMD(object)
    object <- removeFarFromDiagonal(object)
    l <- computeInsertion(object = object) %>% drop_na()
    p <- plotInsertions1(l, object@chromosome)
    return(list(data = l, plot = p))
}

checkInsertions <- function(object) {
    message("Splitting matrix.")
    objects <- splitByRef(object)
    insertions <- bplapply(objects, normalizeAndInsertion)
    return(list(data = map_dfr(insertions, "data", .id = "ref") %>%
                    mutate(ref = factor(ref)), 
                plot = map(insertions, "plot")))
}
