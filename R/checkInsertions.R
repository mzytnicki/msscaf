checkInsertion <- function(object) {
    rowSums <- object@interactionMatrix %>%
        makeSymmetric() %>%
        group_by(bin1) %>%
        summarise(rowSum = sum(count)) %>%
        ungroup()
    diagonalValue <- object@interactionMatrix %>%
        filter(bin1 == bin2) %>%
        select(bin1, count) %>%
        rename(diagCount = count)
    inner_join(rowSums, diagonalValue, by = "bin1") %>%
        mutate(ratio = diagonalValue / rowSums) %>%
        filter(ratio >= 0.9)
}

normalizeAndInsertion <- function(object) {
    message(paste0("  Working on ", object@chromosome, "."))
    object <- normalizeKR(object)
    checkInsertion(object)
}

checkInsertions <- function(object) {
    message("Splitting matrix.")
    objects <- splitByRef(object)
    breaks <- lapply(objects, normalizeAndInsertion)
    tibble(ref   = map_chr(breaks, "ref"),
           bin   = map_dbl(breaks, "bin"), 
           plot1 = map(breaks, "plot1"),
           plot2 = map(breaks, "plot2"),
           plot3 = map(breaks, "plot3"))
}