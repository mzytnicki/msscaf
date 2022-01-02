# Are the number of bins consistent with the sequence size?
checkSizeDifference <- function(object) {
    sizeDiff <- tibble(binSize = unlist(object@sizes),
        strSize = purrr::map_int(object@sequences, str_length),
        binnedStrSize = strSize %/% object@binSize -
            dplyr::if_else(strSize %% object@binSize == 0, 1, 0)) %>%
        dplyr::filter(binSize != binnedStrSize)
    if (nrow(sizeDiff) != 0) {
        stop(paste("Error! Number of bins and ref. sizes differ.",
            str(sizeDiff), sep = "\n"))
    }
}

checkBinDifference <- function(object, sizes) {
    dplyr::bind_rows(
        object@interactionMatrix %>%
            dplyr::select(ref1, bin1) %>%
            dplyr::rename(ref = ref1, bin = bin1) %>%
            dplyr::distinct(),
        object@interactionMatrix %>%
            dplyr::select(ref2, bin2) %>%
            dplyr::rename(ref = ref2, bin = bin2) %>%
            dplyr::distinct()) %>%
        dplyr::distinct() %>%
        dplyr::mutate(size = sizes[ref]) %>%
        dplyr::filter(bin > size)
}

# Are the bins in the count matrix consistent with the sequence size?
checkAllBinDifference <- function(object) {
    differences <- purrr::map(object@data, checkBinDifference, sizes = object@sizes)
    empty <- differences %>%
        purrr::map(~ nrow(.x) == 0) %>%
        purrr::every(isTRUE)
    if (! empty) {
        stop(paste("Error! Bin of interaction matrix exceeds sizes.",
            str(differences), sep = "\n"))
    }
}

checkMatrix <- function(object) {
    object@interactionMatrix %>%
      dplyr::filter((as.integer(ref1) < as.integer(ref2)) | ((as.integer(ref1) == as.integer(ref2)) & (bin1 < bin2))) %>%
      dplyr::mutate(source = object@name)
}

# Are the matrices up-triangular?
checkMatrices <- function(object) {
    matrices <- purrr::map_dfr(object@data, checkMatrix)
    if (nrow(matrices) > 0) {
        stop(paste("Error! Matrices are not up-triangular.",
            str(matrices), sep = "\n"))
    }
}
