isMatrixEmpty <- function(data) {
    return(nrow(data) == 0)
}

makeSymmetric <- function(data) {
    if ("ref1" %in% colnames(data)) {
        data %<>%
            dplyr::bind_rows(data %>%
                          dplyr::rename(tmp = ref1) %>%
                          dplyr::rename(ref1 = ref2) %>%
                          dplyr::rename(ref2 = tmp) %>%
                          dplyr::rename(tmp = bin1) %>%
                          dplyr::rename(bin1 = bin2) %>%
                          dplyr::rename(bin2 = tmp)) %>%
                          dplyr::distinct()
        return(data)
    }
    data %>%
        makeSymmetricRefCpp() %>%
        tibble::as_tibble()
}

# Set all the outlier bins to NA
# Possibly overwrite current data, and create new cells
# Data: a symmetric tibble: bin1, bin2, count
# Outliers: a 1-row tibble: bin
# Size: size of the reference
# MinLim: min. bin to consider
# MaxLim: max. bin to consider
removeOutliersRef <- function(data, outliers, size, minLim = -1, maxLim = -1) {
    removeOutliersRefCpp(data, outliers, size, minLim, maxLim) %>% tibble::as_tibble()
}

makeFullMatrix <- function(data) {
    if (nrow(data) == 0) {
        stop("Matrix is empty.")
    }
    mat <- sparseMatrix(data$bin1 + 1, data$bin2 + 1, x = data$count, symmetric = TRUE)
    return(mat)
}

makeFullMatrixGenome <- function(data, sizes) {
    message(str(data %>% dplyr::filter(ref1 != ref2)))
    message(str(data %>% dplyr::filter(as.integer(ref1) > as.integer(ref2))))
    if (data %>% dplyr::filter(as.integer(ref1) < as.integer(ref2)) %>% nrow() > 0) stop("Input matrix is not upper: problem with references.")
    if (data %>% dplyr::filter(ref1 == ref2, bin1 < bin2) %>% nrow() > 0) stop("Input matrix is not upper: problem with bins.")
    if (data %>% dplyr::mutate(size1 = sizes[ref1], size2 = sizes[ref2]) %>% dplyr::filter(bin1 > size1 | bin2 > size2) %>% nrow() > 0) stop("Problem with ref sizes.")
    cumulatedSizes <- cumsum(sizes)
    names(cumulatedSizes) <- NULL
    cumulatedSizes <- c(1, cumulatedSizes)
    message(str(cumulatedSizes))
    data <- data %>%
        dplyr::mutate(bin1 = bin1 + cumulatedSizes[ref1]) %>%
        dplyr::mutate(bin2 = bin2 + cumulatedSizes[ref2]) %>%
        dplyr::select(-c(ref1, ref2))
    message(str(data %>% dplyr::filter(bin1 < bin2)))
    if (data %>% dplyr::filter(bin1 < bin2) %>% nrow() > 0) stop("Matrix is not upper.")
    mat <- sparseMatrix(data$bin1, data$bin2, x = data$count, symmetric = TRUE)
    return(mat)
}


makeSparseMatrix <- function(data) {
    data <- summary(data)
    tibble(bin1  = data$i - 1, bin2  = data$j - 1, count = data$x) %>% filter(bin1 >= bin2)
}


makeSparseMatrixGenome <- function(data, sizes) {
    refs           <- names(sizes)
    cumulatedSizes <- cumsum(sizes)
    names(cumulatedSizes) <- NULL
    cumulatedSizes <- c(1, cumulatedSizes)
    sumSizes       <- sum(sizes)
    refSizes       <- enframe(cumulatedSizes, name = "ref", value = "bin") %>%
        right_join(tibble(bin = seq.int(sumSizes)), by = "bin") %>%
        dplyr::mutate(ref = factor(refs[ref], levels = refs)) %>%
        dplyr::arrange(bin) %>%
        tidyr::fill(ref)
    data           <- summary(data)
    tibble(bin1  = data$i, bin2  = data$j, count = data$x) %>%
        dplyr::left_join(refSizes, by = c("bin1" = "bin")) %>%
        dplyr::rename(ref1 = ref) %>%
        dplyr::left_join(refSizes, by = c("bin2" = "bin")) %>%
        dplyr::rename(ref2 = ref) %>%
        dplyr::mutate(bin1 = bin1 - cumulatedSizes[ref1]) %>%
        dplyr::mutate(bin2 = bin2 - cumulatedSizes[ref2]) %>%
        dplyr::filter(as.integer(ref1) >= as.integer(ref2)) %>%
        dplyr::filter((as.integer(ref1) != as.integer(ref2)) || (bin1 >= bin2)) %>%
        dplyr::select(ref1, bin1, ref2, bin2, count)
}


makeFullTibble <- function(data) {
    data %>%
        makeSymmetric() %>%
        complete(bin1, bin2, fill = list(count = 0))
}

makeFullTibbleNotSquare <- function(data, size1, size2) {
    list(bin1 = seq.int(size1), bin2 = seq.int(size2)) %>%
        cross_df() %>%
        left_join(data, by = c("bin1", "bin2")) %>%
        mutate(count = replace_na(count, 0))
}

makeTibbleFromList <- function(data, n) {
    tibble(bin1 = rep(seq(n), each = n),
           bin2 = rep(seq(n), times = n),
           count = data) %>%
        filter(count != 0) %>%
        filter(bin1 <= bin2)
}
