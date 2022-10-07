parseHicFile <- function(fileName, binSize) {
    parseHicCpp(fileName, binSize) %>%
        tibble::as_tibble() %>%
        dplyr::arrange(ref1, bin1, ref2, bin2) %>%
        tidyr::drop_na() %>%
        msscafData()
}
