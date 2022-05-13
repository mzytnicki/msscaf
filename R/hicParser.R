parseHicFile <- function(fileName, binSize) {
    l <- parseHicCpp(fileName, binSize)
    data <- tibble::as_tibble(l$data) %>%
        dplyr::arrange(ref1, bin1, ref2, bin2) %>%
        tidyr::drop_na()
    objectData <- msscafData(inputMatrix = data)
    return(objectData)
}
