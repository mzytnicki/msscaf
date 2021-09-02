parseHicFile <- function(fileName, binSize) {
    # as_tibble(parseHicCpp(fileName, binSize))
    l <- parseHicCpp(fileName, binSize)
    data <- as_tibble(l$data) %>%
        arrange(ref1, bin1, ref2, bin2) %>%
        drop_na()
    objectData <- tenxcheckerData(inputMatrix = data)
    return(objectData)
}
