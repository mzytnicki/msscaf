parseBamFile <- function(fileName, binSize) {
    l          <- parseBamFileCpp(fileName, binSize)
    objectData <- tenxcheckerData(inputMatrix = as_tibble(l$data), maxLinkRange = l$size)
}
