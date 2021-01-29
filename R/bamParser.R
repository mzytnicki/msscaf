parseBamFile <- function(fileName, binSize) {
    l          <- parseBamFileCpp(fileName, binSize)
    objectData <- tenxcheckerData(inputMatrix = as_tibble(l$data), binSize = binSize, sizes = l$sizes)
}