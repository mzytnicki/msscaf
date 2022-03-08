parseBamFile <- function(fileName, binSize) {
    tenxcheckerData(inputMatrix = as_tibble(parseBamFileCpp(fileName, binSize)))
}
