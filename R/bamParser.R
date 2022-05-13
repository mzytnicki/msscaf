parseBamFile <- function(fileName, binSize) {
    msscafData(inputMatrix = tibble::as_tibble(parseBamFileCpp(fileName, binSize)))
}
