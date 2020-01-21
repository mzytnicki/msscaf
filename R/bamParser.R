parseBamFile <- function(fileName, binSize) {
    as_tibble(parseBamFileCpp(fileName, binSize))
}