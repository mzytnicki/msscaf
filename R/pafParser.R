parsePafFile <- function(fileName, binSize, minAlnLen = 500, minCount = 10, minNCells = 10) {
    msscafData(inputMatrix = tibble::as_tibble(parsePafCpp(fileName, binSize, minAlnLen, minCount, minNCells)))
}
