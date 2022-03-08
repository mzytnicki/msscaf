parsePafFile <- function(fileName, binSize, minAlnLen, minCount, minNCells) {
    tenxcheckerData(inputMatrix = as_tibble(parsePafCpp(fileName, binSize, minAlnLen, minCount, minNCells)))
}
