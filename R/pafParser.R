parsePafFile <- function(fileName, minAlnLen, minCount, minNCells) {
    #as_tibble(parsePafCpp(fileName, binSize, minAlnLen, minCount, minNCells))
    l          <- parsePafCpp(fileName, binSize, minAlnLen, minCount, minNCells)
    objectData <- tenxcheckerData(inputMatrix = as_tibble(l$data), maxLinkRange = l$size)
    return(objectData)
}
