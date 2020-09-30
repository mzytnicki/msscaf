parsePafFile <- function(fileName, binSize, minAlnLen, minCount, minNCells) {
    #as_tibble(parsePafCpp(fileName, binSize, minAlnLen, minCount, minNCells))
    l          <- parsePafCpp(fileName, binSize, minAlnLen, minCount, minNCells)
    objectData <- tenxcheckerData(inputMatrix = as_tibble(l$data), binSize = binSize, maxLinkRange = l$size, sizes = l$sizes)
    return(objectData)
}
