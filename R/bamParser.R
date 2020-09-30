parseBamFile <- function(fileName, binSize, nThreads) {
    #RcppParallel::setThreadOptions(numThreads = nThreads);
    l          <- parseBamFileCpp(fileName, binSize, nThreads)
    objectData <- tenxcheckerData(inputMatrix = as_tibble(l$data), binSize = binSize, sizes = l$sizes)
}