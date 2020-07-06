parseBamFile <- function(fileName, binSize, nThreads) {
    #RcppParallel::setThreadOptions(numThreads = nThreads);
    as_tibble(parseBamFileCpp(fileName, binSize, nThreads))
}