splitByRef <- function(object) {
    data <- object@interactionMatrix
    data %<>%
        filter(ref1 == ref2) %>%
        select(-ref2) %>%
        split(.$ref1) %>%
        map(~ select(.x, -ref1))
    map2(data,
         object@chromosomes,
         tenxcheckerRefExp,
         parameters = object@parameters)
}

makeSymmetric <- function(data) {
    if ("ref1" %in% colnames(data)) {
        data %<>%
            bind_rows(data %>%
                          rename(tmp = ref1) %>%
                          rename(ref1 = ref2) %>%
                          rename(ref2 = tmp))
    }
    data %<>%
        bind_rows(data %>%
                      rename(tmp = bin1) %>%
                      rename(bin1 = bin2) %>%
                      rename(bin2 = tmp))
    return(data)
}

makeFullMatrix <- function(data) {
    data <- select(data, c(bin1, bin2, count))
    data$bin1 <- data$bin1 + 1
    data$bin2 <- data$bin2 + 1
    n <- max(data$bin1, data$bin2)
    if (n == -Inf) {
        stop("Matrix is empty.")
    }
    mat <- matrix(0, nrow = n, ncol = n)
    tmp <- as.matrix(data)
    mat[ tmp[, 1:2] ] <- tmp[, 3]
    mat <- mat + t(mat) - diag(diag(mat))
    if (!isSymmetric(mat)) {
        stop("Matrix is not symmetric.")
    }
    return(mat)
}

makeFullTibble <- function(data) {
    data %>%
        makeSymmetric() %>%
        complete(bin1, bin2, fill = list(count = 0))
}

makeTibbleFromList <- function(data, n) {
    tibble(bin1 = rep(seq(n), each = n),
           bin2 = rep(seq(n), times = n),
           count = data) %>%
        filter(count != 0) %>%
        filter(bin1 <= bin2)
}