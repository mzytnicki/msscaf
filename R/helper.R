splitByRef <- function(object) {
    data <- object@interactionMatrix %>%
        filter(ref1 == ref2) %>%
        select(-ref2) %>%
        split(.$ref1) %>%
        map(~ select(.x, -ref1))
    pmap(data,
         object@chromosomes,
         object@sizes,
         tenxcheckerRefExp,
         parameters = object@parameters)
}

create2Ref <- function(data, object = object) {
    ref1  <- as.character(data$ref1[[1]])
    ref2  <- as.character(data$ref2[[1]])
    size1 <- object@sizes[[ref1]]
    size2 <- object@sizes[[ref2]]
    return(tenxchecker2RefExp(data %>% select(-c(ref1, ref2)),
                              ref1,
                              ref2,
                              size1,
                              size2,
                              object@parameters))
}

splitBy2Ref <- function(object) {
    data <- object@interactionMatrix %>%
        filter(ref1 != ref2) %>%
        unite("ref1_ref2", ref1, ref2, remove = FALSE) %>%
        group_by(ref1_ref2) %>%
        group_split() %>%
        map(~ select(.x, -ref1_ref2))
    lapply(data, create2Ref, object = object)
}

computeRefSizes <- function(object) {
    bind_rows(object@interactionMatrix %>%
                  select(c(ref1, bin1)) %>%
                  rename(ref = ref1, bin = bin1),
              object@interactionMatrix %>%
                  select(c(ref2, bin2)) %>%
                  rename(ref = ref2, bin = bin2)) %>%
        group_by(ref) %>%
        summarise(size = max(bin)) %>%
        deframe()
}

makeSymmetric <- function(data) {
    if ("ref1" %in% colnames(data)) {
        data %<>%
            bind_rows(data %>%
                          rename(tmp = ref1) %>%
                          rename(ref1 = ref2) %>%
                          rename(ref2 = tmp) %>%
                          rename(tmp = bin1) %>%
                          rename(bin1 = bin2) %>%
                          rename(bin2 = tmp))
        return(data)
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