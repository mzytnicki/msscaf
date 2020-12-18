splitByRef <- function(object) {
    data <- object@interactionMatrix %>%
        filter(ref1 == ref2) %>%
        dplyr::select(-ref2) %>%
        split(.$ref1) %>%
        map(~ dplyr::select(.x, -ref1))
    pmap(list(data, object@chromosomes, object@sizes),
         tenxcheckerRefExp,
         parameters = object@parameters)
}

extractRef <- function(object, ref) {
    data <- object@interactionMatrix %>%
        filter(ref1 == ref2) %>%
        filter(ref1 == ref) %>%
        dplyr::select(-c(ref1, ref2))
    return(tenxcheckerRefExp(data, ref, object@sizes[[ref]], object@parameters))
}

create2Ref <- function(data, object = object) {
    ref1  <- as.character(data$ref1[[1]])
    ref2  <- as.character(data$ref2[[1]])
    size1 <- object@sizes[[ref1]]
    size2 <- object@sizes[[ref2]]
    return(tenxchecker2RefExp(data %>% dplyr::select(-c(ref1, ref2)),
                              ref1,
                              ref2,
                              size1,
                              size2,
                              object@parameters))
}

splitBy2Ref <- function(object) {
    pairs <- object@interactionMatrix %>%
        filter(ref1 != ref2) %>%
        unite("ref1_ref2", ref1, ref2, remove = FALSE) %>%
        group_by(ref1_ref2) %>%
        summarise(nCounts = n(), maxCounts = max(count)) %>%
        #filter(maxCounts >= object@parameters@maxLinkRange) %>%
        filter(nCounts >= object@parameters@breakNCells) %>%
        pull(ref1_ref2)
    message(paste0("Keeping ", length(pairs), " pairs of references."))
    data <- object@interactionMatrix %>%
        filter(ref1 != ref2) %>%
        unite("ref1_ref2", ref1, ref2, remove = FALSE) %>%
        filter(ref1_ref2 %in% pairs) %>%
        group_by(ref1_ref2) %>%
        group_split() %>%
        map(~ dplyr::select(.x, -ref1_ref2))
    names(data) <- pairs
    lapply(data, create2Ref, object = object)
}

computeRefSizes <- function(object) {
    bind_rows(object@interactionMatrix %>%
                  dplyr::select(c(ref1, bin1)) %>%
                  group_by(ref1) %>%
                  summarise(size = max(bin1)) %>%
                  ungroup() %>%
                  dplyr::rename(ref = ref1),
              object@interactionMatrix %>%
                  dplyr::select(c(ref2, bin2)) %>%
                  group_by(ref2) %>%
                  summarise(size = max(bin2)) %>%
                  ungroup() %>%
                  dplyr::rename(ref = ref2)) %>%
        group_by(ref) %>%
        summarise(size = max(size)) %>%
        deframe()
}

isMatrixEmpty <- function(data) {
    return(nrow(data) == 0)
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
    n <- max(data$bin1, data$bin2)
    if (n == -Inf) {
        stop("Matrix is empty.")
    }
    mat <- sparseMatrix(data$bin1 + 1, data$bin2 + 1, x = data$count, symmetric = TRUE)
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

computeScaleFactor <- function(object) {
    if (is(object, "tenxcheckerExp")) {
        length <- sum(object@sizes)
    }
    else if (is(object, "tenxchecker2RefExp")) {
        length <- max(object@size1, object@size2)
    }
    else if (is(object, "tenxcheckerRefExp")) {
        length <- object@size
    }
    else {
        stop("Do not know what to do with object.")
    }
    scaleFactor <- ceiling(log10(length))
    if (scaleFactor <= 3) {
        return(1)
    }
    return(10^(scaleFactor - 3))
}

rescale <- function(data, scale) {
    if (scale == 1) {
        return(data)
    }
    data %<>%
        dplyr::mutate(bin1 = round(bin1 / scale) * scale) %>%
        dplyr::mutate(bin2 = round(bin2 / scale) * scale)
    if ("ref1" %in% colnames(data)) {
        data %<>%
            dplyr::group_by(ref1, ref2, bin1, bin2) %>%
            dplyr::summarise(count = mean(count)) %>%
            dplyr::ungroup()
    }
    else {
        data %<>%
            dplyr::group_by(bin1, bin2) %>%
            dplyr::summarise(count = mean(count)) %>%
            dplyr::ungroup()
    }
    data
}
