# Add the new references to the previously seen references
# Potentially change case
# Return a tibble, where the r1 is the set of the seen refs,
#    r2, the set of new refs
mergeRefs <- function(refs1, refs2) {
    # If one set is contained in the other one, get the biggest
    # Or if the two sets are nearly the same, get the union
    if ((all(refs2 %in% refs1)) |
        (all(refs1 %in% refs2)) |
        (length(intersect(refs1, refs2)) >= 0.8 * length(refs1))) {
        r <- tibble(all = mixedsort(unique(c(refs1, refs2)))) %>%
            dplyr::mutate(name1 = if_else(all %in% refs1, all, NA_character_)) %>%
            dplyr::mutate(name2 = if_else(all %in% refs2, all, NA_character_))
        return(r)
    }
    # Tweak the case
    r1 <- tibble(name1 = refs1, key = str_to_lower(refs1))
    r2 <- tibble(name2 = refs2, key = str_to_lower(refs2))
    r <- full_join(r1, r2, by = "key")
    if (nrow(r) > 1.2 * length(refs1)) {
        stop("Cannot merge references.\nPlease check that chromosomes are similar.")
    }
    r %>%
        dplyr::mutate(all = if_else(is.na(name1), name2, name1)) %>%
        dplyr::select(-key)
}

# Add the ref sizes to the previously ref sizes
# Return the new merged set
mergeSizes <- function(sizes1, sizes2, refs) {
    if (is.null(sizes2)) {
        return(sizes1)
    }
    sizes <- refs %>%
        left_join(enframe(sizes1, name = "name1", value = "size1"), by = "name1") %>%
        left_join(enframe(sizes2, name = "name2", value = "size2"), by = "name2") %>%
        dplyr::mutate(size = pmax(size1, size2, na.rm = TRUE)) %>%
	dplyr::select(c("all", "size")) %>%
        tibble::deframe()
    sizes <- sizes[mixedsort(names(sizes))]
    return(sizes)
}

# Update refs in experiments, so that it matches the genome.
updateRefs <- function(object, data) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Object should be a 'tenxcheckerClass'.")
    }
    if (! is(data, "tenxcheckerData")) {
        stop("Object should be a 'tenxcheckerData'.")
    }
    oldLevels <- levels(data@inputMatrix$ref1)
    oldLevels <- tibble(old = oldLevels, key = str_to_lower(oldLevels))
    newLevels <- object@chromosomes
    newLevels <- tibble(new = newLevels, key = str_to_lower(newLevels))
    if (! all(unlist(map(oldLevels$key, ~ .x %in% newLevels$key)))) {
        stop("Not all references are known references.")
    }
    transLevels <- left_join(oldLevels, newLevels, by = "key") %>% pull(new)
    addLevels   <- anti_join(newLevels, oldLevels, by = "key") %>% pull(new)
    levels(data@inputMatrix$ref1)         <- transLevels
    data@inputMatrix$ref1                 <- fct_expand(data@inputMatrix$ref1, addLevels)
    data@inputMatrix$ref1                 <- fct_relevel(data@inputMatrix$ref1, object@chromosomes)
    levels(data@inputMatrix$ref2)         <- transLevels
    data@inputMatrix$ref2                 <- fct_expand(data@inputMatrix$ref2, addLevels)
    data@inputMatrix$ref2                 <- fct_relevel(data@inputMatrix$ref2, object@chromosomes)
    # Modify the data, so that ref1 >= ref2
    inverted <- data@inputMatrix %>%
        dplyr::filter(as.integer(ref1) < as.integer(ref2)) %>%
        dplyr::rename(ref1 = ref2, ref2 = ref1, bin1 = bin2, bin2 = bin1)
    data@inputMatrix <- data@inputMatrix %>%
        dplyr::filter(as.integer(ref1) >= as.integer(ref2)) %>%
        dplyr::bind_rows(inverted)
    return(data)
}

normalizeRefs <- function(object) {
    levels(object@interactionMatrix$ref1) <- stringr::str_to_lower(levels(object@interactionMatrix$ref1))
    levels(object@interactionMatrix$ref2) <- stringr::str_to_lower(levels(object@interactionMatrix$ref2))
    object@chromosomes <- stringr::str_to_lower(object@chromosomes)
    names(object@sizes) <- stringr::str_to_lower(names(object@sizes))
    return(object)
}

splitByRef <- function(object, chromosomes, sizes) {
    data <- object@interactionMatrix %>%
        filter(ref1 == ref2) %>%
        dplyr::select(-ref2) %>%
        split(.$ref1) %>%
        map(~ dplyr::select(.x, -ref1))
    # There might be fewer chromosomes than expected.
    chromosomes <- names(data)
    sizes <- sizes[chromosomes]
    pmap(list(data, chromosomes, sizes),
         tenxcheckerRefExp,
         parameters = object@parameters)
}

extractRef <- function(object, ref, size) {
    data <- object@interactionMatrix %>%
        dplyr::filter(ref1 == ref2) %>%
        dplyr::filter(ref1 == ref) %>%
        dplyr::select(-c(ref1, ref2))
    return(tenxcheckerRefExp(data, ref, size, object@parameters))
}

extract2Ref <- function(object, r1, r2, size1, size2) {
    data <- object@interactionMatrix %>%
        filter(ref1 == r1, ref2 == r2) %>%
        dplyr::select(-c(ref1, ref2))
    return(tenxchecker2RefExp(data, r1, r2, size1, size2, object@parameters))
}

create2Ref <- function(data, object, sizes) {
    if (! is(object, "tenxcheckerExp")) {
        stop("Object should be a 'tenxcheckerExp'.")
    }
    ref1  <- as.character(data$ref1[[1]])
    ref2  <- as.character(data$ref2[[1]])
    size1 <- sizes[[ref1]]
    size2 <- sizes[[ref2]]
    return(tenxchecker2RefExp(data %>% dplyr::select(-c(ref1, ref2)),
                              ref1,
                              ref2,
                              size1,
                              size2,
                              object@parameters))
}

splitBy2Ref <- function(object, sizes) {
    if (! is(object, "tenxcheckerExp")) {
        stop("Object should be a 'tenxcheckerExp'.")
    }
    breakNCells <- ifelse(is.null(object@parameters@breakNCells), 0, object@parameters@breakNCells)
    pairs <- object@interactionMatrix %>%
        filter(ref1 != ref2) %>%
        unite("ref1_ref2", ref1, ref2, remove = FALSE) %>%
        group_by(ref1_ref2) %>%
        summarise(nCounts = n(), maxCounts = max(count)) %>%
        #filter(maxCounts >= object@parameters@maxLinkRange) %>%
        filter(nCounts >= breakNCells) %>%
        pull(ref1_ref2)
    message(paste0("\tKeeping ", length(pairs), " pairs of references."))
    data <- object@interactionMatrix %>%
        filter(ref1 != ref2) %>%
        unite("ref1_ref2", ref1, ref2, remove = FALSE) %>%
        filter(ref1_ref2 %in% pairs) %>%
        group_by(ref1_ref2) %>%
        group_split() %>%
        map(~ dplyr::select(.x, -ref1_ref2))
    names(data) <- pairs
    lapply(data, create2Ref, object = object, sizes = sizes)
}

computeRefSizes <- function(object) {
    return(computeRefSizesCpp(object@interactionMatrix))
}

isMatrixEmpty <- function(data) {
    return(nrow(data) == 0)
}

.keepScaffolds <- function(object, chromosomes) {
    if (! is(object, "tenxcheckerExp")) {
        stop("Parameter should be a tenxcheckerExp.")
    }
    object@interactionMatrix <- as_tibble(keepScaffoldsCpp(object@interactionMatrix, chromosomes))
    return(object)
}

keepScaffolds <- function(object, chromosomes) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    chromosomes        <- mixedsort(unique(chromosomes))
    object@chromosomes <- chromosomes
    object@sizes       <- object@sizes[chromosomes]
    object@data        <- map(object@data, .keepScaffolds, chromosomes = chromosomes)
    return(object)
}

makeSymmetric <- function(data) {
    if ("ref1" %in% colnames(data)) {
        data %<>%
            dplyr::bind_rows(data %>%
                          dplyr::rename(tmp = ref1) %>%
                          dplyr::rename(ref1 = ref2) %>%
                          dplyr::rename(ref2 = tmp) %>%
                          dplyr::rename(tmp = bin1) %>%
                          dplyr::rename(bin1 = bin2) %>%
                          dplyr::rename(bin2 = tmp))
        return(data)
    }
    data %<>%
        dplyr::bind_rows(data %>%
                      dplyr::rename(tmp = bin1) %>%
                      dplyr::rename(bin1 = bin2) %>%
                      dplyr::rename(bin2 = tmp))
    return(data)
}

makeFullMatrix <- function(data) {
    if (nrow(data) == 0) {
        stop("Matrix is empty.")
    }
    mat <- sparseMatrix(data$bin1 + 1, data$bin2 + 1, x = data$count, symmetric = TRUE)
    return(mat)
}

makeFullMatrixGenome <- function(data, sizes) {
    message(str(data %>% dplyr::filter(ref1 != ref2)))
    message(str(data %>% dplyr::filter(as.integer(ref1) > as.integer(ref2))))
    if (data %>% dplyr::filter(as.integer(ref1) < as.integer(ref2)) %>% nrow() > 0) stop("Input matrix is not upper: problem with references.")
    if (data %>% dplyr::filter(ref1 == ref2, bin1 < bin2) %>% nrow() > 0) stop("Input matrix is not upper: problem with bins.")
    if (data %>% dplyr::mutate(size1 = sizes[ref1], size2 = sizes[ref2]) %>% dplyr::filter(bin1 > size1 | bin2 > size2) %>% nrow() > 0) stop("Problem with ref sizes.")
    cumulatedSizes <- cumsum(sizes)
    names(cumulatedSizes) <- NULL
    cumulatedSizes <- c(1, cumulatedSizes)
    message(str(cumulatedSizes))
    data <- data %>%
        dplyr::mutate(bin1 = bin1 + cumulatedSizes[ref1]) %>%
        dplyr::mutate(bin2 = bin2 + cumulatedSizes[ref2]) %>%
        dplyr::select(-c(ref1, ref2))
    message(str(data %>% dplyr::filter(bin1 < bin2)))
    if (data %>% dplyr::filter(bin1 < bin2) %>% nrow() > 0) stop("Matrix is not upper.")
    mat <- sparseMatrix(data$bin1, data$bin2, x = data$count, symmetric = TRUE)
    return(mat)
}


makeSparseMatrix <- function(data) {
    data <- summary(data)
    tibble(bin1  = data$i - 1, bin2  = data$j - 1, count = data$x) %>% filter(bin1 >= bin2)
}


makeSparseMatrixGenome <- function(data, sizes) {
    refs           <- names(sizes)
    cumulatedSizes <- cumsum(sizes)
    names(cumulatedSizes) <- NULL
    cumulatedSizes <- c(1, cumulatedSizes)
    sumSizes       <- sum(sizes)
    refSizes       <- enframe(cumulatedSizes, name = "ref", value = "bin") %>%
        right_join(tibble(bin = seq.int(sumSizes)), by = "bin") %>%
        dplyr::mutate(ref = factor(refs[ref], levels = refs)) %>%
        dplyr::arrange(bin) %>%
        tidyr::fill(ref)
    data           <- summary(data)
    tibble(bin1  = data$i, bin2  = data$j, count = data$x) %>%
        dplyr::left_join(refSizes, by = c("bin1" = "bin")) %>%
        dplyr::rename(ref1 = ref) %>%
        dplyr::left_join(refSizes, by = c("bin2" = "bin")) %>%
        dplyr::rename(ref2 = ref) %>%
        dplyr::mutate(bin1 = bin1 - cumulatedSizes[ref1]) %>%
        dplyr::mutate(bin2 = bin2 - cumulatedSizes[ref2]) %>%
        dplyr::filter(as.integer(ref1) >= as.integer(ref2)) %>%
        dplyr::filter((as.integer(ref1) != as.integer(ref2)) || (bin1 >= bin2)) %>%
        dplyr::select(ref1, bin1, ref2, bin2, count)
}


makeFullTibble <- function(data) {
    data %>%
        makeSymmetric() %>%
        complete(bin1, bin2, fill = list(count = 0))
}

makeFullTibbleNotSquare <- function(data, size1, size2) {
    list(bin1 = seq.int(size1), bin2 = seq.int(size2)) %>%
        cross_df() %>%
        left_join(data, by = c("bin1", "bin2")) %>%
        mutate(count = replace_na(count, 0))
}

makeTibbleFromList <- function(data, n) {
    tibble(bin1 = rep(seq(n), each = n),
           bin2 = rep(seq(n), times = n),
           count = data) %>%
        filter(count != 0) %>%
        filter(bin1 <= bin2)
}

computeScaleFactor <- function(object, sizes) {
    if (is.null(object)) {
        length <- sizes[[2]] - sizes[[1]]
    }
    else if (is(object, "tenxcheckerExp")) {
        length <- sum(sizes)
    }
    else if (is(object, "tenxchecker2RefExp")) {
        length <- min(object@size1, object@size2)
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
        dplyr::mutate(bin1 = as.integer(round(bin1 / scale) * scale)) %>%
        dplyr::mutate(bin2 = as.integer(round(bin2 / scale) * scale))
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

rescaleValue <- function(value, scale) {
    as.integer(round(value / scale) * scale)
}

computeDistanceCount <- function(interactionMatrix, distance) {
    bins <- seq.int(from = 0, to = distance, by = 1)
    emptyMatrix <- tibble(bin1 = bins, bin2 = bins) %>%
        tidyr::expand(bin1, bin2) %>%
        dplyr::filter(bin1 >= bin2) %>%
        dplyr::filter(bin1 - bin2 <= distance)
    emptyMatrix %>%
        dplyr::left_join(interactionMatrix, by = c("bin1", "bin2")) %>%
        tidyr::replace_na(list(count = 0)) %>%
        dplyr::mutate(distance = abs(bin1 - bin2)) %>%
        dplyr::select(distance, count)
}
