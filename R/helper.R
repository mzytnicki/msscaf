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
        r <- tibble::tibble(all = gtools::mixedsort(unique(c(refs1, refs2)))) %>%
            dplyr::mutate(name1 = dplyr::if_else(all %in% refs1, all, NA_character_)) %>%
            dplyr::mutate(name2 = dplyr::if_else(all %in% refs2, all, NA_character_))
        return(r)
    }
    # Tweak the case
    r1 <- tibble::tibble(name1 = refs1, key = stringr::str_to_lower(refs1))
    r2 <- tibble::tibble(name2 = refs2, key = stringr::str_to_lower(refs2))
    r <- dplyr::full_join(r1, r2, by = "key")
    if (nrow(r) > 1.2 * length(refs1)) {
        stop("Cannot merge references.\nPlease check that chromosomes are similar.")
    }
    r %>%
        dplyr::mutate(all = dplyr::if_else(is.na(name1), name2, name1)) %>%
        dplyr::select(-key)
}

# Add the ref sizes to the previously ref sizes
# Return the new merged set
mergeSizes <- function(sizes1, sizes2, refs) {
    if (is.null(sizes2)) {
        return(sizes1)
    }
    sizes <- refs %>%
        dplyr::left_join(tibble::enframe(sizes1, name = "name1", value = "size1"), by = "name1") %>%
        dplyr::left_join(tibble::enframe(sizes2, name = "name2", value = "size2"), by = "name2") %>%
        dplyr::mutate(size = pmax(size1, size2, na.rm = TRUE)) %>%
	dplyr::select(c("all", "size")) %>%
        tibble::deframe()
    sizes <- sizes[gtools::mixedsort(names(sizes))]
    return(sizes)
}

# Update refs in experiments, so that it matches the genome.
updateRefs <- function(object, data) {
    if (! is(object, "msscafClass")) {
        stop("Object should be a 'msscafClass'.")
    }
    if (! is(data, "msscafData")) {
        stop("Object should be a 'msscafData'.")
    }
    oldLevels <- levels(data@inputMatrix$ref1)
    oldLevels <- tibble::tibble(old = oldLevels, key = stringr::str_to_lower(oldLevels))
    newLevels <- object@chromosomes
    newLevels <- tibble::tibble(new = newLevels, key = stringr::str_to_lower(newLevels))
    if (! all(unlist(purrr::map(oldLevels$key, ~ .x %in% newLevels$key)))) {
        stop("Not all references are known references.")
    }
    transLevels <- dplyr::left_join(oldLevels, newLevels, by = "key") %>% dplyr::pull(new)
    addLevels   <- dplyr::anti_join(newLevels, oldLevels, by = "key") %>% dplyr::pull(new)
    levels(data@inputMatrix$ref1)         <- transLevels
    data@inputMatrix$ref1                 <- forcats::fct_expand(data@inputMatrix$ref1, addLevels)
    data@inputMatrix$ref1                 <- forcats::fct_relevel(data@inputMatrix$ref1, object@chromosomes)
    levels(data@inputMatrix$ref2)         <- transLevels
    data@inputMatrix$ref2                 <- forcats::fct_expand(data@inputMatrix$ref2, addLevels)
    data@inputMatrix$ref2                 <- forcats::fct_relevel(data@inputMatrix$ref2, object@chromosomes)
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

computeRefSizes <- function(object) {
    return(computeRefSizesCpp(object@interactionMatrix))
}

.keepScaffolds <- function(object, chromosomes) {
    if (! is(object, "msscafExp")) {
        stop("Parameter should be a msscafExp.")
    }
    object@outlierBins <- object@outlierBins %>%
        dplyr::filter(ref %in% chromosomes) %>%
        dplyr::mutate(ref = as.character(ref)) %>%
        dplyr::mutate(ref = factor(ref, levels = chromosomes))
    object@interactionMatrix <- tibble::as_tibble(keepScaffoldsCpp(object@interactionMatrix, chromosomes))
    return(object)
}

keepScaffolds <- function(object, chromosomes) {
    if (! is(object, "msscafClass")) {
        stop("Parameter should be a msscafClass.")
    }
    chromosomes        <- gtools::mixedsort(unique(as.character(chromosomes)))
    object@chromosomes <- chromosomes
    object@sizes       <- object@sizes[chromosomes]
    object@data        <- purrr::map(object@data, .keepScaffolds, chromosomes = chromosomes)
    object@sequences   <- object@sequences[chromosomes]
    return(object)
}

# Transforms a sparse matrix (bin1, bin2) to full (bin1, distance)
#   maxDistance: maximum distance wrt to the diagonal or corner
fillCorner <- function(interactionMatrix, maxDistance, isCorner = FALSE) {
    bins <- seq.int(from = 0, to = maxDistance, by = 1)
    d    <- function(bin1, bin2, isCorner) {
        if (isCorner) {
            return(bin1 + bin2)
        }
        return(bin1 - bin2)
    }
    tibble::tibble(bin1 = bins, bin2 = bins) %>%
        tidyr::expand(bin1, bin2) %>%
        dplyr::mutate(distance = d(bin1, bin2, isCorner)) %>%
        dplyr::filter(distance >= 0) %>%
        dplyr::filter(distance <= maxDistance) %>%
        dplyr::left_join(interactionMatrix, by = c("bin1", "bin2")) %>%
        tidyr::replace_na(list(count = 0))
}
