splitByRef <- function(object, chromosomes, sizes) {
    data <- object@interactionMatrix %>%
        dplyr::filter(ref1 == ref2) %>%
        dplyr::select(-ref2) %>%
        split(.$ref1) %>%
        purrr::map(~ dplyr::select(.x, -ref1))
    # There might be fewer chromosomes than expected.
    chromosomes <- names(data)
    sizes <- sizes[chromosomes]
    names <- rep.int(object@name, length(sizes))
    purrr::pmap(list(data, chromosomes, sizes, names),
         msscafRefExp,
         parameters = object@parameters)
}

extractRef <- function(object, keptRef, size) {
    data <- object@interactionMatrix %>%
        dplyr::filter(ref1 == ref2) %>%
        dplyr::filter(ref1 == keptRef) %>%
        dplyr::select(-c(ref1, ref2))
    outlierBins <- object@outlierBins %>%
        dplyr::filter(ref == keptRef) %>%
        dplyr::select(bin)
    return(msscafRefExp(data, keptRef, size, object@name, outlierBins, object@parameters))
}

extract2Ref <- function(object, r1, r2, size1, size2) {
    data <- object@interactionMatrix %>%
#       dplyr::mutate(ref1 = as.integer(ref1)) %>%
#       dplyr::mutate(ref2 = as.integer(ref2)) %>%
        dplyr::filter(ref1 == r1, ref2 == r2) %>%
        dplyr::select(-c(ref1, ref2))
    return(msscaf2RefExp(data, r1, r2, size1, size2, object@name, object@parameters))
}

create2Ref <- function(data, object, sizes) {
    if (! is(object, "msscafExp")) {
        stop("Object should be a 'msscafExp'.")
    }
    ref1  <- as.character(data$ref1[[1]])
    ref2  <- as.character(data$ref2[[1]])
    size1 <- sizes[[ref1]]
    size2 <- sizes[[ref2]]
    return(msscaf2RefExp(data %>% dplyr::select(-c(ref1, ref2)),
                              ref1,
                              ref2,
                              size1,
                              size2,
                              object@name,
                              object@parameters))
}

splitBy2Ref <- function(object, sizes) {
    if (! is(object, "msscafExp")) {
        stop("Object should be a 'msscafExp'.")
    }
    breakNCells <- ifelse(is.null(object@parameters@breakNCells), 0, object@parameters@breakNCells)
    pairs <- object@interactionMatrix %>%
        dplyr::filter(ref1 != ref2) %>%
        tidyr::unite("ref1_ref2", ref1, ref2, remove = FALSE) %>%
        dplyr::group_by(ref1_ref2) %>%
        dplyr::summarise(nCounts = n()) %>%
        dplyr::filter(nCounts >= breakNCells) %>%
        dplyr::pull(ref1_ref2)
    #message(paste0("\tKeeping ", length(pairs), " pairs of references."))
    data <- object@interactionMatrix %>%
        dplyr::filter(ref1 != ref2) %>%
        tidyr::unite("ref1_ref2", ref1, ref2, remove = FALSE) %>%
        dplyr::filter(ref1_ref2 %in% pairs) %>%
        dplyr::group_by(ref1_ref2) %>%
        dplyr::group_split() %>%
        purrr::map(~ dplyr::select(.x, -ref1_ref2))
    names(data) <- pairs
    lapply(data, create2Ref, object = object, sizes = sizes)
}

splitBy2RefFromMatrix <- function(object, interactionMatrix, sizes) {
    if (! is(object, "msscafExp")) {
        stop("Object should be a 'msscafExp'.")
    }
    data <- interactionMatrix %>%
        dplyr::filter(ref1 != ref2) %>%
        tidyr::unite("ref1_ref2", ref1, ref2, remove = FALSE) %>%
        dplyr::group_by(ref1_ref2)
    splitData <- dplyr::group_split(data)
    names(splitData) <- dplyr::group_keys(data) %>% dplyr::pull(ref1_ref2)
    purrr::map(splitData, create2Ref, object = object, sizes = sizes)
}

getDataset <- function(object, name) {
    if (! is(object, "msscafClass")) {
        stop("Object should be a 'msscafClass'.")
    }
    index <- purrr::detect_index(object@data, function(d) d@name == name)
    if (index == 0) {
        stop(paste0("Dataset name '", name, "' not found."))
    }
    return(object@data[[index]])
}

getDatasetRef <- function(object, datasetName, ref1, ref2 = NULL) {
    if (! is(object, "msscafClass")) {
        stop("Object should be a 'msscafClass'.")
    }
    dataset <- getDataset(object, datasetName)
    if (is.null(ref2)) {
        return(extractRef(dataset, ref1, object@sizes[[ref1]], object@name))
    }
    return(extract2Ref(dataset, ref1, ref2, object@sizes[[ref1]], object@sizes[[ref2]], object@name))
}
