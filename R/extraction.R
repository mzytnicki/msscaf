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
        dplyr::mutate(ref1 = as.integer(ref1)) %>%
        dplyr::mutate(ref2 = as.integer(ref2)) %>%
        dplyr::filter(ref1 == r1, ref2 == r2) %>%
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
    if (! is(object, "tenxcheckerExp")) {
        stop("Object should be a 'tenxcheckerExp'.")
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
    if (! is(object, "tenxcheckerClass")) {
        stop("Object should be a 'tenxcheckerClass'.")
    }
    index <- purrr::detect_index(object@data, function(d) d@name == name)
    if (index == 0) {
        stop(paste0("Dataset name '", name, "' not found."))
    }
    return(object@data[[index]])
}

getDatasetRef <- function(object, datasetName, ref1, ref2 = NULL) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Object should be a 'tenxcheckerClass'.")
    }
    dataset <- getDataset(object, datasetName)
    if (is.null(ref2)) {
        return(extractRef(dataset, ref1, object@sizes[[ref1]]))
    }
    return(extract2Ref(dataset, ref1, ref2, object@sizes[[ref1]], object@sizes[[ref2]]))
}
