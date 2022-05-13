getNInteractionsPerRef <- function(object) {
    object@interactionMatrix %>%
        dplyr::filter(ref1 == ref2) %>%
        dplyr::rename(ref = ref1) %>%
        dplyr::count(ref, name = "nInteractions")
}

# Isolate one chromosome.
# Generate a (spurious) split.
# Return the new object (with 2 chromosomes)
generateSplit <- function(object) {
    minNInteractions <- 100
    minSize          <- 100
    # Select reference with a sufficient number of interactions
    chromosomes <- object@data %>%
        purrr::map_dfr(getNInteractionsPerRef, .id = "id") %>%
        dplyr::group_by(ref) %>%
        dplyr::summarize(nInteractions = min(nInteractions), .groups = "drop") %>%
        dplyr::filter(nInteractions > minNInteractions) %>%
        dplyr::mutate(size = object@sizes[ref]) %>%
        dplyr::filter(size > 2 * minSize) %>%
        dplyr::pull(ref)
    selectedChromosome <- sample(chromosomes, 1)
    size               <- object@sizes[[selectedChromosome]]
    object             <- keepScaffolds(object, selectedChromosome)
    selectedSplit      <- round(runif(1, min = minSize, max = size - minSize))
    breaks             <- tibble::tibble(ref = factor(selectedChromosome), bin = selectedSplit)
    object@breaks      <- breaks
    object             <- splitChromosomes(object)
    return(invisible(object))
}


# Isolate two chromosomes.
# Generate a (spurious) join.
# Return the new object (with 1 chromosome)
generateJoin <- function(object) {
    minNInteractions <- 100
    minSize          <- 100
    # Select reference with a sufficient number of interactions
    chromosomes <- object@data %>%
        purrr::map_dfr(getNInteractionsPerRef, .id = "id") %>%
        dplyr::group_by(ref) %>%
        dplyr::summarize(nInteractions = min(nInteractions), .groups = "drop") %>%
        dplyr::filter(nInteractions > minNInteractions) %>%
        dplyr::mutate(size = object@sizes[ref]) %>%
        dplyr::filter(size > 2 * minSize) %>%
        dplyr::pull(ref)
    selectedChromosomes <- sample(chromosomes, 2)
    message(str(paste0("Selected refs ", paste(selectedChromosomes, collapse = ", "), ", of size ", paste(object@sizes[selectedChromosomes], collapse = ", "))))
    object              <- keepScaffolds(object, selectedChromosomes)
    object@joins        <- tibble::tibble(ref1 = selectedChromosomes[[1]], ref2 = selectedChromosomes[[2]], after1 = TRUE, after2 = FALSE, pvalue = 0.0, pvalueCorner = 0.0)
    object              <- scaffold(object)
    return(invisible(object))
}

