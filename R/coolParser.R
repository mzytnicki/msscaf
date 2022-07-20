parseCoolFile <- function(fileName, binSize) {
    uri <- function(path) {
        if (!is.numeric(binSize)) return(path)
        return(
            paste(
                "resolutions",
                format(binSize, scientific = FALSE),
                path,
                sep = "/"
            )
        )
    }
    bins <- tibble::tibble(
        ref   = rhdf5::h5read(file = fileName, name = uri("bins/chrom")),
        start = rhdf5::h5read(file = fileName, name = uri("bins/start"))) %>%
        dplyr::mutate(ref = factor(ref)) %>%
        dplyr::mutate(start = as.integer(start)) %>%
        dplyr::mutate(bin = start / binSize) %>%
        tibble::rowid_to_column(var = "id") %>%
        dplyr::mutate(id = id - 1) %>%
        dplyr::select(id, ref, bin)
    tibble::tibble(
        id1   = rhdf5::h5read(file = fileName, name = uri("pixels/bin1_id")),
        id2   = rhdf5::h5read(file = fileName, name = uri("pixels/bin2_id")),
        count = rhdf5::h5read(file = fileName, name = uri("pixels/count"))) %>%
        dplyr::mutate(id1 = as.integer(id1)) %>%
        dplyr::mutate(id2 = as.integer(id2)) %>%
        dplyr::left_join(bins, by = c("id1" = "id")) %>%
        dplyr::rename(ref1 = ref, bin1 = bin) %>%
        dplyr::left_join(bins, by = c("id2" = "id"), suffix = c("", "2")) %>%
        dplyr::rename(ref2 = ref, bin2 = bin) %>%
        dplyr::select(ref1, bin1, ref2, bin2, count) %>%
        msscafData()
}