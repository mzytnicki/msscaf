parseHicFile <- function(fileName, binSize) {
    # as_tibble(parseHicCpp(fileName, binSize))
    as_tibble(parseHicCpp(fileName, binSize)) %>%
        mutate(ref1 = fct_relevel(ref1, str_sort(levels(ref1)))) %>%
        mutate(ref2 = fct_relevel(ref2, str_sort(levels(ref2)))) %>%
        arrange(ref1, bin1, ref2, bin2) %>%
        drop_na()
}