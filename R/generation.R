# Isolate one chromosome.
# Generate a (spurious) split.
# Return the new object (with 2 chromosomes)
generateSplit <- function(object) {
    selectedChromosome <- sample(object@chromosomes, 1)
    object             <- keepScaffolds(object, selectedChromosome)
    selectedSplit      <- round(runif(1, min = 0, max = object@sizes[[1]]))
    breaks             <- tibble(ref = factor(selectedChromosome), bin = selectedSplit)
    object@breaks      <- breaks
    object             <- splitChromosomes(object)
    return(invisible(object))
}

