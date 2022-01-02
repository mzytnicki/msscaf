# .splitChromosome <- function(object, refs, prevRef, newRef, shiftedRef, splitPoint, firstPart) {
#     if (! is(object, "tenxcheckerExp")) {
# 	stop("First parameter should be a tenxcheckerExp.")
#     }
#     #newRef <- factor(newRef, levels = refs)
#     levels(object@interactionMatrix$ref1) <- refs
#     levels(object@interactionMatrix$ref2) <- refs
# #   message(str(prevRef))
# #   message(str(newRef))
# #   message(str(shiftedRef))
# #   message(str(firstPart))
# #   message(str(object@interactionMatrix))
#     splitChromosomeCpp(object@interactionMatrix, prevRef, newRef, shiftedRef, splitPoint, firstPart)
# #   object@interactionMatrix <- object@interactionMatrix %>% 
# #       dplyr::mutate(ref1 = if_else((ref1 == prevRef) & (!comparator(bin1)), newRef, ref1)) %>%
# #       dplyr::mutate(ref2 = if_else((ref2 == prevRef) & (!comparator(bin2)), newRef, ref2)) %>%
# #       dplyr::mutate(bin1 = if_else(ref1 == shiftedRef, as.integer(bin1 - splitPoint + 1), bin1)) %>%
# #       dplyr::mutate(bin2 = if_else(ref2 == shiftedRef, as.integer(bin2 - splitPoint + 1), bin2))
#     return(object)
# }

# splitChromosome <- function(object, parameters) {
#     if (! is(object, "tenxcheckerClass")) {
# 	stop("Parameter should be a tenxcheckerClass.")
#     }
#     prevRef    <- object@chromosomes[[parameters$ref]]
#     splitPoint <- parameters$bin
#     newRef     <- parameters$newRef
# #   message(str(parameters))
# #   message(str(prevRef))
#     prevSize   <- object@sizes[[prevRef]]
#     if (splitPoint >= prevSize) {
# 	stop(paste0("Error while splitting chromosome '", prevRef, "' of size ", prevSize, " at point ", splitPoint, "."))
#     }
#     newSize    <- prevSize - splitPoint
#     firstPart  <- (splitPoint > prevSize / 2)
#     shiftedRef <- if (firstPart) newRef else prevRef
#     if (firstPart) {
# 	comparator <- function(x) { (x < splitPoint) }
#     } else {
# 	comparator <- function(x) { (x > splitPoint) }
#     }
#     prevRefId    <- which(object@chromosomes == prevRef)
#     newRefId     <- which(object@chromosomes == newRef)
#     shiftedRefId <- which(object@chromosomes == shiftedRef)
#     object@data <- map(object@data, .splitChromosome, object@chromosomes, prevRef = prevRefId, newRef = newRefId, shiftedRef = shiftedRefId, splitPoint = splitPoint, firstPart = firstPart)
#     currentRef <- object@sequences[[prevRef]]
#     if (firstPart) {
# 	object@sizes[[prevRef]]     <- splitPoint
# 	object@sizes[[newRef]]      <- newSize
# 	object@sequences[[prevRef]] <- subseq(currentRef, 1, splitPoint * object@binSize)
# 	object@sequences[[newRef]]  <- subseq(currentRef, (splitPoint+1) * object@binSize, length(currentRef))
#     }
#     else {
# 	object@sizes[[newRef]]      <- splitPoint
# 	object@sizes[[prevRef]]     <- newSize
# 	object@sequences[[newRef]]  <- subseq(currentRef, 1, splitPoint * object@binSize)
# 	object@sequences[[prevRef]] <- subseq(currentRef, (splitPoint+1) * object@binSize, length(currentRef))
#     }
#     return(object)
# }

splitChromosomes <- function(object) {
    if (! is(object, "tenxcheckerClass")) {
        stop("Parameter should be a tenxcheckerClass.")
    }
    nSplits <- nrow(object@breaks)
    if (nSplits == 0) {
        return(invisible(object))
    }
    object <- splitCpp(object)
    for (i in seq_along(object@data)) {
        object@data[[i]]@interactionMatrix <- as_tibble(object@data[[i]]@interactionMatrix)
    }
    checkSizeDifference(object)
    checkAllBinDifference(object)
    checkMatrices(object)
#   newRefs            <- paste0("new_ref_", seq.int(nrow(object@breaks)))
#   object@chromosomes <- c(object@chromosomes, newRefs)
#   parameters <- object@breaks %>%
#       arrange(ref, desc(bin)) %>%
#       mutate(newRef = newRefs) %>%
#       transpose()
#   pb <- progress_bar$new(total = nSplits)
#   for (param in parameters) {
#       object <- splitChromosome(object, param)
#       pb$tick()
#   }
    gc(verbose = FALSE)
    return(invisible(object))
}
