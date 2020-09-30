###############################################################################
### tenxcheckerParameters S4 class definition
###############################################################################
#' Infrastructure for tenxchecker experiment and differential interaction
#'
#' \code{tenxcheckerParameters} is an S4 class providing the infrastructure (slots)
#' to store the input data, methods parameters, intermediate calculations
#' and results of a differential interaction pipeline
#'
#' @details \code{tenxcheckerParameters} does this and that...
#' TODO
#'
#' @name tenxcheckerExp
#' @rdname tenxcheckerExp
#' @docType class
#' @aliases tenxcheckerParameters tenxcheckerParameters-class
#'
#' @slot inputMatrix  The input matrix
#' @slot parameters   An named \code{list}. The parameters for the
#'                    segmentation methods. See \code{\link{parameters}}.
#'
#' @export
setClass("tenxcheckerParameters", slots = c(minNBins        = "ANY",
                                            minCount        = "ANY",
                                            minRowCount     = "ANY",
                                            binSize         = "ANY",
                                            sampleSize      = "ANY",
                                            loessSpan       = "ANY",
                                            breakThreshold  = "ANY", 
                                            breakNCells     = "ANY", 
                                            maxLinkRange    = "ANY", 
                                            nRandomizations = "ANY", 
                                            pvalueThreshold = "ANY", 
                                            nBinZoom        = "ANY")
)


###############################################################################
### tenxcheckerData S4 class definition
###############################################################################
#' Infrastructure for tenxchecker data
#'
#' \code{tenxcheckerData} is an S4 class providing the ...
#'
#' @details \code{tenxcheckerData} does this and that...
#' TODO
#'
#' @name tenxcheckerData
#' @rdname tenxcheckerData
#' @docType class
#' @aliases tenxcheckerData tenxcheckerData-class
#'
#' @slot inputMatrix  The input matrix
#' @slot binSize      The resolution
#' @slot maxLinkRange The maximum size of a molecule
#' @slot sizes        The reference sizes
#'
#' @export
setClass("tenxcheckerData", slots = c(inputMatrix  = "ANY",
                                      binSize      = "ANY",
                                      maxLinkRange = "ANY", 
                                      sizes        = "ANY")
)


##- tenxcheckerData S4 class constructor --------------------------------------#
##----------------------------------------------------------------------------#
#' @rdname tenxcheckerData
#' @docType class
#'
#' @param inputMatrix A matrix with the data.
#'
#' @return \code{tenxcheckerExp} constructor returns an \code{tenxcheckerExp}
#'         object of class S4.
#'
#' @examples
#'
#' @export
tenxcheckerData <- function(inputMatrix  = NULL,
                            binSize      = NULL,
                            maxLinkRange = NULL,
                            sizes        = NULL) {
    
    ##- checking general input arguments -------------------------------------#
    ##------------------------------------------------------------------------#
    
    ##- dataSet
    if (is.null(inputMatrix)) {
        stop("'inputMatrix' must be specified", call. = FALSE)
    }
    if (!is_tibble(inputMatrix)) {
        stop("'inputMatrix' should be a tibble", call. = FALSE)
    }
    if (ncol(inputMatrix) != 5) {
        stop("'inputMatrix' should have 5 columns", call. = FALSE)
    }
    if (any(colnames(inputMatrix) != c("ref1", "bin1", "ref2", "bin2", "count"))) {
        stop("'inputMatrix' names are incorrect", call. = FALSE)
    }
    
    ##- parameters
    
    
    ##- end checking ---------------------------------------------------------#
    
    object <- new("tenxcheckerData")
    
    object@inputMatrix  <- inputMatrix
    object@binSize      <- binSize
    object@maxLinkRange <- maxLinkRange
    object@sizes        <- sizes
    
    return(invisible(object))
}

###############################################################################
### tenxchecker S4 class definition
###############################################################################
#' Infrastructure for tenxchecker experiment and differential interaction
#'
#' \code{tenxchecker} is an S4 class providing the infrastructure (slots)
#' to store the input data, methods parameters, intermediate calculations
#' and results of a differential interaction pipeline
#'
#' @details \code{tenxchecker} does this and that...
#' TODO
#'
#' @name tenxcheckerExp
#' @rdname tenxcheckerExp
#' @docType class
#' @aliases tenxcheckerExp tenxcheckerExp-class
#'
#' @slot data  The input matrix
#' @slot parameters   An named \code{list}. The parameters for the
#'                    segmentation methods. See \code{\link{parameters}}.
#'
#' @export
setClass("tenxcheckerExp", slots = c(interactionMatrix  = "ANY",
                                     chromosomes        = "ANY",
                                     sizes              = "ANY",
                                     lowCounts          = "ANY",
                                     newChromosomes     = "ANY",
                                     mergedChromosomes  = "ANY",
                                     parameters         = "ANY")
)


##- tenxcheckerExp S4 class constructor --------------------------------------#
##----------------------------------------------------------------------------#
#' @rdname tenxcheckerExp
#' @docType class
#'
#' @param inputMatrix A matrix with the data.
#'
#' @return \code{tenxcheckerExp} constructor returns an \code{tenxcheckerExp}
#'         object of class S4.
#'
#' @examples
#'
#' @export
tenxcheckerExp <- function(data = NULL) {
    
    ##- checking general input arguments -------------------------------------#
    ##------------------------------------------------------------------------#
    
    ##- dataSet
    
    ##- parameters
    
    
    ##- end checking ---------------------------------------------------------#
    
    object <- new("tenxcheckerExp")
    object@parameters <- new("tenxcheckerParameters")
    object@interactionMatrix       <- data@inputMatrix
    object@parameters@binSize      <- data@binSize
    object@parameters@maxLinkRange <- data@maxLinkRange
    object@sizes                   <- data@sizes
    
    object@chromosomes <- mixedsort(unique(c(as.vector(object@interactionMatrix$ref1), as.vector(object@interactionMatrix$ref2))))
    object@interactionMatrix <- object@interactionMatrix %>%
        mutate(ref1 = factor(ref1, levels = object@chromosomes)) %>%
        mutate(ref2 = factor(ref2, levels = object@chromosomes))
    object@newChromosomes    <- c()
    object@mergedChromosomes <- c()
    
    object@parameters@sampleSize      <- 10000
    object@parameters@loessSpan       <- 0.5
    object@parameters@minNBins        <- 20
    object@parameters@minCount        <- 3
    object@parameters@minRowCount     <- 100
    object@parameters@breakThreshold  <- -1
    object@parameters@breakNCells     <- 100
    object@parameters@nRandomizations <- 10
    object@parameters@pvalueThreshold <- 1e-4
    object@parameters@nBinZoom        <- 100

    if (! is.null(object@sizes)) {
        object@sizes <- computeRefSizes(object)
    }
    
    return(invisible(object))
}

###############################################################################
### tenxcheckerRef S4 class definition
###############################################################################
#' Infrastructure for 1 reference
#'
#' \code{tenxchecker} is an S4 class providing the infrastructure (slots)
#' to store the input data, methods parameters, intermediate calculations
#' and results of a differential interaction pipeline
#'
#' @details \code{tenxchecker} does this and that...
#' TODO
#'
#' @name tenxcheckerExp
#' @rdname tenxcheckerExp
#' @docType class
#' @aliases tenxcheckerRefExp tenxcheckerRefExp-class
#'
#' @slot inputMatrix  The input matrix
#' @slot parameters   An named \code{list}. The parameters for the
#'                    segmentation methods. See \code{\link{parameters}}.
#'
#' @export
setClass("tenxcheckerRefExp", slots = c(interactionMatrix  = "ANY",
                                        chromosome         = "ANY",
                                        size               = "ANY",
                                        parameters         = "ANY")
)


##- tenxchecker S4 class constructor -----------------------------------------#
##----------------------------------------------------------------------------#
#' @rdname tenxchecker
#' @docType class
#'
#' @param inputMatrix A matrix with the data.
#'
#' @return \code{tenxchecker} constructor returns an \code{tenxcheckerExp}
#'         object of class S4.
#'
#' @examples
#'
#' @export
tenxcheckerRefExp <- function(matrix     = NULL,
                              chromosome = NULL,
                              size       = NULL,
                              parameters = NULL) {
    
    ##- checking general input arguments -------------------------------------#
    ##------------------------------------------------------------------------#
    
    ##- dataSet
    if (is.null(matrix)) {
        stop("'matrix' must be specified", call. = FALSE)
    }
    if (!is_tibble(matrix)) {
        stop("'matrix' should be a tibble", call. = FALSE)
    }
    
    ##- parameters
    
    
    ##- end checking ---------------------------------------------------------#
    
    object <- new("tenxcheckerRefExp")
    
    object@chromosome        <- chromosome
    object@size              <- size
    object@interactionMatrix <- matrix
    object@parameters        <- parameters
    
    return(invisible(object))
}


###############################################################################
### tenxchecker2Ref S4 class definition
###############################################################################
#' Infrastructure for 1 reference
#'
#' \code{tenxchecker} is an S4 class providing the infrastructure (slots)
#' to store the input data, methods parameters, intermediate calculations
#' and results of a differential interaction pipeline
#'
#' @details \code{tenxchecker} does this and that...
#' TODO
#'
#' @name tenxcheckerExp
#' @rdname tenxcheckerExp
#' @docType class
#' @aliases tenxcheckerRefExp tenxcheckerRefExp-class
#'
#' @slot inputMatrix  The input matrix
#' @slot parameters   An named \code{list}. The parameters for the
#'                    segmentation methods. See \code{\link{parameters}}.
#'
#' @export
setClass("tenxchecker2RefExp", slots = c(interactionMatrix  = "ANY",
                                         chromosome1        = "ANY",
                                         chromosome2        = "ANY",
                                         size1              = "ANY",
                                         size2              = "ANY",
                                         parameters         = "ANY")
)


##- tenxcheckerExp S4 class constructor --------------------------------------#
##----------------------------------------------------------------------------#
#' @rdname tenxcheckerExp
#' @docType class
#'
#' @param inputMatrix A matrix with the data.
#'
#' @return \code{tenxcheckerExp} constructor returns an \code{tenxcheckerExp}
#'         object of class S4.
#'
#' @examples
#'
#' @export
tenxchecker2RefExp <- function(matrix      = NULL,
                               chromosome1 = NULL,
                               chromosome2 = NULL,
                               size1       = NULL,
                               size2       = NULL,
                               parameters  = NULL) {
    
    ##- checking general input arguments -------------------------------------#
    ##------------------------------------------------------------------------#
    
    ##- dataSet
    if (is.null(matrix)) {
        stop("'matrix' must be specified", call. = FALSE)
    }
    if (!is_tibble(matrix)) {
        stop("'matrix' should be a tibble", call. = FALSE)
    }
    
    ##- parameters
    
    
    ##- end checking ---------------------------------------------------------#
    
    object <- new("tenxchecker2RefExp")
    
    object@chromosome1       <- chromosome1
    object@chromosome2       <- chromosome2
    object@size1             <- size1
    object@size2             <- size2
    object@interactionMatrix <- matrix
    object@parameters        <- parameters
    
    return(invisible(object))
}
