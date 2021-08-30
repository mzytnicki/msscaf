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
setClass("tenxcheckerParameters", slots = c(minCount        = "ANY",
                                            minRowCount     = "ANY",
                                            binSize         = "ANY",
                                            sampleSize      = "ANY",
                                            loessSpan       = "ANY",
                                            distanceCount   = "ANY",
                                            cornerScores    = "ANY",
                                            cornerLimit     = "ANY",
                                            breakThreshold  = "ANY", 
                                            breakNCells     = "ANY", 
                                            maxLinkRange    = "ANY", 
                                            nRandomizations = "ANY", 
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
                                      maxLinkRange = "ANY") 
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
                            maxLinkRange = NULL) {
    
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
    correctMatrixColNames <- c("ref1", "bin1", "ref2", "bin2", "count")
    if (any(colnames(inputMatrix) != correctMatrixColNames)) {
        stop(paste0("'inputMatrix' names are incorrect.  Should be ", paste0(correctMatrixColNames, collapse = ", "), ", found ", paste0(colnames(inputMatrix), collapse = ", "), "."), call. = FALSE)
    }
    if ((! is.factor(inputMatrix$ref1)) | (! is.factor(inputMatrix$ref2))) {
        stop("References should be factors.", call. = FALSE)
    }
    
    ##- parameters
    
    
    ##- end checking ---------------------------------------------------------#

    if (! identical(levels(inputMatrix$ref1), levels(inputMatrix$ref2))) {
        stop("Levels of the references should be identical.")
    }

    # Convert "diagonal matrices" to upper matrices
    diagonalData <- inputMatrix %>%
        dplyr::filter(ref1 == ref2) %>%
        dplyr::mutate(maxBin = pmax(bin1, bin2)) %>%
        dplyr::mutate(minBin = pmin(bin1, bin2)) %>%
        dplyr::select(-c(bin1, bin2))            %>%
        dplyr::rename(bin1 = maxBin)             %>%
        dplyr::rename(bin2 = minBin)
    inputMatrix <- inputMatrix %>% 
        dplyr::filter(ref1 != ref2) %>%
        dplyr::bind_rows(diagonalData) %>%
        dplyr::arrange(ref1, ref2, bin1, bin2)
    
    object <- new("tenxcheckerData")
    
    object@inputMatrix  <- inputMatrix
    object@maxLinkRange <- maxLinkRange

    chromosomes <- mixedsort(levels(object@inputMatrix$ref1))
    object@inputMatrix <- object@inputMatrix %>%
        mutate(ref1 = fct_relevel(ref1, chromosomes)) %>%
        mutate(ref2 = fct_relevel(ref2, chromosomes))
    
    return(invisible(object))
}

setClass("tenxcheckerBreaks", slots = c(data            = "ANY",
                                        filteredData    = "ANY",
                                        changePlots     = "ANY",
                                        changeDistPlots = "ANY",
                                        mapPlots        = "ANY")
)

setClass("tenxcheckerJoins", slots = c(data      = "ANY",
                                       testPlots = "ANY",
                                       mapPlots  = "ANY")
)


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
setClass("tenxcheckerExp", slots = c(interactionMatrix = "ANY",
                                     name              = "ANY",
                                     parameters        = "ANY",
                                     breaks            = "ANY",
                                     joins             = "ANY")
)


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
setClass("tenxcheckerClass", slots = c(data               = "ANY",
                                       sequences          = "ANY",
                                       binSize            = "ANY",
                                       minNBins           = "ANY",
                                       chromosomes        = "ANY",
                                       sizes              = "ANY",
                                       breaks             = "ANY",
                                       breakPlots         = "ANY",
                                       joins              = "ANY",
                                       joinPlots          = "ANY",
                                       newChromosomes     = "ANY",
                                       mergedChromosomes  = "ANY")
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
tenxchecker <- function(sequenceFileName, binSize, minNBins = 20) {
    
    ##- checking general input arguments -------------------------------------#
    ##------------------------------------------------------------------------#
    if (!is.character(sequenceFileName)) {
        stop("The sequence file should be a character.")
    }
    if (!is.numeric(binSize)) {
        stop("'bin size' should be numeric.")
    }
    
    ##- end checking ---------------------------------------------------------#
    
    object                   <- new("tenxcheckerClass")
    object@binSize           <- binSize
    object@minNBins          <- minNBins
    object@data              <- c()
    object@breaks            <- NULL
    object@breakPlots        <- NULL
    object@joins             <- NULL
    object@joinPlots         <- NULL
    object@newChromosomes    <- c()
    object@mergedChromosomes <- c()

    object@sequences         <- readDNAStringSet(sequenceFileName) %>% as.character()
    names(object@sequences)  <- unlist(map(str_split(names(object@sequences), " "), 1))
    object@chromosomes       <- mixedsort(names(object@sequences))
    object@sizes             <- ceiling(nchar(object@sequences) / object@binSize) - 1
    names(object@sizes)      <- object@chromosomes
    return(invisible(object))
}

##- add experiment -----------------------------------------------------------#
##----------------------------------------------------------------------------#
#' @param data A \code{tenxcheckerData}
#'
#' @return An \code{tenxcheckerExp}
#'
#' @examples
#'
#' @export
addExp <- function(object, data, expName) {
    
    ##- checking general input arguments -------------------------------------#
    ##------------------------------------------------------------------------#
    if (! is(object, "tenxcheckerClass")) {
        stop("Input should be 'tenxcheckerClass'.")
    }
    if (! is(data, "tenxcheckerData")) {
        stop("Input should be 'tenxcheckerData'.")
    }
    if (! is.character(expName)) {
        stop("Name should be character.")
    }
    if (length(expName) != 1) {
        stop("Name should have size 1.")
    }
    if (expName %in% names(object@data)) {
        stop("Name is already used.  Please choose another one.")
    }
    ##- end checking ---------------------------------------------------------#


    data <- updateRefs(object, data)

    parameters <- new("tenxcheckerParameters")
    parameters@binSize         <- object@binSize
    parameters@maxLinkRange    <- data@maxLinkRange
    parameters@sampleSize      <- 10000
    parameters@loessSpan       <- 0.5
    parameters@minCount        <- 3
    parameters@minRowCount     <- 100
    parameters@breakThreshold  <- NULL
    parameters@breakNCells     <- NULL
    parameters@distanceCount   <- NULL
    parameters@cornerScores    <- NULL
    parameters@cornerLimit     <- NULL
    parameters@nRandomizations <- 10
    parameters@nBinZoom        <- 100

    newData <- new("tenxcheckerExp")
    newData@interactionMatrix <- data@inputMatrix
    newData@name              <- expName
    newData@parameters        <- parameters
    newData@breaks            <- NULL
    newData@joins             <- NULL

    object@data <- c(object@data, newData)

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
    if (is.null(chromosome) | is.na(chromosome)) {
        stop("'chromosome' must be specified", call. = FALSE)
    }
    if (is.null(size) | is.na(size)) {
        stop(paste0("'size' must be specified for ref ", chromosome), call. = FALSE)
    }
    if (!is(parameters, "tenxcheckerParameters")) {
        stop("'parameters' should be 'tenxcheckerParameters' class instance", call. = FALSE)
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
