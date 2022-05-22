###############################################################################
### msscafParameters S4 class definition
###############################################################################
#' Infrastructure for msscaf experiment and differential interaction
#'
#' \code{msscafParameters} is an S4 class providing the infrastructure (slots)
#' to store the input data, methods parameters, intermediate calculations
#' and results of a differential interaction pipeline
#'
#' @details \code{msscafParameters} does this and that...
#' TODO
#'
#' @name msscafExp
#' @rdname msscafExp
#' @docType class
#' @aliases msscafParameters msscafParameters-class
#'
#' @slot inputMatrix  The input matrix
#' @slot parameters   An named \code{list}. The parameters for the
#'                    segmentation methods. See \code{\link{parameters}}.
#'
#' @export
setClass("msscafParameters", slots = c(minCount        = "ANY",
                                       minRowCount     = "ANY",
                                       maxRowCount     = "ANY",
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
                                       nBinZoom        = "ANY",
                                       metaSize        = "ANY",
                                       rowCountSize    = "ANY",
                                       rowCountMu      = "ANY")
)


###############################################################################
### msscafData S4 class definition
###############################################################################
#' Infrastructure for msscaf data
#'
#' \code{msscafData} is an S4 class providing the ...
#'
#' @details \code{msscafData} does this and that...
#' TODO
#'
#' @name msscafData
#' @rdname msscafData
#' @docType class
#' @aliases msscafData msscafData-class
#'
#' @slot inputMatrix  The input matrix
#' @slot binSize      The resolution
#' @slot maxLinkRange The maximum size of a molecule
#' @slot sizes        The reference sizes
#'
#' @export
setClass("msscafData", slots = c(inputMatrix  = "ANY",
                                 maxLinkRange = "ANY") 
)


##- msscafData S4 class constructor --------------------------------------#
##----------------------------------------------------------------------------#
#' @rdname msscafData
#' @docType class
#'
#' @param inputMatrix A matrix with the data.
#'
#' @return \code{msscafExp} constructor returns an \code{msscafExp}
#'         object of class S4.
#'
#' @examples
#'
#' @export
msscafData <- function(inputMatrix  = NULL,
                       maxLinkRange = NULL) {
    
    ##- checking general input arguments -------------------------------------#
    ##------------------------------------------------------------------------#
    
    ##- dataSet
    if (is.null(inputMatrix)) {
        stop("'inputMatrix' must be specified", call. = FALSE)
    }
    if (!tibble::is_tibble(inputMatrix)) {
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
    
    object <- new("msscafData")
    
    object@inputMatrix  <- inputMatrix
    object@maxLinkRange <- maxLinkRange

    chromosomes <- gtools::mixedsort(levels(object@inputMatrix$ref1))
    object@inputMatrix <- object@inputMatrix %>%
        dplyr::mutate(ref1 = forcats::fct_relevel(ref1, chromosomes)) %>%
        dplyr::mutate(ref2 = forcats::fct_relevel(ref2, chromosomes))
    
    return(invisible(object))
}

setClass("msscafBreaks", slots = c(data            = "ANY",
                                   filteredData    = "ANY",
                                   changePlots     = "ANY",
                                   changeDistPlots = "ANY",
                                   mapPlots        = "ANY")
)

setClass("msscafJoins", slots = c(data      = "ANY",
                                  testPlots = "ANY",
                                  mapPlots  = "ANY")
)


###############################################################################
### msscaf S4 class definition
###############################################################################
#' Infrastructure for msscaf experiment and differential interaction
#'
#' \code{msscaf} is an S4 class providing the infrastructure (slots)
#' to store the input data, methods parameters, intermediate calculations
#' and results of a differential interaction pipeline
#'
#' @details \code{msscaf} does this and that...
#' TODO
#'
#' @name msscafExp
#' @rdname msscafExp
#' @docType class
#' @aliases msscafExp msscafExp-class
#'
#' @slot data  The input matrix
#' @slot parameters   An named \code{list}. The parameters for the
#'                    segmentation methods. See \code{\link{parameters}}.
#'
#' @export
setClass("msscafExp", slots = c(interactionMatrix = "ANY",
                                name              = "ANY",
                                parameters        = "ANY",
                                breaks            = "ANY",
                                outlierBins       = "ANY",
                                lowCountRows      = "ANY",
                                joins             = "ANY")
)


###############################################################################
### msscaf S4 class definition
###############################################################################
#' Infrastructure for msscaf experiment and differential interaction
#'
#' \code{msscaf} is an S4 class providing the infrastructure (slots)
#' to store the input data, methods parameters, intermediate calculations
#' and results of a differential interaction pipeline
#'
#' @details \code{msscaf} does this and that...
#' TODO
#'
#' @name msscafExp
#' @rdname msscafExp
#' @docType class
#' @aliases msscafExp msscafExp-class
#'
#' @slot data  The input matrix
#' @slot parameters   An named \code{list}. The parameters for the
#'                    segmentation methods. See \code{\link{parameters}}.
#'
#' @export
setClass("msscafClass", slots = c(data               = "ANY",
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


##- msscafExp S4 class constructor --------------------------------------#
##----------------------------------------------------------------------------#
#' @rdname msscafExp
#' @docType class
#'
#' @param inputMatrix A matrix with the data.
#'
#' @return \code{msscafExp} constructor returns an \code{msscafExp}
#'         object of class S4.
#'
#' @examples
#'
#' @export
msscaf <- function(sequenceFileName, binSize, minNBins = 20) {
    
    ##- checking general input arguments -------------------------------------#
    ##------------------------------------------------------------------------#
    if (!is.character(sequenceFileName)) {
        stop("The sequence file should be a character.")
    }
    if (!is.numeric(binSize)) {
        stop("'bin size' should be numeric.")
    }
    
    ##- end checking ---------------------------------------------------------#
    
    object                   <- new("msscafClass")
    object@binSize           <- binSize
    object@minNBins          <- minNBins
    object@data              <- c()
    object@breaks            <- NULL
    object@breakPlots        <- NULL
    object@joins             <- NULL
    object@joinPlots         <- NULL
    object@newChromosomes    <- c()
    object@mergedChromosomes <- c()

    object@sequences         <- Biostrings::readDNAStringSet(sequenceFileName) %>% as.character()
    names(object@sequences)  <- unlist(purrr::map(stringr::str_split(names(object@sequences), " "), 1))
    object@chromosomes       <- gtools::mixedsort(names(object@sequences))
    object@sizes             <- ceiling(nchar(object@sequences) / object@binSize) - 1
    names(object@sizes)      <- object@chromosomes
    return(invisible(object))
}

##- add experiment -----------------------------------------------------------#
##----------------------------------------------------------------------------#
#' @param data A \code{msscafData}
#'
#' @return An \code{msscafExp}
#'
#' @examples
#'
#' @export
addExp <- function(object, data, expName) {
    
    ##- checking general input arguments -------------------------------------#
    ##------------------------------------------------------------------------#
    if (! is(object, "msscafClass")) {
        stop("Input should be 'msscafClass'.")
    }
    if (! is(data, "msscafData")) {
        stop("Input should be 'msscafData'.")
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

    parameters <- new("msscafParameters")
    parameters@binSize         <- object@binSize
    parameters@maxLinkRange    <- data@maxLinkRange
    parameters@sampleSize      <- 10000
    parameters@loessSpan       <- 0.5
    parameters@minCount        <- 3
    parameters@minRowCount     <- 100
    parameters@maxRowCount     <- 1000000
    parameters@breakThreshold  <- NULL
    parameters@breakNCells     <- NULL
    parameters@distanceCount   <- NULL
    parameters@cornerScores    <- NULL
    parameters@cornerLimit     <- NULL
    parameters@nRandomizations <- 10
    parameters@nBinZoom        <- 100
    parameters@metaSize        <- 1

    newData <- new("msscafExp")
    newData@interactionMatrix <- data@inputMatrix
    newData@name              <- expName
    newData@parameters        <- parameters
    newData@breaks            <- NULL
    newData@joins             <- NULL
    newData@outlierBins       <- tibble::tibble(ref = factor(), bin = integer())

    # Check whether some bins are out of range
    overSized <- checkBinDifference(newData, object@sizes)
    if (nrow(overSized) > 0) {
        warning(paste("Some bins are beyond the reference sizes:",
            str(overSized),
            "The reference sequences and bins may be uncompatible.",
            "Removing the bins."))
        newData@interactionMatrix <- newData@interactionMatrix %>%
            dplyr::mutate(size1 = object@sizes[ref1]) %>%
            dplyr::mutate(size2 = object@sizes[ref2]) %>%
            dplyr::filter(bin1 <= size1, bin2 <= size2) %>%
            dplyr::select(-c(size1, size2))
    }

    object@data <- c(object@data, newData)

    return(invisible(object))
}

###############################################################################
### msscafRef S4 class definition
###############################################################################
#' Infrastructure for 1 reference
#'
#' \code{msscaf} is an S4 class providing the infrastructure (slots)
#' to store the input data, methods parameters, intermediate calculations
#' and results of a differential interaction pipeline
#'
#' @details \code{msscaf} does this and that...
#' TODO
#'
#' @name msscafExp
#' @rdname msscafExp
#' @docType class
#' @aliases msscafRefExp msscafRefExp-class
#'
#' @slot inputMatrix  The input matrix
#' @slot parameters   An named \code{list}. The parameters for the
#'                    segmentation methods. See \code{\link{parameters}}.
#'
#' @export
setClass("msscafRefExp", slots = c(interactionMatrix  = "ANY",
                                   chromosome         = "ANY",
                                   size               = "ANY",
                                   name               = "ANY",
                                   outlierBins        = "ANY",
                                   parameters         = "ANY")
)


##- msscaf S4 class constructor -----------------------------------------#
##----------------------------------------------------------------------------#
#' @rdname msscaf
#' @docType class
#'
#' @param inputMatrix A matrix with the data.
#'
#' @return \code{msscaf} constructor returns an \code{msscafExp}
#'         object of class S4.
#'
#' @examples
#'
#' @export
msscafRefExp <- function(matrix      = NULL,
                         chromosome  = NULL,
                         size        = NULL,
                         name        = NULL,
                         outlierBins = NULL,
                         parameters  = NULL) {
    
    ##- checking general input arguments -------------------------------------#
    ##------------------------------------------------------------------------#
    
    ##- dataSet
    if (is.null(matrix)) {
        stop("'matrix' must be specified", call. = FALSE)
    }
    if (!tibble::is_tibble(matrix)) {
        stop("'matrix' should be a tibble", call. = FALSE)
    }
    if (is.null(chromosome) | is.na(chromosome)) {
        stop("'chromosome' must be specified", call. = FALSE)
    }
    if (is.null(size) | is.na(size)) {
        stop(paste0("'size' must be specified for ref ", chromosome), call. = FALSE)
    }
    if (is.null(name) | is.na(name)) {
        stop("'name' must be specified", call. = FALSE)
    }
    if (! tibble::is_tibble(outlierBins)) {
        stop("'outlierBins' should be a tibble", call. = FALSE)
    }
    if (!is(parameters, "msscafParameters")) {
        stop("'parameters' should be 'msscafParameters' class instance", call. = FALSE)
    }
    
    ##- parameters
    
    
    ##- end checking ---------------------------------------------------------#
    
    object <- new("msscafRefExp")

    object@chromosome        <- chromosome
    object@size              <- size
    object@interactionMatrix <- matrix
    object@name              <- name
    object@outlierBins       <- outlierBins
    object@parameters        <- parameters
    
    return(invisible(object))
}


###############################################################################
### msscaf2Ref S4 class definition
###############################################################################
#' Infrastructure for 1 reference
#'
#' \code{msscaf} is an S4 class providing the infrastructure (slots)
#' to store the input data, methods parameters, intermediate calculations
#' and results of a differential interaction pipeline
#'
#' @details \code{msscaf} does this and that...
#' TODO
#'
#' @name msscafExp
#' @rdname msscafExp
#' @docType class
#' @aliases msscafRefExp msscafRefExp-class
#'
#' @slot inputMatrix  The input matrix
#' @slot parameters   An named \code{list}. The parameters for the
#'                    segmentation methods. See \code{\link{parameters}}.
#'
#' @export
setClass("msscaf2RefExp", slots = c(interactionMatrix  = "ANY",
                                    chromosome1        = "ANY",
                                    chromosome2        = "ANY",
                                    size1              = "ANY",
                                    size2              = "ANY",
                                    name               = "ANY",
                                    parameters         = "ANY")
)


##- msscafExp S4 class constructor --------------------------------------#
##----------------------------------------------------------------------------#
#' @rdname msscafExp
#' @docType class
#'
#' @param inputMatrix A matrix with the data.
#'
#' @return \code{msscafExp} constructor returns an \code{msscafExp}
#'         object of class S4.
#'
#' @examples
#'
#' @export
msscaf2RefExp <- function(matrix      = NULL,
                          chromosome1 = NULL,
                          chromosome2 = NULL,
                          size1       = NULL,
                          size2       = NULL,
                          name        = NULL,
                          parameters  = NULL) {
    
    ##- checking general input arguments -------------------------------------#
    ##------------------------------------------------------------------------#
    
    ##- dataSet
    if (is.null(matrix)) {
        stop("'matrix' must be specified", call. = FALSE)
    }
    if (!tibble::is_tibble(matrix)) {
        stop("'matrix' should be a tibble", call. = FALSE)
    }
    if (is.null(name) | is.na(name)) {
        stop("'name' must be specified", call. = FALSE)
    }
    
    ##- parameters
    
    
    ##- end checking ---------------------------------------------------------#
    
    object <- new("msscaf2RefExp")
    
    object@chromosome1       <- chromosome1
    object@chromosome2       <- chromosome2
    object@size1             <- size1
    object@size2             <- size2
    object@interactionMatrix <- matrix
    object@name              <- name
    object@parameters        <- parameters
    
    return(invisible(object))
}
