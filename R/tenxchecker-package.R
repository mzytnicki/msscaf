###############################################################################
#
# tenxchecker-package organization of R files
#
# AllClasses .... class definitions and object constructors
# AllGenerics ... the generics defined in srnadiff-package
# methods ....... the S4 methods (accessors for slots and replace methods)
# helper ........ computeNormFactors, computeLogFoldChange, cvgNormalization,
#                 reconcileRegions, checkParameters, IRlist2GR
# RcppExports ... the R wrappers for the C++ functions (auto)
#
#
# General outline of the internal function calls.
# Note: not all of these functions are exported.
#
# --------------------------------
# tenxchecker
# |- functions
# --------------------------------
#
# --------------------------------
# tenxchecker
# |- tenxcheckercore
#    |- page
#       |- function
# --------------------------------
#
###############################################################################

#' Title...
#'
#' \code{tenxchecker} is a package...
#' TODO

#'
#' @examples
#'
#' @docType package
#' @name tenxchecker
#' @aliases tenxchecker-package
#'
#' @author Matthias Zytnicki
#'
#' @import methods
#' @import gtools
#' @import magrittr
#' @import tidyverse
#' @import ggplot2
#' @import ggforce
#' @import GenomicRanges
#' @import BiocParallel
#' @import progress
#' @import Matrix
#' @import Biostrings
#' @import gridExtra
#' @importFrom BiocManager version
#' @importFrom Rcpp evalCpp sourceCpp
#' @useDynLib tenxchecker
#'
#' @keywords package
NULL
