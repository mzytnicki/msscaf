###############################################################################
#
# msscaf-package organization of R files
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
# msscaf
# |- functions
# --------------------------------
#
# --------------------------------
# msscaf
# |- msscafcore
#    |- page
#       |- function
# --------------------------------
#
###############################################################################

#' Title...
#'
#' \code{msscaf} is a package...
#' TODO

#'
#' @examples
#'
#' @docType package
#' @name msscaf
#' @aliases msscaf-package
#'
#' @author Matthias Zytnicki
#'
#' @import MASS
#' @import methods
#' @importFrom gtools mixedsort
#' @importFrom magrittr %>%
#' @import tidyverse
#' @import ggplot2
#' @import ggforce
#' @import GenomicRanges
#' @import BiocParallel
#' @import progress
#' @import Matrix
#' @import Biostrings
#' @importFrom cowplot plot_grid
#' @importFrom BiocManager version
#' @importFrom Rcpp evalCpp sourceCpp
#' @useDynLib msscaf
#'
#' @keywords package
NULL
