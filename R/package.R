#' rMSIcleanup
#'
#' Anotation and removal of the matrix peaks in MALDI data
#'
#' @docType package
#' @author Gerard Baquer Gómez
#' @exportPattern "^[[:alpha:]]+"
#' @import Rcpp R.utils rMSI rMSIproc gridExtra lattice RColorBrewer
#' @importFrom Rcpp evalCpp
#' @importFrom stats cor cov kmeans prcomp sd
#' @importFrom graphics plot
#' @importFrom grDevices dev.new palette
#' @useDynLib rMSIcleanup, .registration = TRUE
#' @name rMSIcleanup
NULL
