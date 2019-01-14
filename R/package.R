#' rMSIcleanup
#'
#' Anotation and removal of the matrix peaks in MALDI data
#'
#' @docType package
#' @author Gerard Baquer Gómez
#' @exportPattern "^[[:alpha:]]+"
#' @import Rcpp R.utils rMSIproc gridExtra lattice RColorBrewer graphics GlobalOptions XML base64enc
#' @importFrom Rcpp evalCpp
#' @importFrom stats cor cov kmeans prcomp sd
#' @importFrom grDevices dev.new palette
#' @importFrom utils write.table
#' @useDynLib rMSIcleanup, .registration = TRUE
#' @name rMSIcleanup
NULL
