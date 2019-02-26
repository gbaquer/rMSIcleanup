#' rMSIcleanup
#'
#' Anotation and removal of the matrix peaks in MALDI data
#'
#' @docType package
#' @author Gerard Baquer Gómez
#' @exportPattern "^[[:alpha:]]+"
#' @import Rcpp R.utils rMSIproc gridExtra lattice RColorBrewer graphics GlobalOptions XML pack reticulate ggplot2 reshape2 enviPat tsne
#' @importFrom rMSI loadImageSliceFromMass loadImageSliceFromCols
#' @importFrom Rcpp evalCpp
#' @importFrom stats cor cov kmeans prcomp sd weighted.mean dist
#' @importFrom grDevices dev.new dev.off pdf rainbow palette
#' @importFrom utils write.table read.table data packageVersion
#' @importFrom stringr str_pad
#' @useDynLib rMSIcleanup, .registration = TRUE
#' @name rMSIcleanup
NULL
