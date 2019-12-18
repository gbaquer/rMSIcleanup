#' rMSIcleanup
#'
#' Anotation and removal of the matrix peaks in MALDI data
#'
#' @docType package
#' @author Gerard Baquer Gómez
#' @exportPattern "^[[:alpha:]]+"
#' @import R.utils rMSIproc gridExtra  ggrepel graphics pack ggplot2 reshape2 enviPat PRROC Rtsne
#' @importFrom rMSI loadImageSliceFromMass loadImageSliceFromCols
#' @importFrom stats cor cov kmeans prcomp sd weighted.mean dist
#' @importFrom grDevices dev.new dev.off pdf
#' @importFrom utils write.table read.table data packageVersion
#' @importFrom stringr str_pad
#' @name rMSIcleanup
NULL
