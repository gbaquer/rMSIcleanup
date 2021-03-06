% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calssification_module.R
\name{removeMatrix_kMeansTranspose}
\alias{removeMatrix_kMeansTranspose}
\title{Remove Matrix k-Means Transpose}
\usage{
removeMatrix_kMeansTranspose(
  pks_Norharmane,
  correlation = FALSE,
  normalize = TRUE
)
}
\arguments{
\item{pks_Norharmane}{Peak Matrix}

\item{correlation}{Binary valiable determining whether to use the correlation of the raw data or the raw data}

\item{normalize}{Binary variable determining whether to normalize the data or not. TIC normalization is used.}
}
\value{
List of mass indices considered to be endogenous. The rest of the peaks are deamed as matrix related or non-anatomically relevant.
}
\description{
Remove the matrix by performing k-means clustering on the transpose of the peak matrix.
The algorithm identifies similar spectral peaks.
}
