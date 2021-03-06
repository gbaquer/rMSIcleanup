% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_module.R
\name{is_within_scan_tol}
\alias{is_within_scan_tol}
\title{Quiet}
\usage{
is_within_scan_tol(
  experimental_mass,
  calculated_mass,
  full_spectrum_masses,
  tol_scans
)
}
\arguments{
\item{experimental_mass}{Experimental mass}

\item{calculated_mass}{Calculated mass}

\item{full_spectrum_masses}{Full spectrum mass vector}

\item{tol_scans}{Tolerance specified in scans}
}
\value{
Boolean value indicating whether the experimental mass is within "tol_scans" of the calculated_mass
}
\description{
Forces a function to execute without printing any messages
}
