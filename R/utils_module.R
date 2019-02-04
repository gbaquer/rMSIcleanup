#########################################################################
#
#     UTILS MODULE
#
#########################################################################
#     rMSIcleanup - R package for MSI matrix removal
#     Copyright (C) 2019 Gerard Baquer Gómez
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
############################################################################


#' Package Global Options
#'
#' Global options used by the package
#'
#' @section verbose_level:
#' Verbose level determining which level of information is printed or plotted
#' \itemize{
#'  \item  1: Silent
#'  \item  0: Application messages
#'  \item -1: Code Sections
#'  \item -2: Code subsections
#'  \item -3: Debug mode
#' }

#'
#' @param ... Options to change
#' @param RESET Reset the package options to the default
#' @param READ.ONLY Return only the READ.ONLY options
#' @param LOCAL Change to from global to local mode when TRUE and change from local to global mode when FALSE. In the local mode a copy of the options is created and the copy is modified.
#' @param ADD Add a new option on the fly

pkg_opt= set_opt(
  verbose_level=list(.value=1,
                     .read.only=FALSE,
                     .validate= function(x) is.numeric(x) && x%%1==0 && x>=-3 && x<=1,
                     .failed_msg = "Verbose should be an integer in the range [-3,1]",
                     .description = "Verbose level determining which level of information is printed or plotted"
  )
)

#' Specify decimal
#'
#' Exports the results to a text file which can be read by mmass for easy interpretation and validation of the results.
#'
#' @param x Input double
#' @param k Precision
#'
#' @return x with k precision
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))
