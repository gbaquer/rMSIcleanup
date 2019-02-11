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

#' Cosine similarity two vectors
#'
#' Returns the cosine similarity between two vectors
#'
#' @param ix Input vectors
#'
cos_sim_vectors <- function(ix)
{
  A = X[ix[1],]
  B = X[ix[2],]
  return( sum(A*B)/sqrt(sum(A^2)*sum(B^2)) )
}

#' Cosine similarity two vectors
#'
#' Returns the cosine similarity between two vectors
#'
#' @param X Input matrix
#'
cos_sim <- function(X)
{
  n <- nrow(X)
  cmb <- expand.grid(i=1:n, j=1:n)
  C <- matrix(apply(cmb,1,cos_sim_vectors),n,n)
  return(C)
}
#' Plot text
#'
#' Plots text with white background
#'
#' @param text Text to be plotted
#' @param size Size of the text
#'
plot_text <- function(text,size=8)
{
  p=ggplot() + annotate("text", x = 0, y = 0, size=size, label = text) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) + theme(axis.line=element_blank(),axis.text.x=element_blank(),
                                              axis.text.y=element_blank(),axis.ticks=element_blank(),
                                              axis.title.x=element_blank(),
                                              axis.title.y=element_blank(),legend.position="none",
                                              panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                                              panel.grid.minor=element_blank(),plot.background=element_blank())
  return(p)
}
#' Specify decimal
#'
#' Exports the results to a text file which can be read by mmass for easy interpretation and validation of the results.
#'
#' @param x Input double
#' @param k Precision
#'
#' @return x with k precision
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))
