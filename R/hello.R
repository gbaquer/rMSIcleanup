# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
#   Load function:             'Ctrl + Shift + L'
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Run Roxygen2:              'Ctrl + Shift + D'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#LOAD LIBRARIES
library(R.utils)

#GLOBAL VARIABLES

#' Print verbose
#'
#' Print \code{message} if the verbose \code{level} is higher than the specified \code{threshold}.
#'
#'  1: Silent
#'  0: Application messages
#' -1: Code Sections
#' -2: Code subsections
#' -3: Debug mode
#'
#' @param message String to print out
#' @param level Verbose level
#'
#' @return None
#'
#' @examples
#' printv("Hello World!",0)
#'
#' @export
printv <- function(message,level) {
  threshold=0
  if(level>=threshold)
  {
    print(message)
  }
}

#START FUNCTIONS
#' Hello world
#'
#' Hello world function
#'
#' @return None
#'
#' @examples
#' hello()
#'
#' @export
hello <- function() {
  print("Hello, world 2!")
}


#' Identify matrix correlation
#'
#' Identifies de matrix peaks by following the method presented by Fonville et al. 2012.
#'
#' @return None
#'
#' @examples
#' identify_matrix_correlation()
#'
#' @export
identify_matrix_correlation <- function () {
  #LOAD IMAGE
  printv("2",2)
  printv("1",1)
  printv("0",0)
  printv("-1",-1)
  printv("-2",-2)

  #PROCESS IMAGE

  #GET SUBSETS (TISSUE & OUTSIDE)

  #DETERMINE CORRELATION MATRIX
}



