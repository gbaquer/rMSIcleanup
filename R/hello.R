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

#' Verbose threshold:
#' 0: Application messages
#' -1: Code Sections
#' -2: Code subsections
#' -3: Debug mode
v=Verbose(threshold=0)

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
  printf(v,level=2,"2")
  printf(v,level=1,"1")
  printf(v,level=0,"0")
  printf(v,level=-1,"-1")
  printf(v,level=-2,"-2")

  #PROCESS IMAGE

  #GET SUBSETS (TISSUE & OUTSIDE)

  #DETERMINE CORRELATION MATRIX
}



