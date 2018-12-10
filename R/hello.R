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
library(rMSI)

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
  #SECTION A: LOAD & PREPROCESS IMAGE
  printv("SECTION A: Loading and preprocessing image",0)
  #img_in <- rMSI::LoadMsiData("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/comparativa_matrix_au/20170913_Norharmane_Fanny_RP_a.imzML")
  #info_in <- getrMSIdataInfo(img_in)
  #lapply(img_in$data, function(x){ ff::close.ff(x) }) #close ff files
  #FullImageSmoothing(fileNames = info_in$filenames,
  #                   massChannels = info_in$masschannels,
  #                   numRows = info_in$nrows,
  #                   dataType = info_in$datatype,
  #                   numOfThreads = parallel::detectCores(),
  #                   SmoothingKernelSize = 5)

  #SECTION B: DETERMINE MATRIX PEAKS
  printv("SECTION B: Determining matrix peaks",0)
  #B.1:Determine regions (inside and outside)
  #B.2:Determine correlation matrix
  #B.3:Determine where are the matrix peaks located

  #Dummy application draw a circle inside of the image
  #img_out=img_in
  #img_out$name="Output Image"
  #center_xy=ceiling(img_out$size/2)
  #radius=min(center_xy)/2

  #SECTION C: REMOVE THE MATRIX PEAKS
  printv("SECTION C: Removing matrix peaks",0)

  #SECTION D: STORE CLEAN IMAGE
  printv("SECTION D: Storing clean image",0)
  #Plot 2 images
  #rMSI::MSIWindow(img_in,img_out)
}



