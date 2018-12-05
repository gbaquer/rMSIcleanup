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
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#LOAD LIBRARIES
library(R.utils)
#Set Verbose threshold:
# 0: Application messages
# -1: Code Sections
# -2: Code subsections
# -3: Debug mode
v=Verbose(threshold=0)

#START FUNCTIONS
hello <- function() {
  print("Hello, world 2!")
}

#FONVILLE et al. 2012
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



