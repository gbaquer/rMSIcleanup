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
#
#   Comment/Uncomment:         'Ctrl + Shift + C'

#LOAD LIBRARIES
library("R.utils")
library("rMSI")
library("rMSIproc")
library("lattice")
library("gridExtra")

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
#' find_matrix_correlation()
#'
#' @export
find_matrix_correlation <- function () {
  #SECTION 0 :: Load peak matrix
  pks_Norharmane <- rMSIproc::LoadPeakMatrix("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/comparativa_matrix_au/peak_matrix_norharmane/mergeddata-peaks.zip")
  pks_Au <- rMSIproc::LoadPeakMatrix("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/comparativa_matrix_au/peak_matrix_au/mergeddata-peaks.zip")

  #SECTION 1:: Peak Selection

  #Approach A: Use of matrix peaks to remove nonbiological variables
  #A.1. Load manually identified peaks outside of tissue
  exo_peaks <- c(753.7274,672.5297,678.5946)
  exo_is=double()
  for (p in exo_peaks)
    exo_is=append(exo_is,match(TRUE,pks_Norharmane$mass>=p))

  #A.2. Get peaks with highest correlation to the manually identified peaks (to find other likely exogenous peaks)

  #Compute and plot covariance and correlation
  covMAT=cov(pks_Norharmane$intensity)
  corMAT=cor(pks_Norharmane$intensity)
  grid.arrange(levelplot(covMAT))
  grid.arrange(levelplot(corMAT))

  #Rank correlation of exo indices [IMPROVE: get the "max_exo" elements that correlate the best to the three of them]
  max_exo=10
  top_exo_mat=integer()
  for (i in exo_is)
    top_exo_mat=cbind(top_exo_mat,sort(corMAT[i,],decreasing=TRUE,index.return=TRUE)$ix)
  top_exo=unique(as.vector(t(top_exo_mat)))[1:max_exo]

  #A.3. Choose likely biological peaks (the ones that correlate negatively)
  #Compute mean correlation to the top exo peaks
  mean_exo_cor=apply(corMAT[top_exo,],2,mean)
  bio_peaks=which(mean_exo_cor<0)
  nonbio_peaks=which(mean_exo_cor>=0)

  #A.4. Remove the peaks from the analysis (In the paper they just set the whole peak to 0, what else could be done?)
  pks_Norharmane_A=pks_Norharmane
  pks_Norharmane_A$mass=pks_Norharmane_A$mass[bio_peaks]
  pks_Norharmane_A$intensity=pks_Norharmane_A$intensity[,bio_peaks]
  pks_Norharmane_A$area=pks_Norharmane_A$area[,bio_peaks]
  pks_Norharmane_A$SNR=pks_Norharmane_A$SNR[,bio_peaks]

  #Approach B: Use Image Anatomy to Find relevant variables
  #B.1. variance explained (VE) in first singular value
  #Dividing the square of the first value of the diagonal matrix by the sum of all its squared elements
  mean_centered_intensity=scale(pks_Norharmane_A$intensity,scale=F)#pks_Norharmane_processed$intensity-mean(pks_Norharmane_processed$intensity)
  ve=double()
  svd_total=sum(svd(mean_centered_intensity,0,0)$d^2)^0.5
  for (i in 1:length(pks_Norharmane_A$mass))
    ve=append(ve,(svd(mean_centered_intensity[,i],0,0)$d)/svd_total)
  #B.2. VE-threshold for variable selection
  ve_threshold=0.05   #HOW TO DETERMINE THIS THRESHOLD heuristic parameter based on pareto-efficiency considerations (sum of all explained variances divided by the number of variables)
  anatomical_peaks=which(ve>ve_threshold)
  nonanatomical_peaks=which(ve<=ve_threshold)


  #B.4. Remove the peaks from the analysis (In the paper they just set the whole peak to 0, what else could be done?)
  pks_Norharmane_AB=pks_Norharmane_A
  pks_Norharmane_AB$mass=pks_Norharmane_AB$mass[anatomical_peaks]
  pks_Norharmane_AB$intensity=pks_Norharmane_AB$intensity[,anatomical_peaks]
  pks_Norharmane_AB$area=pks_Norharmane_AB$area[,anatomical_peaks]
  pks_Norharmane_AB$SNR=pks_Norharmane_AB$SNR[,anatomical_peaks]

  #SECTION 3: Pixel selection
  #Determine inner pixels
  TIC_AB=apply(pks_Norharmane_AB$intensity,1,sum)
  TIC=apply(pks_Norharmane$intensity,1,sum)
  TIC_ratio=log10(TIC_AB/TIC)
  pixels_out=which(TIC_ratio<(-0.5))
  pixels_in=which(TIC_ratio>=(-0.5))


  #Remove outside peaks
  pks_Norharmane_Final=pks_Norharmane_AB
  pks_Norharmane_Final$numPixels=length(pixels_in)
  pks_Norharmane_Final$intensity=pks_Norharmane_Final$intensity[pixels_in,]
  pks_Norharmane_Final$area=pks_Norharmane_Final$area[pixels_in,]
  pks_Norharmane_Final$SNR=pks_Norharmane_Final$SNR[pixels_in,]
  pks_Norharmane_Final$pos=pks_Norharmane_Final$pos[pixels_in,]
  pks_Norharmane_Final$posMotors=pks_Norharmane_Final$posMotors[pixels_in,]
  pks_Norharmane_Final$normalizations=pks_Norharmane_Final$normalizations[pixels_in,]

  #Smoothing
  #[To Be Implemented]

  #Print image
  rMSIproc::plotPeakImage(pks_Norharmane_Final,c=2)
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
  #printv("SECTION A: Loading and preprocessing image",0)
  #info_in <- getrMSIdataInfo(img_in)
  #lapply(img_in$data, function(x){ ff::close.ff(x) }) #close ff files
  # FullImageSmoothing(fileNames = info_in$filenames,
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
  img_out=img_in
  #img_out$name="Output Image"
  #center_xy=ceiling(img_out$size/2)
  #radius=min(center_xy)/2

  #SECTION C: REMOVE THE MATRIX PEAKS
  printv("SECTION C: Removing matrix peaks",0)

  #SECTION D: STORE CLEAN IMAGE
  printv("SECTION D: Storing clean image",0)
  #Plot 2 images
  rMSI::MSIWindow(img_in,img_out)
}



