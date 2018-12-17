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

#' Remove Matrix Padded
#'
#' Remove the matrix following the method presented by Fonville et al. 2012.
#' It takes a padded image where there is a wide enough region of pixels outside of the tissue.
#' It consists of three basic steps.
#'
#' @return None
#'
#' @examples
#' removeMatrix_padded()
#'
#' @export
removeMatrix_padded <- function () {
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

#' Remove Matrix Compare Au
#'
#' Remove the matrix by comparing anatomically similar regions of the tissue to a reference image taken
#' with a gold matrix. Assesses the degree of similarity between peaks
#' It takes a padded image where there is a wide enough region of pixels outside of the tissue.
#' It consists of three basic steps.
#'
#' @return None
#'
#' @examples
#' removeMatrix_compareAu()
#'
#' @export
removeMatrix_compareAu <- function () {
  #SECTION 0 :: Load peak matrix
  pks_Norharmane <- rMSIproc::LoadPeakMatrix("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/comparativa_matrix_au/peak_matrix_norharmane/mergeddata-peaks.zip")
  pks_Au <- rMSIproc::LoadPeakMatrix("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/comparativa_matrix_au/peak_matrix_au/mergeddata-peaks.zip")

  #SECTION 1:: Get regions of similarity
  regions <- list(
    num = 0,
    centers_Norharmane = t(matrix(c(130,54,61,60,111,92,121,13,14,71,119,30,67,46), nrow=2)),
    centers_Au =t(matrix(c(103,54,54,79,78,85,109,18,15,45,84,33,45,36), nrow=2)),
    radial_point_Norharmane =t(matrix(c(137,47,56,85,112,87,121,4,5,71,114,27,68,41), nrow=2)),
    radial_point_Au=t(matrix(c(106,63,56,83,76,79,110,9,8,41,80,33,43,30), nrow=2)),
    radius_Norharmane =matrix(),
    radius_Au =matrix(),
    pixels_Norharmane =list(),
    pixels_Au=list(),
    mean_Norharmane=list(),
    mean_Au=list()
  )

  regions$num=dim(regions$centers_Norharmane)[1]
  regions$radius_Norharmane=sqrt(apply((regions$centers_Norharmane - regions$radial_point_Norharmane) ^ 2,1,sum))
  regions$radius_Au=sqrt(apply((regions$centers_Au - regions$radial_point_Au) ^ 2,1,sum))

  for (i in 1:regions$nu)
  {
    distance_Norharmane = sqrt(apply((regions$centers_Norharmane[i,] - pks_Norharmane$pos) ^ 2,1,sum))
    regions$pixels_Norharmane[[paste("r",i,sep="_")]]=which(distance_Norharmane<regions$radius_Norharmane[i])

    distance_Au = sqrt(apply((regions$centers_Au[i,] - pks_Au$pos) ^ 2,1,sum))
    regions$pixels_Au[[paste("r",i,sep="_")]]=which(distance_Au<regions$radius_Au[i])
  }

  #SECTION 2:: Calculate, median and average spectrum
  for (i in 1:regions$num)
  {
    regions$mean_Norharmane[[i]]=apply(pks_Norharmane$intensity[regions$pixels_Norharmane[[i]],],2,mean)
    regions$mean_Au[[i]]=apply(pks_Au$intensity[regions$pixels_Au[[i]],],2,mean)
  }
  #In order to compare correlations they need to have the same size. Standardize and align spectra.
  mass_threshold=0.2
  masses=sort(unique(append(pks_Norharmane$mass,pks_Au$mass)))
  masses_diff=abs(diff(masses))
  masses_to_merge=which(masses_diff<mass_threshold)
  merged_masses=masses
  merged_masses[masses_to_merge]=(masses[masses_to_merge]+masses[masses_to_merge+1])/2
  merged_masses[masses_to_merge+1]=(masses[masses_to_merge]+masses[masses_to_merge+1])/2
  merged_masses=unique(merged_masses)
  #The solution should be recursive
  #create spectrum

  #SECTION 3:: Calculate correlation of each peak.

  #SECTION 4:: Remove peaks with negative correlation (assumed exogenous)

  #Print image
  rMSIproc::plotPeakImage(pks_Norharmane,c=2)
}


#CODE TO INCORPORATE
#clus <- kmeans(pks_Norharmane$intensity, centers = 2)
#rMSIproc::plotClusterImage(pks_Au, clus$cluster)

