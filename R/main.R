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
library("RColorBrewer")

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
#' Plot image summary
#'
#' Plot the eman spectra and a few representative images of the two MSI images used.
#'
#' @return None
#'
#'
#' @export
plot_TIC <- function () {
  pks_Norharmane <- rMSIproc::LoadPeakMatrix("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/comparativa_matrix_au/peak_matrix_norharmane/mergeddata-peaks.zip")
  pks_Au <- rMSIproc::LoadPeakMatrix("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/comparativa_matrix_au/peak_matrix_au/mergeddata-peaks.zip")

  # TIC_Norharmane=apply(pks_Norharmane$intensity,2,sum)
  # TIC_Au=apply(pks_Au$intensity,2,sum)

  # plot(pks_Norharmane$mass,mean_Norharmane,type="h",col=1)
  # lines(pks_Au$mass,mean_Au,type="h",col=2)


  pks_Norharmane_plot=pks_Norharmane
  pks_Au_plot=pks_Au

  pks_Norharmane_plot$intensity[,1]=(apply(pks_Norharmane$intensity,1,sum))
  pks_Au_plot$intensity[,1]=(apply(pks_Au$intensity,1,sum))

  dev.new()
  rMSIproc::plotPeakImage(pks_Norharmane_plot,c=1)
  title("Norharmane")

  dev.new()
  rMSIproc::plotPeakImage(pks_Au_plot,c=1)
  title("Au")
}
#' Plot image summary
#'
#' Plot the eman spectra and a few representative images of the two MSI images used.
#'
#' @return None
#'
#'
#' @export
plot_image_summary <- function () {
  pks_Norharmane <- rMSIproc::LoadPeakMatrix("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/comparativa_matrix_au/peak_matrix_norharmane/mergeddata-peaks.zip")
  pks_Au <- rMSIproc::LoadPeakMatrix("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/comparativa_matrix_au/peak_matrix_au/mergeddata-peaks.zip")

  mean_Norharmane=apply(pks_Norharmane$intensity,2,mean,lwd = 3,lend=1)
  mean_Au=apply(pks_Au$intensity,2,mean,lwd = 3,lend=1)

  plot(pks_Norharmane$mass,mean_Norharmane,type="h",col=1)
  lines(pks_Au$mass,mean_Au,type="h",col=2)


  pks_Norharmane_plot=pks_Norharmane
  pks_Au_plot=pks_Au

  pks_Norharmane_plot$intensity[,1]=(apply(pks_Norharmane$intensity,1,mean))
  pks_Au_plot$intensity[,1]=(apply(pks_Au$intensity,1,mean))

  dev.new()
  rMSIproc::plotPeakImage(pks_Norharmane_plot,c=1)
  title("Norharmane")

  dev.new()
  rMSIproc::plotPeakImage(pks_Au_plot,c=1)
  title("Au")
}

#' Remove Matrix Padded
#'
#' Remove the matrix following the method presented by Fonville et al. 2012.
#' It takes a padded image where there is a wide enough region of pixels outside of the tissue.
#' It consists of three basic steps.
#'
#' @param pks_Norharmane Peak Matrix
#'
#' @return List of mass indices corresponding to the matrix peaks
#'
#'
#' @export
removeMatrix_padded <- function (pks_Norharmane) {
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
  # grid.arrange(levelplot(covMAT))
  # grid.arrange(levelplot(corMAT))

  #Rank correlation of exo indices [IMPROVE: get the "max_exo" elements that correlate the best to the three of them]
  max_exo=5
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

  centers=10
  palette(brewer.pal(n=centers,name = "Paired"))

  #Plot average spectra
  dev.new()
  mean_spectra=apply(pks_Norharmane$intensity,2,mean)
  colors=rep(1,length(mean_spectra))
  colors[nonbio_peaks]=2
  colors[bio_peaks[nonanatomical_peaks]]=3

  plot(pks_Norharmane$mass,mean_spectra,type="h",col=colors,lwd = 3,lend=1)
  legend(1000, 30, legend=c("Endogenous", "Matrix","Non-anatomical"),
         col=colors)


  #Print image
  pks_Norharmane_plot=pks_Norharmane
  pks_Norharmane_plot$intensity[,1]=(apply(pks_Norharmane$intensity[,bio_peaks[anatomical_peaks]],1,mean))
  pks_Norharmane_plot$intensity[,2]=(apply(pks_Norharmane$intensity[,nonbio_peaks],1,mean))
  pks_Norharmane_plot$intensity[,3]=(apply(pks_Norharmane$intensity[,bio_peaks[nonanatomical_peaks]],1,mean))
  pks_Norharmane_plot$intensity[,4]=(apply(pks_Norharmane$intensity[,union(nonbio_peaks,bio_peaks[nonanatomical_peaks])],1,mean))
  dev.new()
  rMSIproc::plotPeakImage(pks_Norharmane_plot,c=1)
  title("Endogenous peaks")
  dev.new()
  rMSIproc::plotPeakImage(pks_Norharmane_plot,c=2)
  title("Matrix peaks")
  dev.new()
  rMSIproc::plotPeakImage(pks_Norharmane_plot,c=3)
  title("Non-anaotomical peaks")
  dev.new()
  rMSIproc::plotPeakImage(pks_Norharmane_plot,c=4)
  title("Matrix + Non-anaotomical peaks")

  #Return the peaks corresponding to the matrix
  return(union(nonbio_peaks,bio_peaks[nonanatomical_peaks]))
}

#' Remove Matrix Compare Au
#'
#' Remove the matrix by comparing anatomically similar regions of the tissue to a reference image taken
#' with a gold matrix. Assesses the degree of similarity between peaks
#' It takes a padded image where there is a wide enough region of pixels outside of the tissue.
#' It consists of three basic steps.
#'
#' @param pks_Au Peak Matrix with gold coating
#' @inherit removeMatrix_padded
#'
#'
#'
#' @export
removeMatrix_compareAu <- function (pks_Norharmane,pks_Au) {
  #SECTION 1:: Get regions of similarity
  regions <- list(
    num = 0,
    centers_Norharmane = t(matrix(c(130,54,61,60,111,92,14,71,119,30,67,46), nrow=2)),
    centers_Au =t(matrix(c(103,54,54,79,78,85,15,45,84,33,45,36), nrow=2)),
    radial_point_Norharmane =t(matrix(c(137,47,56,85,112,87,5,71,114,27,68,41), nrow=2)),
    radial_point_Au=t(matrix(c(106,63,56,83,76,79,8,41,80,33,43,30), nrow=2)),
    radius_Norharmane =matrix(),
    radius_Au =matrix(),
    pixels_Norharmane =list(),
    pixels_Au=list(),
    mean_Norharmane=array(),
    mean_Au=array()
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

  # SECTION 2: Set the two images to have the same masses, Calibration, alignment
  #In order to compare correlations they need to have the same size. Standardize and align spectra.
  pks_Au_new=pks_Au
  pks_Norharmane_new=pks_Norharmane

  mass_threshold=2
  masses= sort(union(pks_Norharmane$mass,pks_Au$mass))


  indices_Au=rep(1,length(pks_Au$mass))
  for(i in 1:length(pks_Au$mass))
  {
    indices_Au[i]=which(masses==pks_Au$mass[i])
  }

  indices_Norharmane=rep(1,length(pks_Norharmane$mass))
  for(i in 1:length(pks_Norharmane$mass))
  {
    indices_Norharmane[i]=which(masses==pks_Norharmane$mass[i])
  }

  indices_Norharmane_backwards=rep(1,length(masses))
  for(i in 1:length(masses))
  {
    index=which(pks_Norharmane$mass==masses[i])
    if(length(index)>0)
    {
      indices_Norharmane_backwards[i]=index
    }
    else
    {
      indices_Norharmane_backwards[i]=NA
    }
  }

  pks_Au_new$mass=masses
  pks_Au_new$intensity=matrix(0,pks_Au_new$numPixels,length(masses))
  pks_Au_new$area=matrix(0,pks_Au_new$numPixels,length(masses))
  pks_Au_new$SNR=matrix(0,pks_Au_new$numPixels,length(masses))
  pks_Au_new$intensity[,indices_Au]=pks_Au$intensity
  pks_Au_new$area[,indices_Au]=pks_Au$area
  pks_Au_new$SNR[,indices_Au]=pks_Au$SNR

  pks_Norharmane_new$mass=masses
  pks_Norharmane_new$intensity=matrix(0,pks_Norharmane_new$numPixels,length(masses))
  pks_Norharmane_new$area=matrix(0,pks_Norharmane_new$numPixels,length(masses))
  pks_Norharmane_new$SNR=matrix(0,pks_Norharmane_new$numPixels,length(masses))
  pks_Norharmane_new$intensity[,indices_Norharmane]=pks_Norharmane$intensity
  pks_Norharmane_new$area[,indices_Norharmane]=pks_Norharmane$area
  pks_Norharmane_new$SNR[,indices_Norharmane]=pks_Norharmane$SNR

  #Missing calibration and alignment
  masses_to_merge=which(diff(masses)<mass_threshold)

  pks_Au_new$intensity[,masses_to_merge]=pks_Au_new$intensity[,masses_to_merge]+pks_Au_new$intensity[,masses_to_merge+1]
  pks_Au_new$area[,masses_to_merge]=pks_Au_new$area[,masses_to_merge]+pks_Au_new$area[,masses_to_merge+1]
  pks_Au_new$SNR[,masses_to_merge]=pks_Au_new$SNR[,masses_to_merge]+pks_Au_new$SNR[,masses_to_merge+1]

  pks_Norharmane_new$intensity[,masses_to_merge]=pks_Norharmane_new$intensity[,masses_to_merge]+pks_Norharmane_new$intensity[,masses_to_merge+1]
  pks_Norharmane_new$area[,masses_to_merge]=pks_Norharmane_new$area[,masses_to_merge]+pks_Norharmane_new$area[,masses_to_merge+1]
  pks_Norharmane_new$SNR[,masses_to_merge]=pks_Norharmane_new$SNR[,masses_to_merge]+pks_Norharmane_new$SNR[,masses_to_merge+1]

  #SECTION 3:: Calculate, median and average spectrum
  regions$mean_Norharmane=array(0,c(regions$num,length(masses)))
  regions$mean_Au=array(0,c(regions$num,length(masses)))

  for (i in 1:regions$num)
  {
    # regions$mean_Norharmane=rbind(regions$mean_Norharmane,apply(pks_Norharmane_new$intensity[regions$pixels_Norharmane[[i]],],2,mean))
    # regions$mean_Au=rbind(regions$mean_Au,apply(pks_Au_new$intensity[regions$pixels_Au[[i]],],2,mean))

    regions$mean_Norharmane[i,]=apply(pks_Norharmane_new$intensity[regions$pixels_Norharmane[[i]],],2,mean)
    regions$mean_Au[i,]=apply(pks_Au_new$intensity[regions$pixels_Au[[i]],],2,mean)

    # regions$mean_Norharmane[i,]=pks_Norharmane_new$intensity[sample(regions$pixels_Norharmane[[i]],3),]
    # regions$mean_Au[i,]=pks_Au_new$intensity[sample(regions$pixels_Au[[i]],3),]
  }
  #SECTION 4:: Calculate correlation of each peak.
  # for (i in 1:regions$num)
  # {
  #   #regions$mean_Norharmane[[i]]=apply(pks_Norharmane_new$intensity[regions$pixels_Norharmane[[i]],],2,mean)
  #   #regions$mean_Au[[i]]=apply(pks_Au_new$intensity[regions$pixels_Au[[i]],],2,mean)
  #
  #   regions$mean_Norharmane[[i]]=pks_Norharmane_new$intensity[sample(1:length(regions$pixels_Norharmane[[i]]),3),]
  #   regions$mean_Au[[i]]=pks_Au_new$intensity[sample(1:length(regions$pixels_Au[[i]]),3),]
  # }
  #
  #SECTION 5:: Remove peaks with negative correlation (assumed exogenous)
  # for (i in 1:regions$num)
  # {
  #   dev.new()
  #   levelplot(cor(rbind(regions$mean_Norharmane[[i]],regions$mean_Au[[i]])))
  # }
  correlation_vector=array(0,c(1,length(masses)))
  for(i in 1:length(masses))
  {
    correlation_vector[i]=cor(regions$mean_Norharmane[,i],regions$mean_Au[,i])
  }
  correlation_threshold=0.8
  return(indices_Norharmane_backwards[which(correlation_vector>correlation_threshold)])
  #Print image
  #rMSIproc::plotPeakImage(pks_Norharmane,c=2)

  #Return the peaks corresponding to the matrix according to region 1
  # correlation_results<-cor(rbind(regions$mean_Norharmane[[1]],regions$mean_Au[[1]]))
  # d<-density(correlation_results)
  # dev.new()
  # levelplot(correlation_results)
  # dev.new()
  # plot(d)
  #
  # correlation_threshold=0.9
  # return(union(nonbio_peaks,bio_peaks[nonanatomical_peaks]))
}

#' Remove Matrix k-Means Transpose
#'
#' Remove the matrix by performing k-means clustering on the transpose of the peak matrix.
#' The algorithm identifies similar spectral peaks.
#'
#' @inherit removeMatrix_padded
#'
#'
#'
#' @export
removeMatrix_kMeansTranspose <- function (pks_Norharmane) {
  #pks_Norharmane=pks_Au
  #scale abans del k-means
  centers=10
  data=scale((t((pks_Norharmane$intensity))))
  clus <- kmeans(data, centers = centers)
  pks_cluster_mean=pks_Norharmane
  for(attr in attributes(pks_cluster_mean)$names)
  {
     pks_cluster_mean[[attr]]=double()
  }
  for(i in 1:centers )
  {
    #Compute average
    spectra=matrix(0, dim(pks_Norharmane$intensity)[1],dim(pks_Norharmane$intensity)[2])
    cluster_spectra=pks_Norharmane$intensity[,which(clus$cluster==i)]
    if(!is.null(dim(cluster_spectra)))
      spectra[,1]=apply(cluster_spectra,1,mean)
    else
      spectra[,1]=cluster_spectra
    pks_cluster_mean$intensity=rbind(pks_cluster_mean$intensity,log(spectra))
    for(attr in attributes(pks_cluster_mean)$names)
    {
      if(attr!="intensity")
        pks_cluster_mean[[attr]]=rbind(pks_cluster_mean[[attr]],pks_Norharmane[[attr]])
    }
    #Show image
    #grid.arrange(grob(rMSIproc::plotPeakImage(pks_cluster_mean,c=i)))
  }
  dev.new()
  rMSIproc::plotPeakImage(pks_cluster_mean,c=1)

  #Adjust color palette
  palette(brewer.pal(n=centers,name = "Paired"))

  #Perform & plot PCA
  dev.new()
  pca=prcomp(data)
  #plot(pca$rotation) # loadings
  plot(pca$x[,1:2],col=clus$cluster,xlab=paste("PC1",round(100*pca$sdev[1]/sum(pca$sdev),2),"%"),ylab=paste("PC2",round(100*pca$sdev[2]/sum(pca$sdev)),"%"),lwd=3)

  #Intra-cluster variance
  dev.new()
  intracluster_sd=rep(0,centers)
  for(i in 1:centers)
    intracluster_sd[i]=sd(pks_Norharmane$intensity[,which(clus$cluster==i)])
  plot(1:centers,log(intracluster_sd),type="h",col=1:centers,lwd = 10,lend=1)

  #Cluster assignement
  dev.new()
  mean_spectra=apply(pks_Norharmane$intensity,2,mean)
  plot(pks_Norharmane$mass,mean_spectra,type="h",col=clus$cluster,lwd = 3,lend=1)

  #Volcano plot
  #Return the peaks corresponding to the matrix
  # correlation_threshold=0.7
  # return(union(nonbio_peaks,bio_peaks[nonanatomical_peaks]))
}

#' Remove Matrix k-Means Transpose based on correlation
#'
#' Remove the matrix by performing k-means clustering on the transpose of the correlation of the peak matrix.
#' The algorithm identifies similar spectral peaks.
#'
#' @inherit removeMatrix_padded
#'
#'
#'
#' @export
removeMatrix_kMeansTransposeCor <- function (pks_Norharmane) {
  #Proposed procedure
  #dev.new()
  centers=10
  norm_data=pks_Norharmane$intensity/pks_Norharmane$normalizations$TIC
  cor_data=cor(norm_data)
  scaled_data=scale(cor_data)
  dev.new()
  levelplot(scaled_data)
  #clustering
  clus <- kmeans(scaled_data, centers = centers)

  #Calculate mean spectra of each cluster for plotting
  pks_cluster_mean=pks_Norharmane
  for(attr in attributes(pks_cluster_mean)$names)
  {
    pks_cluster_mean[[attr]]=double()
  }
  for(i in 1:centers )
  {
    #Compute average
    spectra=matrix(0, dim(pks_Norharmane$intensity)[1],dim(pks_Norharmane$intensity)[2])
    cluster_spectra=pks_Norharmane$intensity[,which(clus$cluster==i)]
    if(!is.null(dim(cluster_spectra)))
      spectra[,1]=apply(cluster_spectra,1,mean)
    else
      spectra[,1]=cluster_spectra
    pks_cluster_mean$intensity=rbind(pks_cluster_mean$intensity,log(spectra))
    for(attr in attributes(pks_cluster_mean)$names)
    {
      if(attr!="intensity")
        pks_cluster_mean[[attr]]=rbind(pks_cluster_mean[[attr]],pks_Norharmane[[attr]])
    }
    #Show image
    #grid.arrange(grob(rMSIproc::plotPeakImage(pks_cluster_mean,c=i)))
  }
  dev.new()
  rMSIproc::plotPeakImage(pks_cluster_mean,c=1)

  #Adjust color palette
  palette(brewer.pal(n=centers,name = "Paired"))

  #Perform & plot PCA
  dev.new()
  pca=prcomp(scaled_data)
  #plot(pca$rotation) # loadings
  plot(pca$x[,1:2],col=clus$cluster,xlab=paste("PC1",round(100*pca$sdev[1]/sum(pca$sdev),2),"%"),ylab=paste("PC2",round(100*pca$sdev[2]/sum(pca$sdev)),"%"),lwd=3)

  #Intra-cluster variance
  dev.new()
  intracluster_sd=rep(0,centers)
  for(i in 1:centers)
    intracluster_sd[i]=sd(pks_Norharmane$intensity[,which(clus$cluster==i)])
  plot(1:centers,log(intracluster_sd),type="h",col=1:centers,lwd = 10,lend=1)

  #Cluster assignement
  dev.new()
  mean_spectra=apply(pks_Norharmane$intensity,2,mean)
  plot(pks_Norharmane$mass,mean_spectra,type="h",col=clus$cluster,lwd = 3,lend=1)

  #Volcano plot

  #Return the peaks corresponding to the matrix
  # correlation_threshold=0.7
  # return(union(nonbio_peaks,bio_peaks[nonanatomical_peaks]))
}

#' Cross validation.
#'
#' Compare the results of all methods to assess consistency.
#'
#' @return None
#'
#'
#' @export
cross_validation <- function () {
  #LOAD DATA
  pks_Norharmane <- rMSIproc::LoadPeakMatrix("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/comparativa_matrix_au/peak_matrix_norharmane/mergeddata-peaks.zip")
  pks_Au <- rMSIproc::LoadPeakMatrix("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/comparativa_matrix_au/peak_matrix_au/mergeddata-peaks.zip")

  #METHOD 1
  removeMatrix_padded(pks_Norharmane)
  #METHOD 2
  removeMatrix_compareAu(pks_Norharmane,pks_Au)
  #METHOD 3
  removeMatrix_kMeansTranspose(pks_Norharmane)
  #METHOD 3.1
  removeMatrix_kMeansTransposeCor(pks_Norharmane)
  #COMPARE RESULTS

}

#' Return 2
#'
#' Dummy function
#'
#' @return None
#'
#' @export
return_2 <- function() {
  return(2)
}

# __ RANDOM CHUNCKS OF CODE __

#CODE TO INCORPORATE
#clus <- kmeans(pks_Norharmane$intensity, centers = 2)
#rMSIproc::plotClusterImage(pks_Au, clus$cluster)

# __ TRYING TO PLOT MULTIPLE IMAGES __

# removeMatrix_kMeansTranspose <- function () {
#   pks_Norharmane <- rMSIproc::LoadPeakMatrix("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/comparativa_matrix_au/peak_matrix_norharmane/mergeddata-peaks.zip")
#   #scale abans del k-means
#   centers=10
#   clus <- kmeans(t(pks_Norharmane$intensity), centers = centers)
#   pks_cluster_mean=pks_Norharmane
#   for(i in 1:centers )
#   {
#     #Compute average
#     spectra=matrix(0, dim(pks_Norharmane$intensity)[1],dim(pks_Norharmane$intensity)[2])
#     cluster_spectra=pks_Norharmane$intensity[,which(clus$cluster==i)]
#     if(!is.null(dim(cluster_spectra)))
#       spectra[,1]=apply(cluster_spectra,1,mean)
#     else
#       spectra[,1]=cluster_spectra
#     pks_cluster_mean$intensity=append(pks_cluster_mean$intensity,spectra)
#     for(attr in attributes(pks_cluster_mean)$names)
#     {
#       if(attr!="intensity")
#         pks_cluster_mean[[attr]]=cbind(pks_cluster_mean[[attr]],pks_Norharmane[[attr]])
#     }
#     #Show image
#     #grid.arrange(grob(rMSIproc::plotPeakImage(pks_cluster_mean,c=i)))
#   }
#   rMSIproc::plotPeakImage(pks_cluster_mean,c=1)
# }

# __ CALIBRATION AND ALIGNMENT __
# masses_diff=abs(diff(masses))
# masses_to_merge=which(masses_diff<mass_threshold)
# merged_masses=masses
# merged_masses[masses_to_merge]=(masses[masses_to_merge]+masses[masses_to_merge+1])/2
# merged_masses[masses_to_merge+1]=(masses[masses_to_merge]+masses[masses_to_merge+1])/2
# merged_masses=unique(merged_masses)
#The solution should be recursive
