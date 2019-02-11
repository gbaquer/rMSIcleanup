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

#
#   Load function:             'Ctrl + Shift + L'
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Run Roxygen2:              'Ctrl + Shift + D'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
#
#   Comment/Uncomment:         'Ctrl + Shift + C'

############################################################################



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
#' Plot PCA
#'
#' Plot PCA
#'
#' @param pks Peak Matrix
#' @param normalize Boolean value to determine if normalization is to be performed. If TRUE TIC normalization is performed
#' @param num_clusters Number of clusters
#' @param correlation_order Number of times that correlation needs to be applied.
#'
#' @return Nothing
#'
#'
#' @export
plot_pca <- function (pks,normalize=TRUE,num_clusters=2,correlation_order=0) {
  #SECTION 0 :: Preprocessing
  # Select first image if there are multiple
  if(length(pks$numPixels)>1)
  {
    rows=1:pks$numPixels[1]

    for(attr in attributes(pks)$names)
    {
      if(is.null(dim(pks[[attr]])))
      {
        if(attr!="mass")
        {
          pks[[attr]]=pks[[attr]][1]
        }
      }
      else
      {
        pks[[attr]]=pks[[attr]][rows,]
      }
    }
  }

  # Normalize to TIC
  if(normalize)
    pks$intensity=pks$intensity/pks$normalizations$TIC



  #SECTION 1 :: Identify regions outside of the tissue area

  clus=kmeans(pks$intensity,centers = 2)

  sorted_clusters=sort(clus$size,index.return=TRUE,decreasing = TRUE)$ix
  #[This section belongs in a plotting function]
  #Load ground truth
  annotations_Ag=read.table("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/comparativa_Ag_Au/Ag.ref")
  gt=generate_gt("Ag1",pks)#annotations_Ag[[2]]
  m=pks$mass
  pch_1=4
  pch_2=20
  gt_pch=abs(pch_1-pch_2)*(is.element(m,gt)==(pch_1>pch_2))+min(pch_1,pch_2)

  #Calculate spectral correlation
  corMAT_in=pks$intensity[which(clus$cluster==sorted_clusters[1]),]
  corMAT_out=pks$intensity[which(clus$cluster==sorted_clusters[2]),]
  corMAT=pks$intensity
  if(correlation_order==0)
  {
    corMAT_in=t(corMAT_in)
    corMAT_out=t(corMAT_out)
    corMAT=t(corMAT)
  }
  else if(correlation_order==-1)
  {
    corMAT=sum()
  }
  else
  {
    i=0
    while(i < correlation_order)
    {
      corMAT_in=cor(corMAT_in)
      corMAT_out=cor(corMAT_out)
      corMAT=cor(corMAT)
      i=i+1
    }
  }
  spec_clus_in=kmeans(corMAT_in,centers = 2)
  spec_clus_out=kmeans(corMAT_out,centers = 2)
  spec_clus=kmeans(corMAT,centers = 2)

  #PLOT PCA
  pca_in=prcomp(corMAT_in,center=TRUE,scale=TRUE)
  pca_plot_in=plot(pca_in$x[,1:2],pch=gt_pch,col=spec_clus_in$cluster,xlab=paste("PC1",round(100*pca_in$sdev[1]/sum(pca_in$sdev),2),"%"),ylab=paste("PC2",round(100*pca_in$sdev[2]/sum(pca_in$sdev)),"%"))
  title("In region")
  print(pca_plot_in)

  pca_out=prcomp(corMAT_out,center=TRUE,scale=TRUE)
  pca_plot_out=plot(pca_out$x[,1:2],pch=gt_pch,col=spec_clus_out$cluster,xlab=paste("PC1",round(100*pca_out$sdev[1]/sum(pca_out$sdev),2),"%"),ylab=paste("PC2",round(100*pca_out$sdev[2]/sum(pca_out$sdev)),"%"))
  title("Out region")
  print(pca_plot_out)

  pca=prcomp(corMAT,center=TRUE,scale=TRUE)
  pca_plot=plot(pca$x[,1:2],pch=gt_pch,col=spec_clus$cluster,xlab=paste("PC1",round(100*pca$sdev[1]/sum(pca$sdev),2),"%"),ylab=paste("PC2",round(100*pca$sdev[2]/sum(pca$sdev)),"%"))
  title("In+Out regions")
  print(pca_plot)

  #PLOT TSNE [Not giving valuable results]
  # tsne_in=tsne::tsne(scale(corMAT_in,center=TRUE,scale=TRUE))
  # tsne_plot_in=plot(tsne_in,pch=gt_pch,col=spec_clus_in$cluster)
  # title("In region")
  # print(tsne_plot_in)
  #
  # tsne_out=tsne::tsne(scale(corMAT_out,center=TRUE,scale=TRUE))
  # tsne_plot_out=plot(tsne_out,pch=gt_pch,col=spec_clus_out$cluster)
  # title("Out region")
  # print(tsne_plot_out)
  #
  # tsne_all=tsne::tsne(scale(corMAT,center=TRUE,scale=TRUE))
  # tsne_plot=plot(tsne_all,pch=gt_pch,col=spec_clus$cluster)
  # title("In+Out region")
  # print(tsne_plot)
}





#' Plot scores
#'
#' Returns several scores for performance assessment of a given binary classification result.
#'
#' @param x X variable
#' @param scores_list List of scores for each x
#'
#' @export
plot_scores <- function (x,scores_list) {

}



# __ RANDOM CHUNCKS OF CODE __
#CODE TO INCORPORATE
#clus <- kmeans(pks_Norharmane$intensity, centers = 2)
#rMSIproc::plotClusterImage(pks_Au, clus$cluster)

# __ LOADING THE LIBRARIES DOESN'T SEEM TO BE NECESSARY
#LOAD LIBRARIES
# library("R.utils")
# library("rMSI")
# library("rMSIproc")
# library("lattice")
# library("gridExtra")
# library("RColorBrewer")
# library("GlobalOptions")

# __ CODE FROM METHOD 2 COMPARE AU __
# for (i in 1:regions$num)
# {
#   # regions$mean_Norharmane=rbind(regions$mean_Norharmane,apply(pks_Norharmane_new$intensity[regions$pixels_Norharmane[[i]],],2,mean))
#   # regions$mean_Au=rbind(regions$mean_Au,apply(pks_Au_new$intensity[regions$pixels_Au[[i]],],2,mean))
#
#   regions$mean_Norharmane[i,]=apply(pks_Norharmane_new$intensity[regions$pixels_Norharmane[[i]],],2,mean)
#   regions$mean_Au[i,]=apply(pks_Au_new$intensity[regions$pixels_Au[[i]],],2,mean)
#
#   # regions$mean_Norharmane[i,]=pks_Norharmane_new$intensity[sample(regions$pixels_Norharmane[[i]],3),]
#   # regions$mean_Au[i,]=pks_Au_new$intensity[sample(regions$pixels_Au[[i]],3),]
# }
# for (i in 1:regions$num)
# {
#   #regions$mean_Norharmane[[i]]=apply(pks_Norharmane_new$intensity[regions$pixels_Norharmane[[i]],],2,mean)
#   #regions$mean_Au[[i]]=apply(pks_Au_new$intensity[regions$pixels_Au[[i]],],2,mean)
#
#   regions$mean_Norharmane[[i]]=pks_Norharmane_new$intensity[sample(1:length(regions$pixels_Norharmane[[i]]),3),]
#   regions$mean_Au[[i]]=pks_Au_new$intensity[sample(1:length(regions$pixels_Au[[i]]),3),]
# }
#

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

# __ OLD IMPLEMENTATION OF VERBOSE __
# #' Print verbose
# #'
# #' Print \code{message} if the verbose \code{level} is higher than the specified \code{threshold}.
# #'
# #'  1: Silent
# #'  0: Application messages
# #' -1: Code Sections
# #' -2: Code subsections
# #' -3: Debug mode
# #'
# #' @param message String to print out
# #' @param level Verbose level
# #'
# #' @return None
# #'
# #' @examples
# #' printv("Hello World!",0)
# #'
# #' @export
# printv <- function(message,level) {
#   threshold=0
#   if(level>=threshold)
#   {
#     print(message)
#   }
# }
