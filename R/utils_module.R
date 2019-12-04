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
  ),
  round_digits=list(.value=4,
                     .read.only=FALSE,
                     .validate= function(x) is.numeric(x) && x%%1==0 && x>=0,
                     .failed_msg = "Round digits should be a positive integer",
                     .description = "Number of digits to be used in all the printing functions"
  )
)

#' Plot text
#'
#' Plots text with white background
#'
#' @param text Text to be plotted
#' @param size Size of the text
#'
plot_text <- function(text,size=3)
{
  p=ggplot() + geom_label()+ annotate("text", x = 0:1, y = c(1,0), size=size, label = c(text,NA), vjust=1,hjust=0) +
    #coord_fixed(sqrt(2),expand = F)+
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())
  return(p)
}

#' Plot peak image
#'
#' Plots a peak image
#'
#' @param pks Peak matrix
#' @param values Values to be plotted
#' @param title Title of the plot
#' @param isNA Is NA
#' @param chosen Is chosen
#' @param in_pksMat Is in pksMat
#'
ggplot_peak_image <- function(pks,values=NA,title="",isNA=F,chosen=T,in_pksMat=T)
{
  # zplots<-matrix(NA, nrow=max(pks$pos[,1]), ncol=max(pks$pos[,2]))
  # for( i in 1:nrow(pks$pos))
  # {
  #   zplots[pks$pos[ i , 1 ], pks$pos[ i , 2 ]] <- values[i]
  # }
  #levelplot(zplots)
  x=NULL
  y=NULL
  z=NULL
  coul <- viridis(100)
  if(isNA)
  {
    p=ggplot() + geom_label()+ annotate("text", x = 0, y = 0, size=5, label = "NA") +
      ggtitle(title) +
      theme(axis.line=element_blank(),axis.text.x=element_blank(),
            axis.text.y=element_blank(),axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),legend.position="none",panel.border=element_blank(),panel.grid.major=element_blank(), plot.title = element_text(size=8),
            panel.grid.minor=element_blank(),plot.background=element_blank())
  }
  else
  {
    if(chosen)
      background=element_rect(fill = "#45F442FF")
    else
      background=element_rect(fill = "#F46242FF")
    # if(in_pksMat)
    #   if(chosen)
    #     background=element_rect(fill = "#45F442FF")
    #   else
    #     background=element_rect(fill = "#F46242FF")
    # else
    #   background=element_rect(fill = "#918A88FF")

    df=data.frame(x=pks$pos[,2],y=pks$pos[,1],z=values)
    p=ggplot(df, aes(x, y, fill = z)) + geom_raster() +
      coord_fixed(1,expand = F) +
      ggtitle(title) +
      # scale_fill_gradientn(colours=append("#000000FF",rev(rainbow(200,start=5.5/6,end=4.2/6)))) +
      scale_fill_gradientn(colours=coul)+
      theme(axis.line=element_blank(),axis.text.x=element_blank(),
            axis.text.y=element_blank(),axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),legend.position="none",
            panel.background=element_rect(fill = "#000000FF"),panel.border=element_blank(),panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),plot.background=background,plot.title = element_text(hjust = 0.5,size=9))
  }
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

#' Find First GT
#'
#' Returns the index of the first element in x greater than k
#'
#' @param x Input vector
#' @param k Threshold
#'
#' @return x with k precision
find_first_gt <- function(x,k) if(k>=max(x)) length(x) else match(T,x>k)

#' Find First LT
#'
#' Returns the index of the first element in x lower than k
#'
#' @param x Input vector
#' @param k Threshold
#'
#' @return x with k precision
find_first_lt <- function(x,k) match(T,x<k)



#' Generate report spectral clustering
#'
#' Generates a pdf report with the results of a spectral clustering
#'
#' @param pks Input Peak Matrix
#' @param mass_range Vector containing the minimum and maximum masses to be used in the process. The default is the complete mass range.
#' @param num_clus Number of clusters to use. Defaults to 3.
#' @param f Function to use when plotting the resulting spectral clusters in the spatial domain. The function should take a numerical vector and return a single number. Default is mean
#' @param generate_pdf Boolean value indicating if a pdf report needs to be produced
#' @param normalization String indicating the normalization technique to use. Possible values: "None", "TIC","RMS", "MAX" or "AcqTic"
#' @param legend_pos String indicating the legend position. "bottomright","bottomleft", "topleft" or "topright"
#'
#' @return Matrix containing the masses processed along with the cluster to which each of them is assigned
#'
#' @export

spectral_clustering <- function(pks, mass_range=c(min(pks$mass),max(pks$mass)), num_clus=3, f=mean,generate_pdf=F, normalization="None",legend_pos="topright")
{
  # Normalize
  if(normalization!="None")
    pks$intensity <- pks$intensity/pks$normalizations[[normalization]]
  # Prune masses
  cols=find_first_gt(pks$mass,mass_range[1]):find_first_gt(pks$mass,mass_range[2])
  pks$intensity<-pks$intensity[,cols]
  pks$mass<-pks$mass[cols]

  # PCA of spectral correlation
  pca <- prcomp(cor(pks$intensity))
  # Clustering
  clus <- kmeans(pca$x, centers = num_clus,iter.max = 100)
  labels = paste("Cluster",1:length(clus$size))

  # Open pdf report
  if(generate_pdf)
  {
    #Generate pdf name
    pdf_file=generate_file_name(pks$names[1])
    #open pdf file
    a4_height=8.27
    a4_width=11.69
    pdf(pdf_file,paper='a4r',width=a4_width,height=a4_height)

    #First page of metadata [File name, mean image, matrix formula, adducts list, base forms, S1_threshold]
    text=""
    text=add_entry(text,"#")
    text=add_entry(text,"- Package Version:",packageVersion("rMSIcleanup"))
    text=add_entry(text,"- Time:",as.character(Sys.time()))
    text=add_entry(text,"#")
    text=add_entry(text,"IMAGE INFORMATION")
    text=add_entry(text,"- Peak matrix:",pks$names[1])
    text=add_entry(text,"- Number of peaks:",length(pks$mass))
    text=add_entry(text,"- Number of pixels:",pks$numPixels[1])
    text=add_entry(text,"- Mass Range: [",min(pks$mass),", ", max(pks$mass),"]")
    text=add_entry(text,"#")
    text=add_entry(text,"CLUSTERING INFORMATION")
    text=add_entry(text,"- Number of clusters:",num_clus)
    text=add_entry(text,"#")
    print(plot_text(text))
  }
  # Plot PC1 vs. PC2
  pca_plot <- plot(pca$x[,c(1,2)],
                   col=clus$cluster,
                   xlab=paste("PC1",round(100*pca$sdev[1]/sum(pca$sdev),2),"%"),
                   ylab=paste("PC2",round(100*pca$sdev[2]/sum(pca$sdev)),"%"))
  legend("bottomright", legend=labels,
         col=1:length(clus$size), lty=c(1,1), cex=0.8, lwd=2)


  # Plot spectral clustering results
  mean_spectra=apply(pks$intensity,2,mean)
  plot(pks$mass,mean_spectra,type="h",col=clus$cluster,lwd = 3,lend=1)
  legend(legend_pos, legend=labels,
         col=1:length(clus$size), lty=c(1,1), cex=0.8, lwd=2)

  # Plot cluster images
  # Transform the spectral images of each cluster into a single spectral image to be plotted
  page_layout=rbind(c(1,3),
                    c(1,3),
                    c(2,4))
  plts=list()
  for(i in 1:length(clus$withins))
  {
    image=apply(pks$intensity[,which(clus$cluster==i)],1,f)
    #rMSIproc::plotValuesImage( peakMatrix = pks, values = image,labels = labels[i])
    plts=append(plts,list(ggplot_peak_image(pks,image,labels[i])))
    text=add_entry("",paste(round(pks$mass[which(clus$cluster==i)],pkg_opt("round_digits")),collapse=" ; "))
    plts=append(plts,list(plot_text(text)))
    if(i%%2==0||i==length(clus$withins))
    {
      grid.arrange(grobs=plts,layout_matrix=page_layout)
      plts=list()
    }
  }

  # Close pdf report
  if(generate_pdf)
    dev.off()

  return(list(mass=pks$mass,cluster=clus$cluster))
}



#PDF FILE PRINTING UTILS

#' Generate File Name
#'
#' Generates a file name that has not yet been used
#'
#' @param base_name Base name including
#' @param extension File extension. Defaults to ".pdf"
#' @param folder Name of the folder in which to store the file. Defaults to "output/"
#'
#' @return File name that is not used in the
#'
#' @export

generate_file_name <- function(base_name, extension=".pdf",folder="output/")
{
  file_name=paste(folder,base_name,"_000",extension,sep="")
  i=0
  while(file.exists(file_name))
  {
    i=i+1
    file_name=paste(folder,base_name,"_",str_pad(i, 3, pad = "0"),extension,sep="")
  }
  return(file_name)
}

#' Add entry
#'
#' Adds an entry to a string of text to be witten to pdf
#'
#' @param text Starting text string
#' @param ... Contents of the new entry. They will be pasted together with a space. Any vectors will be separated with a semicolon. Any # will be replaced by a division line
#'
#' @return Updated text variable with the new entry
#'
#' @export

add_entry <- function(text,...)
{
  arguments <- list(...)
  # Change division line
  are_division=which(arguments=="#")
  arguments[are_division]=paste(rep("#",80),collapse="")
  # Convert vectors to strings
  are_vectors=which(lapply(arguments,length)>1)
  arguments[are_vectors]=lapply(arguments[are_vectors],function(x)paste(x,collapse="; "))
  #Append new entry
  text=paste(text,paste(strwrap(paste(arguments,collapse=" ")),collapse = "\n"),"\n",collapse="")
  return(text)
}

#' Get one peak matrix
#'
#' Returns the peak matrix in pks corresponding to the ith image
#'
#' @param pks Peak Matrix
#' @param i Index of the peak matrix to be retieved
#'
#' @return Peak matrix containing only the first image from pks
#'
#' @export
get_one_peakMatrix <- function(pks,i=1)
{
  if(length(pks$numPixels)>1)
  {
    rows=1:pks$numPixels[i]
    if(i>1&i<=length(pks$numPixels))
      rows=rows+sum(pks$numPixels[1:i-1])

    for(attr in attributes(pks)$names)
    {
      if(is.null(dim(pks[[attr]])))
      {
        if(attr!="mass")
        {
          pks[[attr]]=pks[[attr]][i]
        }
      }
      else
      {
        pks[[attr]]=pks[[attr]][rows,]
      }
    }
  }
  return(pks)
}

#' Get one peak matrix
#'
#' Returns the peak matrix in pks corresponding to the ith image
#'
#' @param pks Peak Matrix
#' @param i Index of the peak matrix to be retieved
#'
#' @return Peak matrix containing only the first image from pks
#'
#' @export
get_columns_peakMatrix <- function(pks,cols=seq_along(pks$mass))
{
  pks$mass=pks$mass[cols]
  pks$intensity=pks$intensity[,cols]
  pks$area=pks$area[,cols]
  pks$SNR=pks$SNR[,cols]
  return(pks)
}

#' Get closest peak
#'
#' Returns the indices of "experimental_masses" that are the closest to each element in "calc_masses"
#'
#' @param calc_masses Vector with calculated masses
#' @param experimental_masses Vector with experimental masses
#'
#' @return Peak matrix containing only the first image from pks
#'
#' @export
get_closest_peak <- function(calc_masses,experimental_masses)
  apply(abs(outer(calc_masses,experimental_masses,'-')),1,function(x) sort(x,index.return=TRUE)$ix[1])


#' Quiet
#'
#' Forces a function to execute without printing any messages
#'
#' @param x Call to function
#'
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}
#' Quiet
#'
#' Forces a function to execute without printing any messages
#'
#' @param experimental_mass Experimental mass
#' @param calculated_mass Calculated mass
#' @param full_spectrum_masses Full spectrum mass vector
#' @param tol_scans Tolerance specified in scans
#'
#' @return Boolean value indicating whether the experimental mass is within "tol_scans" of the calculated_mass
is_within_scan_tol <- function(experimental_mass,calculated_mass,full_spectrum_masses,tol_scans)
{
  i=which.min(abs(full_spectrum_masses-calculated_mass))
  cols=(i-tol_scans):(i+tol_scans)
  mass_range=full_spectrum_masses[cols]
  return((experimental_mass>min(mass_range))&(experimental_mass<max(mass_range)))
}

#' Pareto Front
#'
#' Returns the pareto front for a set of three dimensional points
#'
#' @param x Vector of first coordinates
#' @param y Vector of second coordinates
#' @param z Vector of third coordinates
#'
#' @return Boolean value indicating whether the experimental mass is within "tol_scans" of the calculated_mass
get_pareto_front <- function(x,y,z)
{
  d = data.frame(x,y,z)
  D = d[order(d$x,d$y,d$z,decreasing=T),]
  front = D[union(which(!duplicated(cummax(D$y))),which(!duplicated(cummax(D$z)))),]

  return(as.matrix(front))
}

# SIMILARITY MEASURES

#' Exponential Decay Similarity
#'
#' Determines the degree of similarity between a vector a and a vector b as 1/exp(d) where d corresponds to the distance between the two normalized vectors calculated according to the specified method
#'
#' @param a Vector a
#' @param b Vector b
#' @param method The distance measure to be used. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". Any unambiguous substring can be given.
#' @param normalize Boolean value indicating whether the a and b vectors should be normalixzed
#'
#' @return Similarity between a and b
exponential_decay_similarity <- function(a,b,method="euclidian",normalize=T)
{
  #Normalize
  if(normalize&&!is.infinite(max(a))&&!is.infinite(max(b)))
  {
    a=a/max(a)
    b=b/max(b)
  }

  #Calculate distance
  d=dist(rbind(a,b),method=method)

  #Calculate similarity
  s=1/exp(d)


  return(s)
}

# SIMILARITY MEASURES

#' Exponential Decay Similarity
#'
#' Determines the degree of similarity between a vector a and a vector b as 1/exp(d) where d corresponds to the distance between the two normalized vectors calculated according to the specified method
#'
#' @param a Vector a
#' @param b Vector b
#' @param method The distance measure to be used. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski". Any unambiguous substring can be given.
#' @param normalize Boolean value indicating whether the a and b vectors should be normalixzed
#'
#' @return Similarity between a and b
centroid2peakmatrix <- function(file_path)
{
  #Import an imzML using centroid mode (after peak picking).
  pkLst <- rMSIproc::import_imzMLpeakList(file_path)

  #Run the peak-binning to get the peak matrix.
  myPkMat <- rMSIproc::PeakList2PeakMatrix(pkLst$peakList, pkLst$pos, BinTolerance = 5, BinToleranceUsingPPM = T)

  #Save the peak matrix.
  rMSIproc::StorePeakMatrix(paste(tools::file_path_sans_ext(file_path),".zip",sep=""), myPkMat)

  #Print validation
  name_Ag<-c("Ag1_100","Ag1_93","Ag2_54","Ag2_100","Ag2_47","Ag3_36","Ag3_100","Ag3_93","Ag3_29","Ag4_19","Ag4_72","Ag4_100","Ag4_62","Ag4_14","Ag5_12","Ag5_54","Ag5_100","Ag5_93","Ag5_43","Ag6_35","Ag6_81","Ag6_100","Ag6_70","Ag6_26","Ag7_23","Ag7_65","Ag7_100","Ag7_93","Ag7_52","Ag7_16","Ag8_14","Ag8_46","Ag8_86","Ag8_100","Ag8_74","Ag8_35","Ag9_33","Ag9_72","Ag9_100","Ag9_93","Ag9_58","Ag9_23","Ag10_22","Ag10_55","Ag10_90","Ag10_100","Ag10_77","Ag10_41","Ag10_14","Ag11_15","Ag11_41","Ag11_77","Ag11_100","Ag11_93","Ag11_62","Ag11_29","Ag12_30","Ag12_62","Ag12_92","Ag12_100","Ag12_80","Ag12_46","Ag12_19")
  mass_Ag<-c(106.9051,108.90475,213.81019,215.80985,217.80951,320.71529,322.71494,324.7146,326.71426,427.62038,429.62004,431.6197,433.61936,435.61902,534.52548,536.52513,538.52479,540.52445,542.52411,643.43023,645.42989,647.42955,649.42921,651.42887,750.33532,752.33498,754.33464,756.3343,758.33396,760.33362,857.24042,859.24008,861.23974,863.2394,865.23906,867.23871,966.14517,968.14483,970.14449,972.14415,974.14381,976.14347,1073.05027,1075.04993,1077.04959,1079.04925,1081.0489,1083.04856,1085.04822,1179.95536,1181.95502,1183.95468,1185.95434,1187.954,1189.95366,1191.95332,1288.86012,1290.85978,1292.85944,1294.85909,1296.85875,1298.85841,1300.85807)
  pks_mass<-myPkMat$mass
  rel_error=lapply(mass_Ag,function(m) min(abs((pks_mass-m)/pks_mass)))
  tol=5-6
  print(name_Ag[which(rel_error<tol)])
}

# if(is.na(s))
#   s=0

#RANDOM CHUNCKS OF CODE
# Cosine similarity two vectors
#
# Returns the cosine similarity between two vectors
#
# @param ix Input vectors
#
# cos_sim_vectors <- function(ix)
# {
#   A = X[ix[1],]
#   B = X[ix[2],]
#   return( sum(A*B)/sqrt(sum(A^2)*sum(B^2)) )
# }

# Cosine similarity two vectors
#
# Returns the cosine similarity between two vectors
#
# @param X Input matrix
#
# cos_sim <- function(X)
# {
#   n <- nrow(X)
#   cmb <- expand.grid(i=1:n, j=1:n)
#   C <- matrix(apply(cmb,1,cos_sim_vectors),n,n)
#   return(C)
# }

# plot_text <- function(text, size=5)
# {
#   default <- par()
#   par(mar = rep(1, 4))
#   par(oma = rep(0, 4))
#   plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',
#        xaxt='n', yaxt='n', xlab='', ylab='')
#   text(0,1,text, pos=4)
#   par(default)
# }


# plot_text <- function(text,size=5)
# {
#   p=ggplot() + geom_label()+ annotate("text", x = 0:1, y = c(1,0), size=size, label = c(text,NA), vjust=1,hjust=0) + scale_x_continuous(expand=c(0,10))+ scale_y_continuous(expand=c(0,sqrt(2)*10)) #+ coord_fixed(sqrt(2))#+
#     # scale_x_continuous(expand=c(0,0)) +
#     # scale_y_continuous(expand=c(0,0)) #+ theme(axis.line=element_blank(),axis.text.x=element_blank(),
#   #                                             axis.text.y=element_blank(),axis.ticks=element_blank(),
#   #                                             axis.title.x=element_blank(),
#   #                                             axis.title.y=element_blank(),legend.position="none",
#   #                                             panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
#   #                                             panel.grid.minor=element_blank(),plot.background=element_blank())
#   return(p)
# }

#FROM GENERATE GROUND TRUTH

#Select one cluster
# clus <- kmeans(pks[[mag_of_interest]], centers = 3)
# sorted_clusters=sort(clus$size,index.return=TRUE,decreasing = TRUE)$ix
#
# rows=which(clus$cluster==sorted_clusters[clus_num])
# for(attr in attributes(pks)$names)
# {
#   if(is.null(dim(pks[[attr]])))
#   {
#     if(attr!="mass")
#     {
#       pks[[attr]]=pks[[attr]][1]
#     }
#   }
#   else
#   {
#     pks[[attr]]=pks[[attr]][rows,]
#   }
# }

# #First page of metadata [File name, mean image, matrix formula, adducts list, base forms, cor_threshold]
# text=NULL
# text=append(text,as.character(Sys.time()))
# text=append(text,"###################################################################")
# text=append(text,"IMAGE INFORMATION")
# text=append(text,paste(strwrap(paste("- File_name:",pks$names[1])),collapse = "\n"))
# text=append(text,paste(strwrap(paste("- Number of peaks:",length(pks$mass))),collapse = "\n"))
# text=append(text,paste(strwrap(paste("- Number of pixels:",pks$numPixels[1])),collapse = "\n"))
# text=append(text,paste(strwrap(paste("- Mass Range: [",min(pks$mass),", ", max(pks$mass),"]")),collapse = "\n"))
# #text=append(text,paste(strwrap(paste("- Masses:",paste(pks$mass,collapse="; "))),collapse = "\n"))
# text=append(text,"###################################################################")
# text=append(text,"MATRIX INFORMATION")
# text=append(text,paste(strwrap(paste("- Matrix formula:",matrix_formula)),collapse = "\n"))
# text=append(text,paste(strwrap(paste("- Adducts list:",paste(adducts_list,collapse="; "))),collapse = "\n"))
# text=append(text,paste(strwrap(paste("- Base forms:",paste(base_forms,collapse="; "))),collapse = "\n"))
# text=append(text,paste(strwrap(paste("- Correlation threshold:",cor_threshold)),collapse = "\n"))
# text=append(text,"###################################################################")
#
# text=paste(text,collapse = "\n")

# MALDI_resolution=cbind(c(1040.189125,1295.508274,1342.789598,1607.565012,2089.834515,2468.085106,3148.93617,4548.463357),c(26012.14575,34514.17004,36437.24696,41497.97571,44939.27126,44534.41296,42510.12146,37044.53441))
# dimnames(MALDI_resolution)[[2]]=c("m/z","R")
