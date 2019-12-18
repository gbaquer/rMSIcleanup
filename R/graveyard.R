#########################################################################
#
#     GRAVEYARD
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


#' #' Compute S3
#' #'
#' #' Returns the mass tolerance similarity score (S3)
#' #'
#' #' @param calculated_mass Calculated mass vector
#' #' @param experimental_mass Experimental mass vector
#' #' @param similarity_method Similarity method to be used
#' #'
#' #' @return Mass tolerance score S3
#' #'
#' #'
#' compute_s3 <- function (calculated_mass,experimental_mass,similarity_method) {
#'   #Compute S3
#'   # norm_factor=max(experimental_mass)
#'   # s3=exponential_decay_similarity(diff(calculated_mass)*1e4/norm_factor,diff(experimental_mass)*1e4/norm_factor,method=similarity_method,normalize = F)
#'   # if(length(calculated_mass)==1)
#'   # {
#'   #   s3=1
#'   # }
#'   # if(is.na(s3))
#'   #   s3=0
#'   s3=0
#'
#'   return(s3)
#'
#' }

# #4.PRINT SUMMARY REPORT
# gt=NULL#pks$mass[which(results$gt)]
# gt_cluster_names=NULL#results$cluster_names[which(results$gt)]
# not_gt=setdiff(pks$mass,gt)
# include_summary=F
# if(include_summary)
# {
#   #Print images of the ground truth
#   if(length(gt)>0)
#   {
#     plts=list()
#     gt_index=match(gt,pks$mass)
#     for(i in 1:length(gt))
#     {
#       plts=c(plts,list(ggplot_peak_image(pks,pks[[mag_of_interest]][,gt_index[i]],paste("m/z",round(gt[i],pkg_opt("round_digits")),"(",gt_cluster_names[i],")"),chosen=T)))
#       if(i%%length(page_layout)==0||i==length(gt))
#       {
#         grid.arrange(grobs=plts,nrow=nrow(page_layout),ncol=ncol(page_layout))
#         plts=list()
#       }
#     }
#   }
#   #Print images of the not ground truth
#   if(length(not_gt)>0)
#   {
#     plts=list()
#     not_gt_index=match(not_gt,pks$mass)
#     for(i in 1:length(not_gt))
#     {
#       plts=c(plts,list(ggplot_peak_image(pks,pks[[mag_of_interest]][,not_gt_index[i]],paste("m/z",round(not_gt[i],pkg_opt("round_digits"))),chosen=F)))
#       if(i%%length(page_layout)==0||i==length(not_gt))
#       {
#         grid.arrange(grobs=plts,nrow=nrow(page_layout),ncol=ncol(page_layout))
#         plts=list()
#       }
#     }
#   }
# }



#' #' Compute scores
#' #'
#' #' Returns several scores for performance assessment of a given binary classification result.
#' #'
#' #' @param gt Ground truth: Positive values
#' #' @param pos Classified positives
#' #' @param neg Classified negatives
#' #'
#' #' @return List with all the computed scores. It includes F1 score and Brier score. As well as the number of
#' #'
#' #'
#' #'
#' compute_scores <- function (gt,pos,neg) {
#'   #compute tp, tn, fp, fn
#'   tp_list=intersect(gt,pos)
#'   fp_list=setdiff(pos,tp_list)
#'   fn_list=intersect(gt,neg)
#'   tn_list=setdiff(neg,fn_list)
#'
#'   # tp_list=intersect(pos,gt) #true positive
#'   # fp_list=which(!is.element(pos,tp)) #setdiff(pos,tp) #false positive
#'   # fn_list=intersect(neg,gt) #false negative
#'   # tn_list=which(!is.element(neg,fn))#setdiff(neg,fn) #true negative
#'
#'   tp=length(tp_list)
#'   fp=length(fp_list)
#'   fn=length(fn_list)
#'   tn=length(tn_list)
#'
#'   p=tp/(tp+fp) #precision
#'   r=tp/(tp+fn) #recall
#'
#'   scores=list()
#'   scores$f1= 2*(r*p)/(r+p) #F1 score
#'   scores$b= (fp+fn)/(length(pos)+length(neg)) #Brier score
#'   scores$tp=tp
#'   scores$fp=fp
#'   scores$fn=fn
#'   scores$tn=tn
#'   scores$p=p
#'   scores$r=r
#'
#'   if(pkg_opt()$verbose_level<=-1)
#'   {
#'     print(paste("tp:",length(tp),"tn:",length(tn),"fp:",length(fp),"fn:",length(fn)))
#'     print(paste("p:",p,"r:",r))
#'     print(paste("F1 score:",scores$f1,"BS:",scores$b))
#'   }
#'   return(scores)
#' }

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

# print_experiment_names <- function (matrix_formula, base_dirs=c("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/Ag Software Test 1","/home/gbaquer/msidata/Ag Software Test 1"),
#                                    s1_threshold=0.80,s2_threshold=0.80, s3_threshold=0.7, similarity_method="euclidean", correlation_method="pearson", experiment_name="output",
#                                    MALDI_resolution=20000, tol_mode="scans",tol_ppm=200e-6,tol_scans=4,
#                                    mag_of_interest="intensity",normalization="None",
#                                    max_multi=10, add_list=NULL, sub_list=NULL, isobaric_detection=T,
#                                    save_results=T,generate_pdf=T,default_page_layout=NULL,include_summary=F,dataset_indices=NULL) {
#   #0. Prepare directories
#   base_dir=base_dirs[which(dir.exists(base_dirs))][1]
#
#   images_dir=paste(base_dir,"/images",sep="")
#   output_dir=paste(base_dir,"/output",sep="")
#
#   #1. Generate experiment metadata file
#   #[PENDING]
#
#   #2. Run experiment
#   full_spectrum_name_list=list.files(images_dir,pattern = "*-proc.tar",recursive = T)
#   pks_name_list=list.files(images_dir,pattern = "*mergeddata-peaks*",recursive = T,full.names = T)
#   subfolder_list=unlist(lapply(strsplit(pks_name_list,'/'),function(x) paste(x[1:length(x)-1],collapse="/")))
#   num_files=length(pks_name_list)
#
#   results=list(meta=list(s1_threshold=s1_threshold,s2_threshold=s2_threshold,s3_threshold=s3_threshold,file_names=NULL),data=list())
#
#   j=1
#   for(i in 1:num_files)
#   {
#     if(is.null(dataset_indices) || is.element(i,dataset_indices))
#     {
#       pks_name=pks_name_list[i]
#       subfolder=subfolder_list[i]
#       pks=rMSIproc::LoadPeakMatrix(pks_name)
#
#       pks_i=1
#       for(name in pks$names)
#       {
#         full_spectrum_name=paste(subfolder,"/",unlist(strsplit(name,".",fixed=T))[1],"-proc.tar",sep="")
#
#         if(file.exists(full_spectrum_name))
#         {
#           cat(paste(j,":",full_spectrum_name,"\n"))
#           j=j+1
#         }
#         pks_i=pks_i+1
#       }
#     }
#   }
#
# }


#Disabled because of XML incompatibilities

# #' Export to mmass
# #'
# #' Exports the results to a text file which can be read by mmass for easy interpretation and validation of the results.
# #'
# #' @param pks_Matrix Peak matrix
# #' @param matrix_annotation Vector determining whether each peak in the spectrum corresponds to the matrix or the tissue
# #' @param metadata List containing metadata on the matrix annotation process
# #'
# #' @return None
# #'
# #'
# #'
# export_mmass <- function (pks_Matrix,matrix_annotation=rep(1,length(pks_Matrix$mass)),metadata=FALSE) {
#   #IMPORT PYTHON MODULES
#   struct=reticulate::import("struct")
#   zlib=reticulate::import("zlib")
#   base64=reticulate::import("base64")
#
#   #Cluster image
#   centers=3
#   clus <- kmeans((pks_Matrix$intensity/pks_Matrix$normalizations$TIC), centers = centers)
#
#   #Plot cluster image
#   dev.new()
#   rMSIproc::plotClusterImage(pks_Matrix,clus$cluster)
#
#   for(i in 0:centers )
#   {
#     if(i==0){
#       pixels=1:pks_Matrix$numPixels
#       file_suffix="whole_image"
#     }
#     else{
#       pixels=which(clus$cluster==i)
#       file_suffix=paste("cluster_",i,sep="")
#     }
#
#     #Calculate mean spectra
#     mean_spectra=apply(pks_Matrix$intensity[pixels,],2,mean)
#     masses=lapply(pks_Matrix$mass,function(x) x)
#     #Create table
#     export_table=data.frame(mass=pks_Matrix$mass,intensity=mean_spectra)
#     #Write txt file
#     write.table(export_table, file = paste("output/complete_spectrum_",file_suffix,".txt",sep=""), sep = "\t", row.names=FALSE, col.names=FALSE)
#
#     #CONVERT TO BINARY WITH THE PROPER UTF ENCODING
#     intArray=""
#     mzArray=""
#     removed=0
#     matrix_peaks=which(matrix_annotation>0)
#     for (i in matrix_peaks) {
#       tryCatch({
#         tmp1=iconv(struct$pack("f",as.single(mean_spectra[i])))
#         tmp2=iconv(struct$pack("f",as.single(masses[i])))
#       },
#       error=function(cond) {
#         print("Skipped one peak")
#         print(cond)
#         tmp1=0
#         tmp2=0
#         removed=removed+1
#       },
#       finally = {
#         intArray=paste(intArray,tmp1,sep = "")
#         mzArray=paste(mzArray,tmp2,sep = "")
#       }
#       )
#     }
#
#     # COMPRESSION [IT FAILS PROBABLY DUE TO ENCODING ISSUES AGAIN]
#     # intArray=zlib$compress(intArray)
#     # mzArray=zlib$compress(mzArray)
#
#     #BASE 64 ENCODING
#     mzArray = base64$b64encode(mzArray)
#     intArray = base64$b64encode(intArray)
#
#     #WRITE mSD FILE
#     xml <- XML::xmlTree("mSD_file")
#     # names(xml)
#     xml$addNode("mSD", close=FALSE, attrs=c(version="2.2"))
#     #Description
#     xml$addNode("description", close=FALSE)
#     xml$addNode("title","MATRIX_ANNOTATION_REPORT")
#     xml$addNode("date", attrs=c(value=Sys.time()))
#     xml$addNode("operator", attrs=c(value=""))
#     xml$addNode("contact", attrs=c(value="Gerard Baquer Gomez (gerard.baquer@urv.cat)"))
#     xml$addNode("institution", attrs=c(value="MIL@B (URV)"))
#     xml$addNode("instrument", attrs=c(value=""))
#     xml$addNode("notes","THIS SPACE IS RESERVED FOR NOTES")
#     xml$closeTag()
#
#     #Spectrum
#     xml$addNode("spectrum",attrs=c(points=(length(mean_spectra)-removed)), close=FALSE)
#     #xml$addNode("spectrum",attrs=c(points=1), close=FALSE)
#     xml$addNode("mzArray", mzArray, attrs=c(precision="32", endian="little"))
#     xml$addNode("intArray",intArray, attrs=c(precision="32", endian="little"))
#     xml$closeTag()
#
#     #Peaklist
#     #[THE PEAK LIST WILL PROBABLY INCLUDE THE PEAK MATRIX WHEN THE FUNCTION IS UPGRADED TO INCLUDE THE FULL SPECTRUM]
#
#     #Annotations
#     xml$addNode("annotations", close=FALSE)
#     for (i in matrix_peaks) {
#       xml$addNode("annotation", matrix_annotation[i], attrs=c(peakMZ=masses[i], peakIntensity=mean_spectra[i]))
#     }
#     xml$closeTag()
#
#     #Save document
#     xml$closeTag()
#     saveXML(xml,paste("output/matrix_peaks_",file_suffix,".msd",sep=""),prefix='<?xml version="1.0" encoding="utf-8" ?>\n')
#   }
# }
