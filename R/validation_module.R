 #########################################################################
#
#     VALIDATION MODULE
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

# HERE

#' Generate ground truth
#'
#' Given the matrix formula and the matrix peaks present in the image it returns the ground truth defined as the list of peaks present in the image that correspond to the matrix.
#'
#' @param matrix_formula String giving the chemical formula of the matrix in the enviPat notation
#' @param pks Peak Matrix Image produced using rMSIproc
#' @param full_spectrum Full spectrum before performing peak picking. It is used to give a higher degree of confidence to the S1 and S2 computation.
#' @param s1_threshold Correlation between the theoretical and the real spectral pattern above which a given cluster is considered to be present
#' @param s2_threshold Correlation between the spatial images of the peak in a cluster above which the peaks are included in the gold truth (gt)
#' @param s3_threshold Correlation between the spatial images of the peak in a cluster above which the peaks are included in the gold truth (gt)
#' @param MALDI_resolution MALDI resolution that is used to merge nearby peaks in-silico as the equipment would in real life
#' @param tol_mode String determining the tolerance mode to be used when comparing the calculated masses with the experimental ones. "ppm": relative tolerance with respect to the calculated one in parts per million. "scans": Number of scans or datapoints present in the full spectrum -It is only applicable if full_spectrum is provided.
#' @param tol_ppm Tolerance in parts per million. Only used if tol_mode="ppm".
#' @param tol_scans Tolerance in number of scans. Only used if tol_mode="scans".
#' @param mag_of_interest Magnitude of interest to use when performing the study. It can be "intensity" or "area".
#' @param normalization String indicating the normalization technique to use. Possible values: "None", "TIC","RMS", "MAX" or "AcqTic"
#' @param max_multi Maximum Cluster Multiplication. For each adduct -Adduct- it will generate a set of base formulas as follows: M+Adduct, 2M+Adduct, 3M+Adduct ... max_multi*M + Adduct
#' @param add_list List of adducts to be added to the matrix formula in the format defined by enviPat. Example: c("H1","Na1","K1")
#' @param sub_list List of compounds to be substracted to the matrix formula in the format defined by enviPat. Example: c("H1",H2O1")
#' @param generate_pdf Boolean indicating whether to generate a pdf report or not
#' @param generate_pdf Boolean indicating whether to generate a TIFF figures or not
#' @param default_page_layout Page layout to be used in the pdf plotting.
#' @param folder Name of the folder in which to store the pdf report
#' @param similarity_method Similarity method to be used in the computation of the distance
#' @param pks_i Index of the peak matrix
#' @param plot_type String indicating the type of plot desired. Can have two possible values "debug" or "poster"
#' @param include_summary Boolean value indicating if the summary is to be included in the pd
#'
#' @return Ground Truth: List of masses available in the image that correspond to the matrix.
#'
#' @export
generate_gt <- function (matrix_formula,pks,full_spectrum=NULL, folder="output/",
                         s_threshold=0.65,s1_threshold=0.80,s2_threshold=0.80, s3_threshold=0.7, similarity_method="euclidean",correlation_method="pearson",
                         MALDI_resolution=20000, tol_mode="scans",tol_ppm=200e-6,tol_scans=4,
                         mag_of_interest="intensity",normalization="None",pks_i=1,
                         min_multi=1,max_multi=10, add_list=NULL, sub_list=NULL, max_charge=1, isobaric_detection=T,
                         generate_pdf=T,generate_figures=F,default_page_layout=NULL,plot_type="debug",include_summary=F) {
  #SECTION -1 :: Input validation
  #Adjust tolerance to ppm if no full_spectrum is provided
  if(is.null(full_spectrum))
    tol_mode="ppm"
  #Page layout
  if(is.null(default_page_layout))
  {
    # if(plot_type=="poster")
    #   default_page_layout=rbind(c(1,1,1,1),
    #                           c(1,1,1,1),
    #                           c(2,2,2,2),
    #                           c(3,4,5,6),
    #                           c(7,8,9,10),
    #                           c(11,12,13,14))
    # else
    #   default_page_layout=rbind(c(1,1,2,2),
    #                             c(1,1,2,2),
    #                             c(3,4,5,6),
    #                             c(7,8,9,10),
    #                             c(11,12,13,14),
    #                             c(15,16,17,18))
    default_page_layout=rbind(c(1,1,4,5,6),
                              c(1,1,7,8,9),
                              c(2,2,10,11,12),
                              c(3,3,13,14,15),
                              c(3,3,16,17,18),
                              c(3,3,19,20,21))

  }
  #SECTION 0 :: Preprocessing

  # Select first image if there are multiple
  #pks=get_one_peakMatrix(pks,pks_i)

  # Normalize
  if(normalization!="None")
    pks[[mag_of_interest]]=pks[[mag_of_interest]]/pks$normalizations[[normalization]]

  #0. Load data
  isotopes=NULL
  data("isotopes", package = "enviPat", envir = environment())
  adducts=NULL
  data("adducts", package = "enviPat", envir=environment())

  #1. Determine adducts depending on matrix_formula

  #For more complex matrix formulas adducts could be loaded from the library
  #adducts_formula=adducts$Formula_add[1:4]

  #2. [Not implemented yet] Assess which mass range is needed

  #3.Generate list of possible chemical formulas
  if(!is.null(add_list))
    matrix_formula=append(matrix_formula,unlist(lapply(add_list,function(x) enviPat::mergeform(matrix_formula,x))))
  base_forms=NULL
  for(i in min_multi:max_multi)
  {
    base_forms=append(base_forms,multiform(matrix_formula,i))
  }
  if(generate_figures)
    base_forms=c("Ag6")
  # if(!is.null(add_list))
  #   base_forms=append(base_forms,unlist(lapply(add_list,function(x) enviPat::mergeform(base_forms,x))))
  # if(!is.null(sub_list))
  #   base_forms=append(base_forms,unlist(lapply(sub_list,function(x) enviPat::subform(base_forms,x)[which(enviPat::check_ded(base_forms,x)=="FALSE")])))

  sorted_ix=sort(check_chemform(isotopes,base_forms)$monoisotopic_mass,index.return=T)$ix
  base_forms=base_forms[sorted_ix]

  #4. Open pdf report
  if(generate_pdf)
  {
    #Generate pdf name [pks$name_001]
    pdf_file=generate_file_name(pks$names[1],folder = folder)
    #open pdf file
    a4_width=8.27
    a4_height=11.69
    pdf(pdf_file,paper='a4',width=a4_width,height=a4_height)

    #First page of metadata [File name, mean image, matrix formula, adducts list, base forms, S1_threshold]
    text=""
    text=add_entry(text,"#")
    text=add_entry(text,"- Package Version:",paste(packageVersion("rMSIcleanup"),collapse = "."))
    text=add_entry(text,"- Time:",as.character(Sys.time()))
    text=add_entry(text,"#")
    text=add_entry(text,"IMAGE INFORMATION")
    text=add_entry(text,"- Peak matrix:",pks$names[1])
    if(is.null(full_spectrum))
      text=add_entry(text,"- Full spectrum: [NOT PROVIDED]")
    else
      text=add_entry(text,"- Full spectrum:",full_spectrum$names[1])

    text=add_entry(text,"- Number of peaks:",length(pks$mass))
    text=add_entry(text,"- Number of pixels:",pks$numPixels[1])
    text=add_entry(text,"- Mass Range: [",min(pks$mass),", ", max(pks$mass),"]")
    text=add_entry(text,"#")
    text=add_entry(text,"MATRIX INFORMATION")
    text=add_entry(text,"- Matrix formula:",matrix_formula)
    text=add_entry(text,"- Add list:",add_list)
    text=add_entry(text,"- Substract list:",sub_list)
    text=add_entry(text,"- Maximum cluster multiplication:",max_multi)
    text=add_entry(text,"- Base forms:",base_forms)
    text=add_entry(text,"#")
    text=add_entry(text,"PROCESSING INFORMATION")
    text=add_entry(text,"- S1 threshold:",s1_threshold)
    text=add_entry(text,"- S2 threshold:",s2_threshold)
    text=add_entry(text,"- Similarity method:",similarity_method)
    text=add_entry(text,"- Magnitude of interest:",mag_of_interest)
    text=add_entry(text,"- Tolerance mode:",tol_mode)
    if(tol_mode=="ppm")
      text=add_entry(text,"- Tolerance ppm:",tol_ppm)
    else
      text=add_entry(text,"- Tolerance scans:",tol_scans)
    text=add_entry(text,"#")

    text=paste(text,collapse = "\n")
    print(plot_text(text))
  }

  #5.Generate pattern list with enviPat

  patterns_out=NULL
  multiplier=1
  while(TRUE)
  {
    #Multiply
    forms=multiform(base_forms,multiplier)

    #Multiple charge [SOLVE: THE CURRENT IMPLEMENTATION IS MESSING UP WITH THE ADDUCT CLUSTERS]
    # charge=rep(1:max_charge, each=max_charge)
    # forms=rep(forms,max_charge)

    #Check
    checked=check_chemform(isotopes,forms)

    #Remove incorrect formulas
    checked=checked[which(checked$warning==F),]

    #Lower mass limit
    checked=checked[which(checked$monoisotopic_mass>=min(pks$mass)),]
    # checked=checked[which(checked$monoisotopic_mass/charge>=min(pks$mass)),]
    # charge=charge[which(checked$monoisotopic_mass/charge>=min(pks$mass))]
    if(nrow(checked)==0)
      next;

    #Upper mass limit
    checked=checked[which(checked$monoisotopic_mass<=max(pks$mass)),]
    # checked=checked[which(checked$monoisotopic_mass/charge<=max(pks$mass)),]
    # charge=charge[which(checked$monoisotopic_mass/charge<=max(pks$mass))]
    if(nrow(checked)==0)
      break;


    #Compute theoretical patterns
    patterns=quiet(isowrap(isotopes,checked,resmass = FALSE,resolution = MALDI_resolution))
    # patterns=quiet(isowrap(isotopes,checked,charge=charge,resmass = FALSE,resolution = MALDI_resolution))

    #Append to final list
    for(i in 1:length(patterns))
    {
      if(!is.element(attributes(patterns)$names[i],patterns_out$cluster))
      {
        patterns_out$mass=append(patterns_out$mass,patterns[[i]][,1]) #store masses
        patterns_out[[mag_of_interest]]=append(patterns_out[[mag_of_interest]],patterns[[i]][,2]) #store intensities
        patterns_out$cluster=append(patterns_out$cluster,rep(attributes(patterns)$names[i],length(patterns[[i]][,1]))) #store cluster names
        # patterns_out$cluster=append(patterns_out$cluster,rep(paste(attributes(patterns)$names[i],"_",charge[i],sep=""),length(patterns[[i]][,1]))) #store cluster names
      }
    }
    multiplier=multiplier+1
    if(multiplier>max_multi)
      break;
  }

  #6. Determine s1 and s2 scores for each calculated cluster
  mean_image=apply(pks[[mag_of_interest]],2,mean)
  sd_image=apply(pks[[mag_of_interest]],2,sd)/mean_image
  clusters=unique(patterns_out$cluster)

  #Initialize output results list
  num_pattern_elements=length(patterns_out$mass)

  patterns_out$experimental_mass=rep(NA,num_pattern_elements)
  patterns_out$s1_scores=rep(NA,num_pattern_elements)
  patterns_out$s2_scores=rep(NA,num_pattern_elements)
  patterns_out$present=rep(F,num_pattern_elements)

  num_peaks_in_matrix = length(pks$mass)

  results=list(
    s1_scores=rep(NA,num_peaks_in_matrix),
    s2_scores=rep(NA,num_peaks_in_matrix),
    s3_scores=rep(NA,num_peaks_in_matrix),
    cluster_names=rep(NA,num_peaks_in_matrix),
    gt=rep(F,num_peaks_in_matrix),
    count=rep(0,num_peaks_in_matrix),
    patterns_out=patterns_out)
  for(cluster in clusters)
  {
    # Calculated cluster
    calculated_index=which(patterns_out$cluster==cluster)
    calculated_mass=patterns_out$mass[calculated_index]
    calculated_magnitude=patterns_out[[mag_of_interest]][calculated_index]


    # Experimental cluster
    experimental_index=get_closest_peak(calculated_mass,pks$mass)
    experimental_mass=pks$mass[experimental_index]
    experimental_magnitude=mean_image[experimental_index]
    experimental_sd=sd_image[experimental_index]

    # Number of peaks in cluster
    num_peaks=length(calculated_index)

    #Determine mass relative error
    rel_error=(experimental_mass-calculated_mass)/calculated_mass
    if(tol_mode=="ppm")
    {
      index_pks=which(abs(rel_error)<=tol_ppm) #Peaks to be taken from the peak matrix
      index_full_spectrum=which(abs(rel_error)>tol_ppm) #Peaks to be taken from the full spectrum
    }
    else
    {
      index_pks=which(unlist(lapply(1:num_peaks,function(i)is_within_scan_tol(experimental_mass[i],calculated_mass[i],full_spectrum$mass,tol_scans))))#Peaks to be taken from the peak matrix
      index_full_spectrum=setdiff(1:num_peaks,index_pks) #Peaks to be taken from the full spectrum
    }
    rel_error[index_full_spectrum]=0
    #Generate images for each peak in the cluster
    final_image=NULL
    for(i in 1:num_peaks)
    {
      if(is.element(i,index_full_spectrum))
      {
        if(!is.null(full_spectrum))
        {
          #Get image from full_spectrum (Processed spectrum before peak picking)
          mass=calculated_mass[i]
          if(tol_mode=="ppm")
            data=rMSI::loadImageSliceFromMass(full_spectrum, mass,mass*tol_ppm)$data
          else
          {
            col=which.min(abs(full_spectrum$mass-mass))
            cols=(col-tol_scans/10):(col+tol_scans/10)
            data=rMSI::loadImageSliceFromCols(full_spectrum,cols)
          }

          if(normalization!="None")
            data=data/pks$normalizations[[normalization]]
          max_col=which.max(apply(data,2,mean))
          data_picture=data[,max_col]
          experimental_magnitude[i]=mean(data_picture) #Update the experimental_magnitude (the previous value corresponds to the closest peak in the peak matrix which in this case is further than the specified tolerance and should thus be dropped)
        }
        else
        {
          #Set image to NA
          data_picture=rep(NA,pks$numPixels[1])
          experimental_magnitude[i]=0
        }
      }
      else
      {
        #Get image from pks (Peak matrix)
        data_picture=pks[[mag_of_interest]][,experimental_index[i]]
      }
      final_image=cbind(final_image,data_picture) #Update list of final images
    }
    colnames(final_image)<-NULL #Remove column names

    image_correl=cor(final_image,method=correlation_method)
    image_correl[which(is.na(image_correl))]=0

    #Compute similarity scores using all peaks in the cluster
    s1=compute_s1(calculated_magnitude,experimental_magnitude,similarity_method)
    s2=compute_s2(calculated_magnitude,image_correl)
    s3=compute_s3(calculated_mass[index_pks],experimental_mass[index_pks],similarity_method)

    #Compute similarity scores for each possible merged clusters
    index_selected=1:num_peaks
    s1_clusters=NULL
    s2_clusters=NULL
    s3_clusters=NULL

    a=kmeans(image_correl,centers=min(2,nrow(image_correl)-1))
    b=kmeans(rel_error,centers=min(2,nrow(image_correl)-1,length(unique(rel_error))))
    ion_clusters=paste(a$clus,b$clus)

    classification_record=rbind(rep(1,num_peaks))
    classes_list=rbind(c(1,s1,s2,s3))
    # if(isobaric_detection&((s1*s2)<s_threshold))
    if(isobaric_detection&(s1>s1_threshold)&(s2>s2_threshold) )
    {
      found=F
      while(!found & any(!is.na(classification_record[1,])))
      {
        parent_classes=classification_record[1,]
        child_classes=rep(NA,num_peaks)
        child_classes_priority=rep(NA,num_peaks)
        for(parent_class in unique(parent_classes[which(!is.na(parent_classes))]))
        {
          index_class=which(parent_classes==parent_class)
          if(length(index_class)>2 & NROW(unique(image_correl[index_class,index_class]))>2 & length(intersect(index_class,index_pks))>0)
            child_classes[index_class]=max(parent_classes,na.rm = T)+(parent_class-min(parent_classes,na.rm = T))+kmeans(image_correl[index_class,index_class],centers=2)$clus
        }

        classification_record=rbind(child_classes,classification_record)

        for(parent_class in unique(child_classes[which(!is.na(child_classes))]))
        {
          index_class=which(child_classes==parent_class)
          child_classes_priority[index_class]=sum(calculated_magnitude[index_class])
        }

        # child_classes=child_classes[which(!is.na(child_classes))]
        #for(child_class in unique(child_classes[which(!is.na(child_classes))]))
        for(child_class_priority in unique(sort(child_classes_priority,decreasing = T)))
        {
          index_class=which(child_classes_priority==child_class_priority)
          index_not_class=which(child_classes_priority!=child_class_priority)
          child_class=child_classes[index_class][1]
          #Skip if cluster only contains 1 peak
          #Skip if other cluster peaks are not higher than 90% of the theoretical
          # if((length(index_class)<=1) | (experimental_magnitude[index_not_class]/max(experimental_magnitude[index_class])<0.9*calculated_magnitude[index_not_class]))
          if((length(index_class)<=1)|length(index_class)/num_peaks<0.5)
          {
            #classification_record[1,index_class]=NA
            next
          }
          #Compute scores
          s1_tmp=compute_s1(calculated_magnitude[index_class],experimental_magnitude[index_class],similarity_method)
          s2_tmp=compute_s2(calculated_magnitude[index_class],image_correl[index_class,index_class])
          s3_tmp=compute_s3(calculated_mass[intersect(index_pks,index_class)],experimental_mass[intersect(index_pks,index_class)],similarity_method)

          #Strore results
          classes_list=rbind(classes_list,c(child_class,s1_tmp,s2_tmp,s3_tmp))
          #Finish if similarity conditions met
          if((s1_tmp>s1_threshold)&(s2_tmp>s2_threshold))
          {
            s1=s1_tmp
            s2=s2_tmp
            s3=s3_tmp
            index_selected=index_class
            found=T
            break
          }
          # #Update pareto front
          # if((s1_tmp>=s1) & (s2_tmp>=s2) & (s3_tmp>=s3))
          # {
          #   s1=s1_tmp
          #   s2=s2_tmp
          #   s3=s3_tmp
          #   index_selected=index_class
          # }
        }
      }

      # s1_clusters=rep(0,length(a$size))
      # s2_clusters=rep(0,length(a$size))
      # s3_clusters=rep(0,length(a$size))
      # for(i in a$clus)
      # {
      #   index_cluster=which(a$clus==i)
      #   if(length(index_cluster)>1)
      #   {
      #     s1_clusters[i]=compute_s1(calculated_magnitude[index_cluster],experimental_magnitude[index_cluster],similarity_method)
      #     s2_clusters[i]=compute_s2(calculated_magnitude[index_cluster],image_correl[index_cluster,index_cluster])
      #     s3_clusters[i]=compute_s3(calculated_mass[intersect(index_pks,index_cluster)],experimental_mass[intersect(index_pks,index_cluster)],similarity_method)
      #   }
      #   if(s1_clusters[i]>s1_threshold&s2_clusters[i]>s2_threshold&s3_clusters[i]>s3_threshold)
      #   {
      #     s1=s1_clusters[i]
      #     s2=s2_clusters[i]
      #     s3=s3_clusters[i]
      #     index_selected=index_cluster
      #   }
      # }
      pareto_front=NA
      # pareto_front=get_pareto_front(append(s1,s1_clusters),append(s2,s2_clusters),append(s3,s3_clusters))
      # s1=pareto_front[1,1]
      # s2=pareto_front[1,2]
      # s3=pareto_front[1,3]
      # if(pareto_front)
      #
      # selected_row=which(pareto_front[,1]>s1_threshold & pareto_front[,2]>s2_threshold & pareto_front[,3]>s3_threshold)[1]
      #
      # selected_row=as.numeric(row.names(r))[selected_row]
    }

    # #Compute S1 using all peaks
    # s1_all=exponential_decay_similarity(calculated_magnitude,experimental_magnitude,method=similarity_method)
    # if(num_peaks==1)
    #   s1=1
    # #Compute S1 using only the peak matrix
    # if(length(index_pks)>1)
    # {
    #   s1_pks=exponential_decay_similarity(calculated_magnitude[index_pks],experimental_magnitude[index_pks],method=similarity_method)
    # }
    # else
    #   s1_pks=NA
    #
    # #Compute S1
    # s1=max(s1_all,s1_pks,na.rm=T)
    #
    # #Compute S2 using all peaks
    # image_correl=cor(final_image)
    # image_correl[which(is.na(image_correl))]=0
    # s2_individual = apply(image_correl,2,weighted.mean,w=calculated_magnitude)
    # s2_individual[which(is.na(s2_individual))]=0
    # s2_all=weighted.mean(s2_individual,calculated_magnitude)
    #
    # #Compute S2 using only the peak matrix
    # if(length(index_pks)>1)
    # {
    #   s2_individual_pks = apply(image_correl[index_pks,index_pks],2,weighted.mean,w=calculated_magnitude[index_pks])
    #   s2_individual_pks[which(is.na(s2_individual_pks))]=0
    #   s2_pks=weighted.mean(s2_individual_pks,calculated_magnitude[index_pks])
    # }
    # else
    #   s2_pks=NA
    #
    # #Compute S2
    # s2=max(s2_all,s2_pks,na.rm=T)
    #
    # #Compute S3
    # a=calculated_mass[index_pks]
    # b=experimental_mass[index_pks]
    # s3=rMSIcleanup::exponential_decay_similarity(diff(a)*1e4/max(b),diff(b)*1e4/max(b),method=similarity_method,normalize = F)
    # if(length(index_pks)==1)
    # {
    #   s3=1
    # }
    #s3=exponential_decay_similarity(calculated_mass[index_pks]/calculated_mass[index_pks],experimental_mass[index_pks]/calculated_mass[index_pks],method=similarity_method)

    #Check if there is an overlap

    # if(s1<s1_threshold & s2<s2_threshold)
    # {
    #   ion_clusters=kmeans(image_correl,centers=min(2,nrow(image_correl)-1))$clus
    #   #Recompute similarity scores
    #   for(i in ion_clusters$clus)
    #   {
    #     ions=which(ion_clusters$clus==i)
    #     num_peaks_tmp=length(ions)
    #     #Compute S1 using all peaks
    #     s1_all=exponential_decay_similarity(calculated_magnitude[ions],experimental_magnitude[ions],method=similarity_method)
    #     if(num_peaks_tmp==1)
    #       s1=1
    #     #Compute S1 using only the peak matrix
    #     if(length(index_pks)>1)
    #     {
    #       s1_pks=exponential_decay_similarity(calculated_magnitude[intersect(index_pks,ions)],experimental_magnitude[intersect(index_pks,ions)],method=similarity_method)
    #     }
    #     else
    #       s1_pks=NA
    #
    #     #Compute S1
    #     s1=max(s1_all,s1_pks,na.rm=T)
    #
    #     #Compute S2 using all peaks
    #     s2_individual_tmp = apply(image_correl[ions,ions],2,weighted.mean,w=calculated_magnitude[ions])
    #     s2_individual_tmp[which(is.na(s2_individual_tmp))]=0
    #     s2_all=weighted.mean(s2_individual_tmp,calculated_magnitude[ions])
    #
    #     #Compute S2 using only the peak matrix
    #     if(length(index_pks)>1)
    #     {
    #       s2_individual_pks_tmp = apply(image_correl[intersect(index_pks,ions),intersect(index_pks,ions)],2,weighted.mean,w=calculated_magnitude[intersect(index_pks,ions)])
    #       s2_individual_pks_tmp[which(is.na(s2_individual_pks_tmp))]=0
    #       s2_pks=weighted.mean(s2_individual_pks_tmp,calculated_magnitude[intersect(index_pks,ions)])
    #     }
    #     else
    #       s2_pks=NA
    #
    #     #Compute S2
    #     s2=max(s2_all,s2_pks,na.rm=T)
    #
    #     #Compute S3
    #     a=calculated_mass[intersect(index_pks,ions)]
    #     b=experimental_mass[intersect(index_pks,ions)]
    #     s3=rMSIcleanup::exponential_decay_similarity(diff(a)*1e4/max(b),diff(b)*1e4/max(b),method=similarity_method,normalize = F)
    #     if(length(intersect(index_pks,ions))==1)
    #       s3=1
    #
    #     if(s1>s1_threshold&s2>s2_threshold&s3>s3_threshold)
    #     {
    #       break;
    #     }
    #   }
    # }


    #Choose which peaks belong in the ground truth (gt)
    chosen=rep(F,num_peaks)
    #chosen[intersect(index_pks,index_selected)]=(!is.na(s1)&!is.na(s2)&!is.na(s3))&(s1*s2>s_threshold)
    # chosen[index_selected]=(!is.na(s1)&!is.na(s2)&!is.na(s3))&(s1*s2>s_threshold)
    chosen[index_selected]=(!is.na(s1)&!is.na(s2)&!is.na(s3))&(s1>s1_threshold)&(s2>s2_threshold)

    #Adjust magnitude in the NA mode for plotting
    if(is.null(full_spectrum))
      experimental_magnitude[index_full_spectrum]=NA


    #Store results
    results$s1_scores[experimental_index[which(chosen)]]=s1
    results$s2_scores[experimental_index[which(chosen)]]=s2
    results$s3_scores[experimental_index[which(chosen)]]=s3
    results$cluster_names[experimental_index[which(chosen)]]=cluster
    results$gt[experimental_index[which(chosen)]]=T
    results$count[experimental_index[which(chosen)]]=results$count[experimental_index[which(chosen)]]+1

    results$patterns_out$experimental_mass[calculated_index]=experimental_mass
    results$patterns_out$s1_scores[calculated_index]=s1
    results$patterns_out$s2_scores[calculated_index]=s2
    results$patterns_out$s3_scores[calculated_index]=s3
    results$patterns_out$present[calculated_index[index_pks]]=T

    # 7. Generate plots
    if(pkg_opt()$verbose_level<=-1 || generate_pdf || generate_figures)
    {
      #Print progress to console
      print(paste(cluster,s1,s2,s3))
      # print(classification_record)
      # print(classes_list)
      # print(index_selected)
      # print(table(a$cluster,b$cluster))
      # print(rbind(s1_clusters,s2_clusters,s3_clusters))
      # print(pareto_front)

      #Normalize
      if(sum(!is.na(experimental_magnitude))!=0 && max(experimental_magnitude,na.rm=T)!=0)
        norm_experimental_magnitude=experimental_magnitude/max(experimental_magnitude,na.rm=T)
      else
        norm_experimental_magnitude=experimental_magnitude

      if(sum(!is.na(experimental_sd))!=0 && max(experimental_sd,na.rm=T)!=0)
        norm_experimental_sd=experimental_sd/max(experimental_sd,na.rm=T)
      else
        norm_experimental_sd=experimental_sd

      if(max(calculated_magnitude)!=0)
        calculated_magnitude=calculated_magnitude/max(calculated_magnitude)

      #Prepare melted dataframe
      df=data.frame(mass=calculated_mass,calculated_magnitude=calculated_magnitude,experimental_magnitude=norm_experimental_magnitude,experimental_sd=norm_experimental_sd)
      vars=c("calculated_magnitude","experimental_magnitude")
      vars_label=c("Calculated","Experimental")
      melted_df=melt(df, id.vars="mass",measure.vars=vars)
      value=NULL
      variable=NULL
      mass=NULL

      #Spectrum comparison plot

      l=num_peaks*length(vars)
      if(plot_type=="poster")
      {
        label=rep("",l)
        linetype=rep("solid",l)
        legend_position="right"
      }
      else
      {
        label=rep("",l)
        label[(num_peaks+1):(num_peaks*2)]=ion_clusters
        linetype=rep("solid",l)
        linetype[index_full_spectrum+num_peaks]="dotted"
        legend_position="bottom"
        #linetype=append(rep("solid",length(s2_individual)),rep("dotted",length(s2_individual)))
      }
      global_title=paste(cluster,"; S1:",round(s1,2),"; S2:",round(s2,2),"; S3:",round(s3,2),"; MAX:", round(max(experimental_magnitude,na.rm=T),2))

      plt_spectrum= ggplot(melted_df, aes(mass,value,color=variable,ymin=0,ymax=value)) +
                    geom_linerange() + geom_point() + theme_bw() +
        # ggtitle(paste(cluster,"; S1:",round(s1,2),"; S2:",round(s2,2),"; S3:",round(s3,2),"; MAX:", round(max(experimental_magnitude,na.rm=T),2))) +
                    xlab("m/z") + ylab("Scaled Intensity") +
                    scale_colour_discrete(name="",breaks=vars,labels=vars_label)+
                    theme(legend.justification = c(1, 1), legend.position = c(1, 1),panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0))+
                    ggtitle("A")

      #Image correlation plot
      coul <- viridis(100)
      w=dim(image_correl)[1]
      x <-rep(seq(1,w),w)
      y <-unlist(lapply(seq(1,w),function(x)rep(x,w)))
      z <-c(image_correl)
      df<-data.frame(x,y,z)
      plt_image_correl <- ggplot(df, aes(x, y, fill = z)) + geom_tile() +
                          xlab("Ion image") + ylab("Ion image") + theme_bw() +
                          scale_fill_gradientn("",limits = c(0,1), colours=coul) +
                          scale_x_continuous(breaks = seq(1, w),expand=c(0,0))+scale_y_continuous(breaks = seq(1, w),expand=c(0,0))+
                          theme(panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),panel.grid = element_blank(), plot.title = element_text(hjust = 0))+
                          ggtitle("D")

      #Raw spectrum plot
      raw_step<-(max(calculated_mass)-min(calculated_mass))/length(calculated_mass)
      raw_max<-max(calculated_mass)+raw_step
      raw_min<-min(calculated_mass)-raw_step
      raw_index<-which((full_spectrum$mass<raw_max)&(full_spectrum$mass>raw_min))
      raw_mass<-full_spectrum$mass[raw_index]
      raw_intensity<-full_spectrum$mean[raw_index]
      df1<-data.frame(x=raw_mass,y=raw_intensity)
      df2<-data.frame(calculated_mass=calculated_mass,calculated_magnitude=calculated_magnitude*(full_spectrum$mean[which.min(abs(full_spectrum$mass-experimental_mass[which.max(calculated_magnitude)]))])/max(calculated_magnitude)+min(raw_intensity),offset=min(raw_intensity))

      plt_raw <- ggplot(df1) + geom_line(aes(x, y,color="1"))+
        geom_vline(xintercept = calculated_mass,linetype="dashed",alpha=0.5)+
        geom_linerange(data=df2,aes(x=calculated_mass,ymin=offset,ymax=calculated_magnitude,color="2")) + geom_point(data=df2,aes(x=calculated_mass,y=calculated_magnitude,color="2"))+
        scale_color_discrete("Mean spectrum", c("Experimental","Calculated for Ag6"),breaks=c(1,2))+
        theme_bw()+xlab("m/z")+ylab("Intensity")+
        #scale_color_discrete("", c(paste("S1 (",round(prc1$auc.integral,2),"AUC )"), paste("S2 (",round(prc2$auc.integral,2),"AUC )"),paste("S1·S2 (",round(prc_product$auc.integral,2),"AUC )")))+
        scale_x_continuous(breaks = round(calculated_mass,4))+
        theme(legend.justification = c(1, 1), legend.position = c(1, 1),legend.background=element_rect(size=0.2,colour=1),panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0))+
        ggtitle("A")

      indices<-c(1,2)
      raw_step<-(max(calculated_mass)-min(calculated_mass))/length(calculated_mass)
      raw_max<-max(calculated_mass[indices])+0.2*raw_step
      raw_min<-min(calculated_mass[indices])-0.2*raw_step
      raw_index<-which((full_spectrum$mass<raw_max)&(full_spectrum$mass>raw_min))
      raw_mass<-full_spectrum$mass[raw_index]
      raw_intensity<-full_spectrum$mean[raw_index]
      df<-data.frame(x=raw_mass,y=raw_intensity)

      plt_zoom <- ggplot(df) + geom_line(aes(x, y,color="2"))+
        theme_bw()+xlab("m/z")+ylab("Intensity")+
        #scale_color_discrete("", c(paste("S1 (",round(prc1$auc.integral,2),"AUC )"), paste("S2 (",round(prc2$auc.integral,2),"AUC )"),paste("S1·S2 (",round(prc_product$auc.integral,2),"AUC )")))+
        scale_x_continuous(breaks = round(calculated_mass[indices],4))+
        theme(legend.position = "none",panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0))+
        ggtitle("C")


      #Recursive exploration plot
      text=""
      text=add_entry(text,"EXPLORATION MAP")
      text=paste(text,paste(apply(classification_record,1,paste,collapse=" "),collapse="\n"),"\n")
      text=add_entry(text,"SIMILARITY SCORES")
      text=paste(text,paste(apply(round(classes_list,2),1,paste,collapse=" "),collapse="\n"),"\n")

      plt_recursive_exploration=plot_text(text,size = 3)


      #Create a list with all the plots

      plts_peaks=list()

      #Print cluster images
      if(generate_figures)
        tiff(paste(pks$names[1],".tiff",sep=""), width = 200, height = 100, units = 'mm', res = 300)

      page_layout=default_page_layout
      offset=match(4,c(t(page_layout)))-1
      for(i in 1:num_peaks)
      {
        title=paste("m/z",round(calculated_mass[i],pkg_opt("round_digits")))
        # if(is.element(i,index_pks))
        #   title=paste(title,round(rel_error[i]*1e6),"ppm")
        # else
        #   title=paste(title,"(NOT IN PEAK MATRIX)")
        plts_peaks=c(plts_peaks,list(ggplot_peak_image(pks,final_image[,i],title,is.na(experimental_magnitude[i]),chosen[i],is.element(i,index_pks))))
        if((i+offset)%%length(page_layout)==0||i==num_peaks)
        {
          plts=c(list(plt_spectrum,plt_image_correl,plt_recursive_exploration),plts_peaks)
          grid.arrange(grobs=plts,layout_matrix=page_layout,top=global_title)
          plts=list()
          offset=0
          page_layout=matrix(1:length(page_layout), ncol = 4)
        }
      }

      if(generate_figures)
        dev.off()
    }
  }



  #Generate putative clusters
  #[Functionality demanded by Maria]

  #8.Summary of the report
  gt=pks$mass[which(results$gt)]
  gt_cluster_names=results$cluster_names[which(results$gt)]
  not_gt=setdiff(pks$mass,gt)
  if(generate_pdf)
  {
    if(include_summary)
    {
      #Print images of the ground truth
      if(length(gt)>0)
      {
        plts=list()
        gt_index=match(gt,pks$mass)
        for(i in 1:length(gt))
        {
          plts=c(plts,list(ggplot_peak_image(pks,pks[[mag_of_interest]][,gt_index[i]],paste("m/z",round(gt[i],pkg_opt("round_digits")),"(",gt_cluster_names[i],")"),chosen=T)))
          if(i%%length(page_layout)==0||i==length(gt))
          {
            grid.arrange(grobs=plts,nrow=nrow(page_layout),ncol=ncol(page_layout))
            plts=list()
          }
        }
      }
      #Print images of the not ground truth
      if(length(not_gt)>0)
      {
        plts=list()
        not_gt_index=match(not_gt,pks$mass)
        for(i in 1:length(not_gt))
        {
          plts=c(plts,list(ggplot_peak_image(pks,pks[[mag_of_interest]][,not_gt_index[i]],paste("m/z",round(not_gt[i],pkg_opt("round_digits"))),chosen=F)))
          if(i%%length(page_layout)==0||i==length(not_gt))
          {
            grid.arrange(grobs=plts,nrow=nrow(page_layout),ncol=ncol(page_layout))
            plts=list()
          }
        }
      }
    }

    #Print calculated clusters
    clusters=unique(results$patterns_out$cluster)
    cluster_indices=match(clusters,results$patterns_out$cluster)
    s1_scores=results$patterns_out$s1_scores[cluster_indices]
    s2_scores=results$patterns_out$s2_scores[cluster_indices]
    s3_scores=results$patterns_out$s3_scores[cluster_indices]

    #S1 vs. S2
    plot(s1_scores,s2_scores)
    text(s1_scores, s2_scores, labels=clusters, cex= 0.7, pos=3)
    abline(v = s1_threshold)
    abline(h = s2_threshold)

    plot(s1_scores,s2_scores,xlim = c(s1_threshold,1),ylim=c(s2_threshold,1))
    text(s1_scores, s2_scores, labels=clusters, cex= 0.7, pos=3)
    abline(v = s1_threshold)
    abline(h = s2_threshold)

    #S1 vs. S3
    plot(s1_scores,s3_scores)
    text(s1_scores, s3_scores, labels=clusters, cex= 0.7, pos=3)
    abline(v = s1_threshold)
    abline(h = s3_threshold)

    plot(s1_scores,s3_scores,xlim = c(s1_threshold,1),ylim=c(s3_threshold,1))
    text(s1_scores, s3_scores, labels=clusters, cex= 0.7, pos=3)
    abline(v = s1_threshold)
    abline(h = s3_threshold)

    #Close pdf file
    dev.off()
  }

  #Temporal
  # plts=c(list(plt_raw,grid.arrange(grobs=plts_peaks,nrow=1,top = grid::textGrob("\t B",x=0,hjust=0,gp=grid::gpar(fontsize=13,font="TT Arial"))),plt_zoom,plt_image_correl))
  plts=c(list(plt_raw,grid.arrange(grobs=plts_peaks,nrow=1,top = grid::textGrob("\t B",x=0,hjust=0,gp=grid::gpar(fontsize=13,font="TT Arial"))),plt_image_correl))

  tiff("/home/gbaquer/msidata/Ag Software Test 1/output/RESULTS/Fig3.tiff", width = 200, height = 300, units = 'mm', res = 300)
  #print(grid.arrange(grobs=plts,layout_matrix=matrix(c(1,2,3,1,2,4),nrow=3)))
  print(grid.arrange(grobs=plts,layout_matrix=matrix(c(1,2,3,1,2,3),nrow=3)))
  dev.off()
  return(results)
}



#' Compute scores
#'
#' Returns several scores for performance assessment of a given binary classification result.
#'
#' @param gt Ground truth: Positive values
#' @param pos Classified positives
#' @param neg Classified negatives
#'
#' @return List with all the computed scores. It includes F1 score and Brier score. As well as the number of
#'
#'
#' @export
compute_scores <- function (gt,pos,neg) {
  #compute tp, tn, fp, fn
  tp_list=intersect(gt,pos)
  fp_list=setdiff(pos,tp_list)
  fn_list=intersect(gt,neg)
  tn_list=setdiff(neg,fn_list)

  # tp_list=intersect(pos,gt) #true positive
  # fp_list=which(!is.element(pos,tp)) #setdiff(pos,tp) #false positive
  # fn_list=intersect(neg,gt) #false negative
  # tn_list=which(!is.element(neg,fn))#setdiff(neg,fn) #true negative

  tp=length(tp_list)
  fp=length(fp_list)
  fn=length(fn_list)
  tn=length(tn_list)

  p=tp/(tp+fp) #precision
  r=tp/(tp+fn) #recall

  scores=list()
  scores$f1= 2*(r*p)/(r+p) #F1 score
  scores$b= (fp+fn)/(length(pos)+length(neg)) #Brier score
  scores$tp=tp
  scores$fp=fp
  scores$fn=fn
  scores$tn=tn
  scores$p=p
  scores$r=r

  if(pkg_opt()$verbose_level<=-1)
  {
    print(paste("tp:",length(tp),"tn:",length(tn),"fp:",length(fp),"fn:",length(fn)))
    print(paste("p:",p,"r:",r))
    print(paste("F1 score:",scores$f1,"BS:",scores$b))
  }
  return(scores)
}


#' Compute S1
#'
#' Returns the spectral similarity score (S1)
#'
#' @param calculated_magnitude Calculated magnitude vector
#' @param experimental_magnitude Experimental magnitude vector
#' @param similarity_method Similarity method to be used
#'
#' @return Spectral similarity score S1
#'
#'
#' @export
compute_s1 <- function (calculated_magnitude,experimental_magnitude,similarity_method) {
  calculated_magnitude=calculated_magnitude/max(calculated_magnitude)
  experimental_magnitude=experimental_magnitude/max(experimental_magnitude)
  #Compute S1 using all peaks
  s1=exponential_decay_similarity(calculated_magnitude,experimental_magnitude,method=similarity_method)
  # if(num_peaks==1)
  #   s1=1
  # #Compute S1 using only the peak matrix
  # if(length(index_pks)>1)
  # {
  #   s1_pks=exponential_decay_similarity(calculated_magnitude[index_pks],experimental_magnitude[index_pks],method=similarity_method)
  # }
  # else
  #   s1_pks=NA

  #Compute S1
  # s1=max(s1_all,s1_pks,na.rm=T)

  #Check if there is an overlap
  # ion_clusters=kmeans(image_correl,centers=min(2,nrow(image_correl)-1))$clus
  if(is.na(s1))
    s1=0
  return(s1)
}

#' Compute S2
#'
#' Returns the intra-cluster morphological similarity score (S2)
#'
#' @param calculated_magnitude Calculated magnitude vector
#' @param image_correl Image correlation matrix of the experimental magnitude vector
#'
#' @return Spectral similarity score S2
#'
#'
#' @export
compute_s2 <- function (calculated_magnitude,image_correl) {
  #Compute S2 using all peaks
  s2_individual = apply(image_correl,2,weighted.mean,w=calculated_magnitude)
  s2_individual[which(is.na(s2_individual))]=0
  s2=weighted.mean(s2_individual,calculated_magnitude)

  # #Compute S2 using only the peak matrix
  # if(length(index_pks)>1)
  # {
  #   s2_individual_pks = apply(image_correl[index_pks,index_pks],2,weighted.mean,w=calculated_magnitude[index_pks])
  #   s2_individual_pks[which(is.na(s2_individual_pks))]=0
  #   s2_pks=weighted.mean(s2_individual_pks,calculated_magnitude[index_pks])
  # }
  # else
  #   s2_pks=NA
  #
  # #Compute S2
  # s2=max(s2_all,s2_pks,na.rm=T)
  if(is.na(s2))
    s2=0
  return(s2)
}

#' Compute S3
#'
#' Returns the mass tolerance similarity score (S3)
#'
#' @param calculated_mass Calculated mass vector
#' @param experimental_mass Experimental mass vector
#' @param similarity_method Similarity method to be used
#'
#' @return Mass tolerance score S3
#'
#'
#' @export
compute_s3 <- function (calculated_mass,experimental_mass,similarity_method) {
  #Compute S3
  norm_factor=max(experimental_mass)
  s3=exponential_decay_similarity(diff(calculated_mass)*1e4/norm_factor,diff(experimental_mass)*1e4/norm_factor,method=similarity_method,normalize = F)
  if(length(calculated_mass)==1)
  {
    s3=1
  }
  #s3=exponential_decay_similarity(calculated_mass[index_pks]/calculated_mass[index_pks],experimental_mass[index_pks]/calculated_mass[index_pks],method=similarity_method)

  #Check if there is an overlap
  # ion_clusters=kmeans(image_correl,centers=min(2,nrow(image_correl)-1))$clus
  if(is.na(s3))
    s3=0

  return(s3)

}
