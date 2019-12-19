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
viridis<-c("#440154FF","#450558FF","#46085CFF","#470D60FF","#471063FF","#481467FF","#481769FF","#481B6DFF","#481E70FF","#482173FF","#482576FF","#482878FF","#472C7AFF","#472F7CFF","#46327EFF","#453581FF","#453882FF","#443B84FF","#433E85FF","#424186FF","#404587FF","#3F4788FF","#3E4A89FF","#3D4D8AFF","#3C508BFF","#3B528BFF","#39558CFF","#38598CFF","#375B8DFF","#355E8DFF","#34608DFF","#33638DFF","#32658EFF","#31688EFF","#2F6B8EFF","#2E6D8EFF","#2D708EFF","#2C718EFF","#2B748EFF","#2A768EFF","#29798EFF","#287C8EFF","#277E8EFF","#26818EFF","#26828EFF","#25858EFF","#24878EFF","#238A8DFF","#228D8DFF","#218F8DFF","#20928CFF","#20938CFF","#1F968BFF","#1F998AFF","#1E9B8AFF","#1F9E89FF","#1FA088FF","#1FA287FF","#20A486FF","#22A785FF","#24AA83FF","#25AC82FF","#28AE80FF","#2BB07FFF","#2EB37CFF","#31B67BFF","#35B779FF","#39BA76FF","#3DBC74FF","#41BE71FF","#47C06FFF","#4CC26CFF","#51C56AFF","#56C667FF","#5BC863FF","#61CA60FF","#67CC5CFF","#6DCD59FF","#73D056FF","#78D152FF","#7FD34EFF","#85D54AFF","#8CD646FF","#92D741FF","#99D83DFF","#A0DA39FF","#A7DB35FF","#ADDC30FF","#B4DE2CFF","#BBDE28FF","#C2DF23FF","#C9E020FF","#D0E11CFF","#D7E219FF","#DDE318FF","#E4E419FF","#EBE51AFF","#F1E51DFF","#F7E620FF","#FDE725FF")
#' Annotate matrix
#'
#' Annotates the MS signals present in 'pks' related to the laser desorption/ionization promoting material used.
#'
#' @param pks Peak Matrix Image produced using rMSIproc
#' @param matrix_formula String giving the chemical formula of the matrix in the enviPat notation#'
#' @param full_spectrum Full spectrum before performing peak picking. It is used to give a higher degree of confidence to the S1 and S2 computation.
#' @param s1_threshold Correlation between the theoretical and the real spectral pattern above which a given cluster is considered to be present
#' @param s2_threshold Correlation between the spatial images of the peak in a cluster above which the peaks are included in the gold truth (gt)
#' @param s3_threshold Correlation between the spatial images of the peak in a cluster above which the peaks are included in the gold truth (gt)
#' @param analyzer_resolution Analyzer resolution that is used to merge nearby peaks in-silico as the equipment would in real life
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
annotate_matrix <- function (pks,matrix_formula,full_spectrum=NULL, normalization="None",
                              add_list=NULL, sub_list=NULL, max_charge=1,
                              analyzer_resolution=20000, tol_mode="scans",tol_ppm=200e-6,tol_scans=4,
                              s_threshold=0.65,s1_threshold=0.80,s2_threshold=0.80,
                              similarity_method="euclidean",correlation_method="pearson",
                              mag_of_interest="intensity",
                              min_multi=1,max_multi=10,  detect_overlapping=T) {

  #SECTION 0 :: Preprocessing
  metadata=list(pks_name=pks$names,full_name=full$name,matrix_formula=matrix_formula, normalization=normalization,
                num_pixels=pks$numPixels,num_peaks=length(pks$mass),min_mass=min(pks$mass),max_mass=max(pks$mass),
                add_list=add_list,sub_list=sub_list,max_charge=max_charge,
                analyzer_resolution=analyzer_resolution, tol_mode=tol_mode,tol_ppm=tol_ppm,tol_scans=tol_scans,
                s_threshold=s_threshold,s1_threshold=s1_threshold,s2_threshold=s2_threshold,similarity_method=similarity_method,correlation_method=correlation_method,
                mag_of_interest=mag_of_interest,
                min_multi=min_multi,max_multi=max_multi, detect_overlapping=detect_overlapping)
  # Normalize
  if(normalization!="None")
    pks[[mag_of_interest]]=pks[[mag_of_interest]]/pks$normalizations[[normalization]]

  #0. Load data
  isotopes=NULL
  data("isotopes", package = "enviPat", envir = environment())
  adducts=NULL
  data("adducts", package = "enviPat", envir= environment())

  #1.Generate list of possible chemical formulas [PENDING: Find a better way of computing all formulas. Wrap it in a function]
  if(!is.null(add_list))
    matrix_formula=append(matrix_formula,unlist(lapply(add_list,function(x) enviPat::mergeform(matrix_formula,x))))
  # if(!is.null(sub_list))
  #   base_forms=append(base_forms,unlist(lapply(sub_list,function(x) enviPat::subform(base_forms,x)[which(enviPat::check_ded(base_forms,x)=="FALSE")])))
  base_forms=NULL
  for(i in min_multi:max_multi)
  {
    base_forms=append(base_forms,multiform(matrix_formula,i))
  }
  sorted_ix=sort(check_chemform(isotopes,base_forms)$monoisotopic_mass,index.return=T)$ix
  base_forms=base_forms[sorted_ix]

  #2.Generate pattern list with enviPat
  cat("Computing Theoretical Patterns\n")

  theoretical_patterns=NULL
  multiplier=1
  while(TRUE)
  {
    #Multiply
    forms=multiform(base_forms,multiplier)

    #Multiple charge [PENDING: Charge implementation THE CURRENT IMPLEMENTATION IS MESSING UP WITH THE ADDUCT CLUSTERS]
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
    patterns=quiet(isowrap(isotopes,checked,resmass = FALSE,resolution = analyzer_resolution))
    # patterns=quiet(isowrap(isotopes,checked,charge=charge,resmass = FALSE,resolution = analyzer_resolution))

    #Append to final list
    for(i in 1:length(patterns))
    {
      if(!is.element(attributes(patterns)$names[i],theoretical_patterns$cluster))
      {
        theoretical_patterns$mass=append(theoretical_patterns$mass,patterns[[i]][,1]) #store masses
        theoretical_patterns[[mag_of_interest]]=append(theoretical_patterns[[mag_of_interest]],patterns[[i]][,2]) #store intensities
        theoretical_patterns$name=append(theoretical_patterns$name,rep(attributes(patterns)$names[i],length(patterns[[i]][,1]))) #store cluster names
        # theoretical_patterns$cluster=append(theoretical_patterns$cluster,rep(paste(attributes(patterns)$names[i],"_",charge[i],sep=""),length(patterns[[i]][,1]))) #store cluster names
      }
    }
    break;
    multiplier=multiplier+1
    if(multiplier>max_multi)
      break;
  }

  #6. Determine s1 and s2 scores for each calculated cluster
  mean_image=apply(pks[[mag_of_interest]],2,mean)
  sd_image=apply(pks[[mag_of_interest]],2,sd)/mean_image
  clusters=unique(theoretical_patterns$name)

  #Initialize output results list
  num_pattern_elements=length(theoretical_patterns$mass)

  theoretical_patterns$pks_index=rep(NA,num_pattern_elements)
  theoretical_patterns$in_pks=rep(T,num_pattern_elements)

  num_peaks_in_matrix = length(pks$mass)

  count=rep(0,num_peaks_in_matrix)

  results=list(
    s1=rep(NA,num_peaks_in_matrix),
    s2=rep(NA,num_peaks_in_matrix),
    s=rep(NA,num_peaks_in_matrix),
    cluster_name=rep(NA,num_peaks_in_matrix),
    type=factor(rep("not_matrix",num_peaks_in_matrix),levels=c("matrix","not_matrix","overlapped_matrix")),
    theoretical_patterns=theoretical_patterns,
    metadata=metadata)

  cat("Computing Similarity Scores\n")
  pb <- ProgressBar(max=100, stepLength=100/length(clusters))
  reset(pb)
  for(cluster in clusters)
  {
    # Calculated cluster
    calculated_index=which(theoretical_patterns$name==cluster)
    calculated_mass=theoretical_patterns$mass[calculated_index]
    calculated_magnitude=theoretical_patterns[[mag_of_interest]][calculated_index]


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
            cols=(col-tol_scans):(col+tol_scans)
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
    #s3=compute_s3(calculated_mass[index_pks],experimental_mass[index_pks],similarity_method)

    #Compute similarity scores for each possible merged clusters
    index_selected=1:num_peaks
    index_overlapped=NULL
    s1_clusters=NULL
    s2_clusters=NULL
    #s3_clusters=NULL

    a=kmeans(image_correl,centers=min(2,nrow(image_correl)-1))
    b=kmeans(rel_error,centers=min(2,nrow(image_correl)-1,length(unique(rel_error))))
    ion_clusters=paste(a$clus,b$clus)

    classification_record=rbind(rep(1,num_peaks))
    classes_list=rbind(c(1,s1,s2))
    # OVERLAPPING DETECTION
    if(detect_overlapping&(s1>s1_threshold)&(s2>s2_threshold) )
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
        for(child_class_priority in unique(sort(child_classes_priority,decreasing = T)))
        {
          index_class=which(child_classes_priority==child_class_priority)
          index_not_class=which(child_classes_priority!=child_class_priority)
          child_class=child_classes[index_class][1]
          #Skip if cluster only contains 1 peak
          #Skip if other cluster peaks are not higher than 50% of the theoretical
          if((length(index_class)<=1)|length(index_class)/num_peaks<0.5)
          {
            next
          }
          #Compute scores
          s1_tmp=compute_s1(calculated_magnitude[index_class],experimental_magnitude[index_class],similarity_method)
          s2_tmp=compute_s2(calculated_magnitude[index_class],image_correl[index_class,index_class])
          #s3_tmp=compute_s3(calculated_mass[intersect(index_pks,index_class)],experimental_mass[intersect(index_pks,index_class)],similarity_method)

          #Strore results
          classes_list=rbind(classes_list,c(child_class,s1_tmp,s2_tmp))
          #Finish if similarity conditions met
          if((s1_tmp>s1_threshold)&(s2_tmp>s2_threshold))
          {
            s1=s1_tmp
            s2=s2_tmp
            #s3=s3_tmp
            index_selected=index_class
            index_overlapped=(1:num_peaks)[-index_class]
            found=T
            break
          }
        }
      }
    }

    #Choose which peaks belong in the ground truth (gt)
    chosen=rep(F,num_peaks)
    chosen[index_selected]=(s1*s2)>s_threshold
    type=factor(rep("not_matrix",num_peaks),levels=c("not_matrix","matrix","overlapped_matrix"))
    if((s1*s2)>s_threshold){
      type[index_selected]="matrix"
      type[index_overlapped]="overlapped_matrix"
    }

    #Adjust magnitude in the NA mode for plotting
    if(is.null(full_spectrum))
      experimental_magnitude[index_full_spectrum]=NA

    #Store results

    results$s1[experimental_index]=s1
    results$s2[experimental_index]=s2
    results$s[experimental_index]=s1*s2
    results$cluster_name[experimental_index]=cluster
    results$type[experimental_index]=type

    results$theoretical_patterns$in_pks[calculated_index[index_full_spectrum]]=F
    results$theoretical_patterns$pks_index[calculated_index]=experimental_index

    count[experimental_index[which(chosen)]]=count[experimental_index[which(chosen)]]+1
    increase(pb)
  }
  class(results)<-"rMSIcleanupAnnotation"

  for(i in which(count>1))
    cat(paste("WARNING: More than one cluster was found for peak m/z=",pks$mass[i]," (i=",i,")\n"))
  cat("DONE\n")
  return(results)
}

#' Annotate matrix
#'
#' Annotates the MS signals present in 'pks' related to the laser desorption/ionization promoting material used.
#'
#' @param pks Peak Matrix Image produced using rMSIproc
#' @param results Annotation results
#'
#' @return Ground Truth: List of masses available in the image that correspond to the matrix.
#'
#' @export
remove_matrix <- function (pks,results) {
  pks_new=get_columns_peakMatrix(pks,which(!(results$type=="matrix")))
  return(pks_new)
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
  if(is.na(s2))
    s2=0
  return(s2)
}


