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
#' @param default_page_layout Page layout to be used in the pdf plotting.
#'
#' @return Ground Truth: List of masses available in the image that correspond to the matrix.
#'
#' @export
generate_gt <- function (matrix_formula,pks,full_spectrum=NULL,
                         s1_threshold=0.80,s2_threshold=0.85, s3_threshold=0.7, similarity_method="euclidean",
                         MALDI_resolution=20000, tol_mode="ppm",tol_ppm=200e-6,tol_scans=4,
                         mag_of_interest="intensity",normalization="None",
                         max_multi=10, add_list=NULL, sub_list=NULL,
                         generate_pdf=F,default_page_layout=NULL) {
  #SECTION -1 :: Input validation
  #Adjust tolerance to ppm if no full_spectrum is provided
  if(is.null(full_spectrum))
    tol_mode="ppm"
  #Page layout
  if(is.null(default_page_layout))
  {
    default_page_layout=rbind(c(1,1,2,2),
                      c(1,1,2,2),
                      c(3,4,5,6),
                      c(7,8,9,10),
                      c(11,12,13,14),
                      c(15,16,17,18))
  }
  #SECTION 0 :: Preprocessing

  # Select first image if there are multiple
  pks=get_one_peakMatrix(pks)

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
  base_forms=NULL
  for(i in 1:max_multi)
  {
    base_forms=append(base_forms,multiform(matrix_formula,i))
  }
  if(!is.null(add_list))
    base_forms=append(base_forms,unlist(lapply(add_list,function(x) enviPat::mergeform(base_forms,x))))
  if(!is.null(sub_list))
    base_forms=append(base_forms,unlist(lapply(sub_list,function(x) enviPat::subform(base_forms,x)[which(enviPat::check_ded(base_forms,x)=="FALSE")])))

  sorted_ix=sort(check_chemform(isotopes,base_forms)$monoisotopic_mass,index.return=T)$ix
  base_forms=base_forms[sorted_ix]

  #4. Open pdf report
  if(generate_pdf)
  {
    #Generate pdf name [pks$name_001]
    pdf_file=generate_file_name(pks$names[1])
    #open pdf file
    a4_width=8.27
    a4_height=11.69
    pdf(pdf_file,paper='a4',width=a4_width,height=a4_height)

    #First page of metadata [File name, mean image, matrix formula, adducts list, base forms, S1_threshold]
    text=""
    text=add_entry(text,"#")
    text=add_entry(text,"- Package Version:",packageVersion("rMSIcleanup"))
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
    checked=check_chemform(isotopes,forms)

    #Remove incorrect formulas
    checked=checked[which(checked$warning==F),]

    #Lower mass limit
    checked=checked[which(checked$monoisotopic_mass>=min(pks$mass)),]
    if(nrow(checked)==0)
      next;

    #Upper mass limit
    checked=checked[which(checked$monoisotopic_mass<=max(pks$mass)),]
    if(nrow(checked)==0)
      break;

    #Compute theoretical patterns
    patterns=quiet(isowrap(isotopes,checked,resmass = FALSE,resolution = MALDI_resolution))

    #Append to final list
    for(i in 1:length(patterns))
    {
      if(!is.element(attributes(patterns)$names[i],patterns_out$cluster))
      {
        patterns_out$mass=append(patterns_out$mass,patterns[[i]][,1]) #store masses
        patterns_out[[mag_of_interest]]=append(patterns_out[[mag_of_interest]],patterns[[i]][,2]) #store intensities
        patterns_out$cluster=append(patterns_out$cluster,rep(attributes(patterns)$names[i],length(patterns[[i]][,1]))) #store cluster names
      }
    }
    multiplier=multiplier+1
  }

  #6. Determine s1 and s2 scores for each calculated cluster
  mean_image=apply(pks[[mag_of_interest]],2,mean)
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

    #Compute S1
    s1=exponential_decay_similarity(calculated_magnitude[index_pks],experimental_magnitude[index_pks],method=similarity_method)
    if(num_peaks==1)
      s1=1

    #Compute S2 using all peaks
    image_correl=cor(final_image)
    image_correl[which(is.na(image_correl))]=0
    s2_individual = apply(image_correl,2,weighted.mean,w=calculated_magnitude)
    s2_individual[which(is.na(s2_individual))]=0
    s2_all=weighted.mean(s2_individual,calculated_magnitude)

    #Compute S2 using only the peak matrix
    if(length(index_pks)>1)
    {
      s2_individual_pks = apply(image_correl[index_pks,index_pks],2,weighted.mean,w=calculated_magnitude[index_pks])
      s2_individual_pks[which(is.na(s2_individual_pks))]=0
      s2_pks=weighted.mean(s2_individual_pks,calculated_magnitude[index_pks])
    }
    else
      s2_pks=NA

    #Compute S2
    s2=max(s2_all,s2_pks,na.rm=T)

    #Compute S3
    a=calculated_mass[index_pks]
    b=experimental_mass[index_pks]
    s3=rMSIcleanup::exponential_decay_similarity(diff(a)*1e4/max(b),diff(b)*1e4/max(b),method=similarity_method,normalize = F)
    #s3=exponential_decay_similarity(calculated_mass[index_pks]/calculated_mass[index_pks],experimental_mass[index_pks]/calculated_mass[index_pks],method=similarity_method)

    #Choose which peaks belong in the ground truth (gt)
    chosen=rep(F,num_peaks)
    chosen[index_pks]=(!is.na(s1)&!is.na(s2)&!is.na(s3))&(s1>s1_threshold & s2>s2_threshold & s3>s3_threshold)


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
    if(pkg_opt()$verbose_level<=-1 || generate_pdf)
    {
      #Print progress to console
      print(paste(cluster,s1,s2,s3))

      #Normalize
      if(sum(!is.na(experimental_magnitude))!=0 && max(experimental_magnitude,na.rm=T)!=0)
        norm_image_intensity=experimental_magnitude/max(experimental_magnitude,na.rm=T)
      else
        norm_image_intensity=experimental_magnitude
      if(max(calculated_magnitude)!=0)
        calculated_magnitude=calculated_magnitude/max(calculated_magnitude)

      #Prepare melted dataframe
      df=data.frame(mass=calculated_mass,calculated_magnitude=calculated_magnitude,experimental_magnitude=norm_image_intensity)
      melted_df=melt(df, id.vars="mass",measure.vars=c("calculated_magnitude","experimental_magnitude"))
      value=NULL
      variable=NULL
      mass=NULL

      #Spectrum comparison plot
      label = append(rep("",length(s2_individual)),round(s2_individual,2))
      linetype=rep("solid",num_peaks*2)
      linetype[index_full_spectrum+num_peaks]="dotted"
      #linetype=append(rep("solid",length(s2_individual)),rep("dotted",length(s2_individual)))
      plt_spectrum= ggplot(melted_df, aes(mass,value,color=variable,ymin=0,ymax=value) ) +
                    geom_linerange(linetype=linetype) + geom_point() + geom_text(aes(label=label,vjust=0)) +
                    ggtitle(paste(cluster,"; S1:",round(s1,2),"; S2:",round(s2,2),"; S3:",round(s3,2),"; MAX:", round(max(experimental_magnitude,na.rm=T),2))) +
                    xlab("m/z") + ylab("Normalized Intensity") +
                    scale_colour_discrete(name="",breaks=c("calculated_magnitude","experimental_magnitude"),labels=c("Calculated m/z","Experimental m/z"))+
                    theme(legend.position="bottom")

      #Image correlation plot
      plt_image_correl=levelplot(image_correl,at=seq(-1,1,length.out = 100))
      plts=list(plt_spectrum,plt_image_correl)


      #Print cluster images
      page_layout=default_page_layout
      offset=match(3,c(t(page_layout)))-1
      for(i in 1:num_peaks)
      {
        title=paste("m/z",round(calculated_mass[i],pkg_opt("round_digits")))
        if(is.element(i,index_pks))
          title=paste(title,round(rel_error[i]*1e6),"ppm")
        else
          title=paste(title,"(NOT IN PEAK MATRIX)")
        plts=c(plts,list(ggplot_peak_image(pks,final_image[,i],title,is.na(experimental_magnitude[i]),chosen[i],is.element(i,index_pks))))
        if((i+offset)%%length(page_layout)==0||i==num_peaks)
        {
          grid.arrange(grobs=plts,layout_matrix=page_layout)
          plts=list()
          offset=0
          page_layout=matrix(1:length(page_layout), ncol = 4)
        }
      }
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
