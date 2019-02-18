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
#' @param pks Peak Matrix Image
#' @param matching_method "loose": All peaks within a given tolerance are considered matrix peaks; "strict": Only peaks within a given tolerance that mantain the abundance ratios for each cluster are considered matrix peaks
#' @param s1_threshold Correlation between the theoretical and the real spectral pattern above which a given cluster is considered to be present and thus included in the gold truth (gt)
#' @param generate_pdf Boolean indicating whether to generate a pdf or not
#' @return Ground Truth: List of masses available in the image that correspond to the matrix.
#'
#' @export
generate_gt <- function (matrix_formula,pks,full_spectrum,
                         s1_threshold=0.85,s2_threshold=0.85,
                         MALDI_resolution=20000, tol_mode="ppm",tol_ppm=200e-6,tol_scans=4,
                         mag_of_interest="intensity",normalization="None",
                         max_multi=10, add_list=NULL, sub_list=NULL,
                         generate_pdf=F,default_page_layout=NULL) {
  #SECTION -1 :: Input validation
  #Adjust tolerance to ppm if no full_spectrum is provided
  if(is.na(full_spectrum))
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
    text=add_entry(text,as.character(Sys.time()))
    text=add_entry(text,"#")
    text=add_entry(text,"IMAGE INFORMATION")
    text=add_entry(text,"- File_name:",pks$names[1])
    text=add_entry(text,"- Number of peaks:",length(pks$mass))
    text=add_entry(text,"- Number of pixels:",pks$numPixels[1])
    text=add_entry(text,"- Mass Range: [",min(pks$mass),", ", max(pks$mass),"]")
    text=add_entry(text,"#")
    text=add_entry(text,"MATRIX INFORMATION")
    text=add_entry(text,"- Matrix formula:",matrix_formula)
    text=add_entry(text,"- Add list:",add_list)
    text=add_entry(text,"- Substract list:",sub_list)
    text=add_entry(text,"- Base forms:",base_forms)
    text=add_entry(text,"- S1 threshold:",s1_threshold)
    text=add_entry(text,"- S2 threshold:",s2_threshold)
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
  #5. Determine matches with masses
  gt=NULL
  mean_image=apply(pks[[mag_of_interest]],2,mean)
  clusters=unique(patterns_out$cluster)
  s1_scores=NULL
  s2_scores=NULL
  for(cluster in clusters)
  {
    #Calculated cluster
    calculated_index=which(patterns_out$cluster==cluster)
    calculated_mass=patterns_out$mass[calculated_index]
    calculated_magnitude=patterns_out[[mag_of_interest]][calculated_index]

    # Experimental cluster
    experimental_index=get_closest_peak(calculated_mass,pks$mass)
    experimental_mass=pks$mass[experimental_index]
    experimental_magnitude=mean_image[experimental_index]

    #Determine mass relative error
    if(tol_mode=="ppm")
    {
      rel_error=abs(experimental_mass-calculated_mass)/calculated_mass
      index_pks=which(rel_error<=tol_ppm) #Peaks to be taken from the peak matrix
      index_full_spectrum=which(rel_error>tol_ppm) #Peaks to be taken from the full spectrum
    }
    else
    {
      rel_error=abs(experimental_mass-calculated_mass)/calculated_mass
      indices=1:length(experimental_mass)
      index_pks=which(unlist(lapply(indices,function(i)is_within_scan_tol(experimental_mass[i],calculated_mass[i],full_spectrum$mass,tol_scans))))#Peaks to be taken from the peak matrix
      index_full_spectrum=setdiff(indices,index_pks) #Peaks to be taken from the full spectrum
    }

    #Generate images for each peak in the cluster
    final_image=NULL
    for(i in 1:length(calculated_index))
    {
      if(is.element(i,index_full_spectrum))
      {
        if(!is.na(full_spectrum))
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
    s1=cor(calculated_magnitude,experimental_magnitude)
    if(length(calculated_index)==1)
      s1=1

    s1_scores=append(s1_scores,s1)

    #Compute S2
    image_correl=cor(final_image)
    image_correl[which(is.na(image_correl))]=0
    s2_individual = apply(image_correl,2,weighted.mean,w=calculated_magnitude)
    s2_individual[which(is.na(s2_individual))]=0
    s2=weighted.mean(s2_individual,calculated_magnitude)

    s2_scores=append(s2_scores,s2)

    #Compute S2 using only the peak matrix
    if(length(index_pks)>1)
    {
      s2_individual_pks = apply(image_correl[index_pks,index_pks],2,weighted.mean,w=calculated_magnitude[index_pks])
      s2_individual_pks[which(is.na(s2_individual_pks))]=0
      s2_pks=weighted.mean(s2_individual_pks,calculated_magnitude[index_pks])
    }
    else
      s2_pks=NA

    #Choose which peaks belong in the ground truth (gt)
    chosen=rep(F,length(experimental_mass))
    if(!is.na(s1)&!(is.na(s2)&is.na(s2_pks)))
    {
      if(s1>s1_threshold & max(s2,s2_pks,na.rm=T)>s2_threshold)
      {
        chosen[index_pks]=T
        new_gt=experimental_mass[index_pks]
        gt=append(gt,new_gt)
      }
    }

    #Adjust magnitude in the NA mode
    if(is.na(full_spectrum))
      experimental_magnitude[index_full_spectrum]=NA

    #Generate plots
    if(pkg_opt()$verbose_level<=-1 || generate_pdf)
    {
      #Print progress to console
      print(paste(cluster,s1,s2,s2_pks))

      #Normalize
      if(sum(!is.na(experimental_magnitude))!=0 && max(experimental_magnitude,na.rm=T)!=0)
        norm_image_intensity=experimental_magnitude/max(experimental_magnitude,na.rm=T)
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
      linetype=rep("solid",length(calculated_mass)*2)
      linetype[index_full_spectrum+length(calculated_mass)]="dotted"
      #linetype=append(rep("solid",length(s2_individual)),rep("dotted",length(s2_individual)))
      plt_spectrum= ggplot(melted_df, aes(mass,value,color=variable,ymin=0,ymax=value) ) +
                    geom_linerange(linetype=linetype) + geom_point() + geom_text(aes(label=label,vjust=0)) +
                    ggtitle(paste(cluster,"; S1:",round(s1,2),"; S2:",round(s2,2),"; S2_p:",round(s2_pks,2),"; MAX:", round(max(experimental_magnitude,na.rm=T),2))) +
                    xlab("m/z") + ylab("Normalized Intensity") +
                    scale_colour_discrete(name="",breaks=c("calculated_magnitude","experimental_magnitude"),labels=c("Calculated m/z","Experimental m/z"))+
                    theme(legend.position="bottom")

      #Image correlation plot
      plt_image_correl=levelplot(image_correl,at=seq(-1,1,length.out = 100))
      plts=list(plt_spectrum,plt_image_correl)


      #Print cluster images
      page_layout=default_page_layout
      offset=match(3,c(t(page_layout)))-1
      for(i in 1:length(experimental_index))
      {
        title=paste("m/z",round(calculated_mass[i],pkg_opt("round_digits")))
        if(is.element(i,index_pks))
          title=paste(title," +/-",round(rel_error[i]*1e6),"ppm")
        else
          title=paste(title,"(NOT IN PEAK MATRIX)")
        plts=c(plts,list(ggplot_peak_image(pks,final_image[,i],title,is.na(experimental_magnitude[i]),chosen[i],is.element(i,index_pks))))
        if((i+offset)%%length(page_layout)==0||i==length(experimental_index))
        {
          grid.arrange(grobs=plts,layout_matrix=page_layout)
          plts=list()
          offset=0
          page_layout=matrix(1:length(page_layout), ncol = 4)
        }
      }
    }
  }

  gt=unique(gt)
  not_gt=setdiff(pks$mass,gt)

  #Generate putative clusters
  #[Functionality demanded by Maria]

  if(generate_pdf)
  {
    #Print images of the ground truth
    if(length(gt)>0)
    {
      plts=list()
      gt_index=match(gt,pks$mass)
      for(i in 1:length(gt))
      {
        plts=c(plts,list(ggplot_peak_image(pks,pks[[mag_of_interest]][,gt_index[i]],paste("m/z",round(gt[i],pkg_opt("round_digits"))),chosen=T)))
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

    #Close pdf file
    dev.off()
  }
  return(gt)
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
