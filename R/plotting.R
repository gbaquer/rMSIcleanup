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

## printing goes somewhere else
# annotate_matrix <- function (folder="output/",generate_pdf=T,generate_figures=F,default_page_layout=NULL,plot_type="debug",include_summary=F) {


#' Annotate matrix
#'
#' Annotates the MS signals present in 'pks' related to the laser desorption/ionization promoting material used.
#'
#' @param pks Peak Matrix Image produced using rMSIproc
#' @param matrix_formula String giving the chemical formula of the matrix in the enviPat notation#'
#' @param full Full spectrum before performing peak picking. It is used to give a higher degree of confidence to the S1 and S2 computation.
#' @param s1_threshold Correlation between the theoretical and the real spectral pattern above which a given cluster is considered to be present
#' @param s2_threshold Correlation between the spatial images of the peak in a cluster above which the peaks are included in the gold truth (gt)
#' @param s3_threshold Correlation between the spatial images of the peak in a cluster above which the peaks are included in the gold truth (gt)
#' @param MALDI_resolution MALDI resolution that is used to merge nearby peaks in-silico as the equipment would in real life
#' @param tol_mode String determining the tolerance mode to be used when comparing the calculated masses with the experimental ones. "ppm": relative tolerance with respect to the calculated one in parts per million. "scans": Number of scans or datapoints present in the full spectrum -It is only applicable if full is provided.
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
generate_pdf_report <- function (results, pks, full=NULL, filename=results$metadata$pks_name,folder="output/",default_page_layout=NULL,include_summary=F, include_PCA=F) {
  #0. Input Validation
  m <-results$metadata
  if((pks$names!=m$pks_name) || (full$name!=m$full_name))
    cat("WARNING: The names of the pks and full files provided do not match with the ones used during matrix annotation\n")
  if(is.null(default_page_layout))
  {
    default_page_layout=rbind(c(1,1,1,2,2),
                              c(1,1,1,2,2),
                              c(1,1,1,2,2),
                              c(3,4,5,6,7),
                              c(8,9,10,11,12),
                              c(13,14,15,16,17))

  }
  #1. Generate filename & open pdf file
  pdf_file=generate_file_name(filename,folder = folder)
  a4_width=8.27
  a4_height=11.69
  pdf(pdf_file,paper='a4r',width=a4_height,height=a4_width)

  #2. PRINT METADATA
  text=""
  text=add_entry(text,"#")
  text=add_entry(text,"- Package Version:",paste(packageVersion("rMSIcleanup"),collapse = "."))
  text=add_entry(text,"- Time:",as.character(Sys.time()))
  text=add_entry(text,"#")
  text=add_entry(text,"IMAGE INFORMATION")
  text=add_entry(text,"- Peak matrix:",m$pks_name)
  if(is.null(m$full_name))
    text=add_entry(text,"- Full spectrum: [NOT PROVIDED]")
  else
    text=add_entry(text,"- Full spectrum:",m$full_name)
  text=add_entry(text,"- Number of peaks:",m$num_peaks)
  text=add_entry(text,"- Number of pixels:",m$num_pixels)
  text=add_entry(text,"- Mass Range: [",m$min_mass,", ", m$max_mass,"]")
  text=add_entry(text,"#")
  text=add_entry(text,"MATRIX INFORMATION")
  text=add_entry(text,"- Matrix formula:",m$matrix_formula)
  text=add_entry(text,"- Add list:",m$add_list)
  text=add_entry(text,"- Substract list:",m$sub_list)
  text=add_entry(text,"- Maximum cluster multiplication:",m$max_multi)
  text=add_entry(text,"- Base forms:",unique(results$theoretical_patterns$name))
  text=add_entry(text,"#")
  text=add_entry(text,"PROCESSING INFORMATION")
  text=add_entry(text,"- S threshold:",m$s_threshold)
  text=add_entry(text,"- S1 threshold:",m$s1_threshold)
  text=add_entry(text,"- S2 threshold:",m$s2_threshold)
  text=add_entry(text,"- Similarity method:",m$similarity_method)
  text=add_entry(text,"- Magnitude of interest:",m$mag_of_interest)
  text=add_entry(text,"- Tolerance mode:",m$tol_mode)
  if(m$tol_mode=="ppm")
    text=add_entry(text,"- Tolerance ppm:",m$tol_ppm)
  else
    text=add_entry(text,"- Tolerance scans:",m$tol_scans)
  text=add_entry(text,"#")

  text=paste(text,collapse = "\n")
  print(plot_text(text))


  #3. PRINT REPORT FOR EACH CLUSTER
  #Load Results data
  theoretical_patterns=results$theoretical_patterns
  num_pattern_elements=length(theoretical_patterns$mass)
  clusters=unique(theoretical_patterns$name)

  mag_of_interest=m$mag_of_interest
  mean_image=apply(pks[[mag_of_interest]],2,mean)
  sd_image=apply(pks[[mag_of_interest]],2,sd)/mean_image

  num_peaks_in_matrix = length(pks$mass)

  #Loop through all clusters
  cat("Plotting Clusters\n")
  pb <- ProgressBar(max=100, stepLength=100/length(clusters))
  reset(pb)
  for(cluster in clusters)
  {
    # Calculated cluster
    calculated_index=which(theoretical_patterns$name==cluster)
    calculated_mass=theoretical_patterns$mass[calculated_index]
    calculated_magnitude=theoretical_patterns[[mag_of_interest]][calculated_index]


    # Experimental cluster
    experimental_index=theoretical_patterns$pks_index[calculated_index]
    experimental_mass=pks$mass[experimental_index]
    experimental_magnitude=mean_image[experimental_index]
    experimental_sd=sd_image[experimental_index]

    # Number of peaks in cluster
    num_peaks=length(calculated_index)

    #Generate images for each peak in the cluster
    final_image=NULL
    index_full_spectrum=which(results$theoretical_patterns$pks_index==0)

    for(i in 1:num_peaks)
    {
      if(!results$theoretical_patterns$in_pks[i])
      {
        if(!is.null(full))
        {
          #Get image from full (Processed spectrum before peak picking)
          mass=calculated_mass[i]
          if(m$tol_mode=="ppm")
            data=rMSI::loadImageSliceFromMass(full, mass,mass*tol_ppm)$data
          else
          {
            col=which.min(abs(full$mass-mass))
            cols=(col-m$tol_scans/10):(col+m$tol_scans/10)
            data=rMSI::loadImageSliceFromCols(full,cols)
          }

          if(m$normalization!="None")
            data=data/pks$normalizations[[m$normalization]]
          max_col=which.max(apply(data,2,mean))
          data_picture=data[,max_col]
        }
        else
        {
          #Set image to NA
          data_picture=rep(NA,pks$numPixels[1])
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

    image_correl=cor(final_image,method=m$correlation_method)
    image_correl[which(is.na(image_correl))]=0

    #Compute similarity scores for each possible merged clusters
    # index_selected=1:num_peaks
    # s1_clusters=NULL
    # s2_clusters=NULL
    # s3_clusters=NULL
    #
    # a=kmeans(image_correl,centers=min(2,nrow(image_correl)-1))
    # b=kmeans(rel_error,centers=min(2,nrow(image_correl)-1,length(unique(rel_error))))
    # ion_clusters=paste(a$clus,b$clus)
    #
    # classification_record=rbind(rep(1,num_peaks))
    # classes_list=rbind(c(1,s1,s2,s3))


    #Choose which peaks belong in the ground truth (gt)
    chosen=rep(F,num_peaks)
    chosen[which(results$type=="matrix")]=T

    #Adjust magnitude in the NA mode for plotting
    if(is.null(full))
      experimental_magnitude[index_full_spectrum]=NA

    #PLOTS
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
      plot_type="poster"
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

      s1=results$s1[experimental_index[1]]
      s2=results$s2[experimental_index[1]]
      global_title=paste(cluster,"; S:",round(s1*s2,2),"; S1:",round(s1,2),"; S2:",round(s2,2),"; MAX:", round(max(experimental_magnitude,na.rm=T),2))

      plt_spectrum= ggplot(melted_df, aes(mass,value,color=variable,ymin=0,ymax=value)) +
        geom_linerange() + geom_point() + theme_bw() +
        # ggtitle(paste(cluster,"; S1:",round(s1,2),"; S2:",round(s2,2),"; S3:",round(s3,2),"; MAX:", round(max(experimental_magnitude,na.rm=T),2))) +
        xlab("m/z") + ylab("Scaled Intensity") +
        scale_colour_discrete(name="",breaks=vars,labels=vars_label)+
        theme(legend.justification = c(1, 1), legend.position = c(1, 1),panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0))

      #Image correlation plot
      coul <- viridis
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
              panel.grid.minor = element_blank(),panel.grid = element_blank(), plot.title = element_text(hjust = 0))

      #Raw spectrum plot
      raw_step<-(max(calculated_mass)-min(calculated_mass))/length(calculated_mass)
      raw_max<-max(calculated_mass)+raw_step
      raw_min<-min(calculated_mass)-raw_step
      raw_index<-which((full$mass<raw_max)&(full$mass>raw_min))
      raw_mass<-full$mass[raw_index]
      raw_intensity<-full$mean[raw_index]
      df1<-data.frame(x=raw_mass,y=raw_intensity)
      df2<-data.frame(calculated_mass=calculated_mass,calculated_magnitude=calculated_magnitude*(full$mean[which.min(abs(full$mass-experimental_mass[which.max(calculated_magnitude)]))])/max(calculated_magnitude)+min(raw_intensity),offset=min(raw_intensity))

      plt_raw <- ggplot(df1) + geom_line(aes(x, y,color="1"))+
        geom_vline(xintercept = calculated_mass,linetype="dashed",alpha=0.5)+
        geom_linerange(data=df2,aes(x=calculated_mass,ymin=offset,ymax=calculated_magnitude,color="2")) + geom_point(data=df2,aes(x=calculated_mass,y=calculated_magnitude,color="2"))+
        scale_color_discrete("Mean spectrum", c("Experimental","Calculated for Ag6"),breaks=c(1,2))+
        theme_bw()+xlab("m/z")+ylab("Intensity")+
        #scale_color_discrete("", c(paste("S1 (",round(prc1$auc.integral,2),"AUC )"), paste("S2 (",round(prc2$auc.integral,2),"AUC )"),paste("S1·S2 (",round(prc_product$auc.integral,2),"AUC )")))+
        scale_x_continuous(breaks = round(calculated_mass,4))+
        theme(legend.justification = c(1, 1), legend.position = c(1, 1),legend.background=element_rect(size=0.2,colour=1),panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0))

      indices<-c(1,2)
      raw_step<-(max(calculated_mass)-min(calculated_mass))/length(calculated_mass)
      raw_max<-max(calculated_mass[indices])+0.2*raw_step
      raw_min<-min(calculated_mass[indices])-0.2*raw_step
      raw_index<-which((full$mass<raw_max)&(full$mass>raw_min))
      raw_mass<-full$mass[raw_index]
      raw_intensity<-full$mean[raw_index]
      df<-data.frame(x=raw_mass,y=raw_intensity)

      plt_zoom <- ggplot(df) + geom_line(aes(x, y,color="2"))+
        theme_bw()+xlab("m/z")+ylab("Intensity")+
        #scale_color_discrete("", c(paste("S1 (",round(prc1$auc.integral,2),"AUC )"), paste("S2 (",round(prc2$auc.integral,2),"AUC )"),paste("S1·S2 (",round(prc_product$auc.integral,2),"AUC )")))+
        scale_x_continuous(breaks = round(calculated_mass[indices],4))+
        theme(legend.position = "none",panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0))


      #Recursive exploration plot
      text=""
      text=add_entry(text,"EXPLORATION MAP")
      # text=paste(text,paste(apply(classification_record,1,paste,collapse=" "),collapse="\n"),"\n")
      # text=add_entry(text,"SIMILARITY SCORES")
      # text=paste(text,paste(apply(round(classes_list,2),1,paste,collapse=" "),collapse="\n"),"\n")

      plt_recursive_exploration=plot_text(text,size = 3)


      #Create a list with all the plots

      plts=c(list(plt_spectrum,plt_image_correl))


      page_layout=default_page_layout
      offset=match(4,c(t(page_layout)))-1
      type=results$type[experimental_index]
      in_pks=results$theoretical_patterns$in_pks[calculated_index]
      for(i in 1:num_peaks)
      {
        title=paste("m/z",round(calculated_mass[i],4))
        # if(is.element(i,index_pks))
        #   title=paste(title,round(rel_error[i]*1e6),"ppm")
        # else
        #   title=paste(title,"(NOT IN PEAK MATRIX)")
        plts=c(plts,list(ggplot_peak_image(pks,final_image[,i],title,type[i],in_pks[i])))
        if((i+offset)%%length(page_layout)==0||i==num_peaks)
        {
          grid.arrange(grobs=plts,layout_matrix=page_layout,top=global_title)
          plts=list()
          offset=0
          page_layout=matrix(1:length(page_layout), ncol = 4)
        }
      }
      increase(pb)
    }

  # 4.PRINT GLOBAL PLOTS
  cat("Plotting Summary\n")
  #[Create own function]
  cluster_indices=match(clusters,results$cluster_name)
  s1_scores=results$s1[cluster_indices]
  s2_scores=results$s2[cluster_indices]

  #S1 vs. S2
  colors=rep(0, length(s1_scores))
  df <- data.frame(x=s1_scores,y=s2_scores,col=as.factor(colors))
  p_a<-ggplot(df,aes(x,y,color=col))+ geom_point() +theme_bw()+xlab("S1")+ylab("S2")+
    #scale_color_discrete("Ground truth", c("Negative class", "Positive class"),breaks=c(0,1))+
    scale_x_continuous(breaks = seq(0, 1, by = 0.2))+scale_y_continuous(breaks = seq(0, 1, by = 0.2))#+
    # theme(legend.justification = c(1, 0), legend.position = c(1, 0),legend.background=element_rect(size=0.2,colour=1),panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(),
          # panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0))

  # S vs cluster number
  scores=s1_scores*s2_scores
  mean_scores=scores
  # for(clus in clusters)
  # {
  #   mean_scores=append(mean_scores,mean(scores[which(clusters==clus)]))
  # }
  sorted_mean_scores=sort(mean_scores,decreasing = T,index.return=T)


  unique_ids=1:length(clusters)
  names(unique_ids)=clusters[sorted_mean_scores$ix]
  lab=rep("",length(clusters))

  df <- data.frame(x=unique_ids[clusters],y=scores,z=as.factor(colors), lab=rep("  \n   ",length(clusters)))
  p_b<-ggplot(df,aes(label=lab))  + geom_point(aes(x, y, color=z))+
    theme_bw()+xlab("Cluster number")+ylab("S1·S2")+
    # scale_color_discrete("Ground truth", c("Negative class", "Positive class"),breaks=c(0,1))+
    scale_x_continuous(breaks = seq(0, 85, by = 5))+scale_y_continuous(breaks = seq(0, 1, by = 0.2))#+
    # theme(legend.justification = c(1, 1), legend.position = c(1, 1),legend.background=element_rect(size=0.2,colour=1),panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(),
    #       panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0))


  plts=c(list(p_a,p_b))
  #Before and after PCA
  if(include_PCA){
    pca_model_before<-prcomp(pks$intensity)
    pca_before<-as.data.frame(pca_model_before$x)

    pks_clean<-rMSIcleanup::remove_matrix(pks,results)
    pca_model_after<-prcomp(pks_clean$intensity)
    pca_after<-as.data.frame(pca_model_after$x)

    rgb_channels=list(c(1,0,0),
                      c(0,1,0),
                      c(0,0,1),
                      c(1,1,1))

    plots_pca_before=lapply(rgb_channels,ggplot_rgb,pos=pks$pos,rgb_data=pca_before)
    plots_pca_after=lapply(rgb_channels,ggplot_rgb,pos=pks_clean$pos,rgb_data=pca_after)
    plots_pca = arrangeGrob(grobs=c(plots_pca_after,plots_pca_after),nrow=2)
    plts=c(plts,list(plots_pca))
  }

  quiet(print(grid.arrange(grobs=plts,layout_matrix=matrix(c(1,3,2,3),nrow=2))))

  #Close pdf file
  dev.off()
  cat("DONE\n")
}

generate_figures_2 <- function (folder="RESULTS") {
  #1. LOAD RESULTS
  results_path = paste("/home/gbaquer/msidata/Ag Software Test 1/output/",folder,"/global_results_000.rds",sep="")
  results = readRDS(results_path)

  output_dir = paste("/home/gbaquer/msidata/Ag Software Test 1/output/",folder,"/",sep="")

  #2. Prepare data

  dataset_list=NULL
  clusters_list=NULL
  classification_list=NULL

  labels=NULL
  colors=NULL

  s1_scores=NULL
  s2_scores=NULL
  s3_scores=NULL

  tp_scores=NULL
  fp_scores=NULL
  tn_scores=NULL
  fn_scores=NULL

  fp_rate_scores=NULL
  tp_rate_scores=NULL

  i=1
  truth=c("Ag1","Ag2","Ag3","Ag4","Ag5","Ag6","Ag7","Ag8","Ag9","Ag10")
  global_clusters=unique(unlist(lapply(results$data,function(x)unique(x$patterns_out$cluster))))
  dataset_names=NULL
  for(r in results$data)
  {
    #Print calculated clusters
    present=r$patterns_out$present
    clusters=r$patterns_out$cluster[present]
    unique_clusters=unique(clusters)

    cluster_indices=match(unique_clusters,clusters)

    s1=r$patterns_out$s1_scores[present]
    s2=r$patterns_out$s2_scores[present]
    s3=r$patterns_out$s3_scores[present]

    s1=s1[cluster_indices]
    s2=s2[cluster_indices]
    s3=s3[cluster_indices]

    #labels=append(labels,paste(i,"_",clusters,sep="")[which(present)])
    dataset_list=append(dataset_list,rep(i,length(s1)))

    is_na=which(is.na(s1)|is.na(s2)|is.na(s3))
    col=is.element(unique_clusters,truth)
    col[is_na]=1
    colors=append(colors,col)

    s1[is.na(s1)]=0
    s2[is.na(s2)]=0
    s3[is.na(s3)]=0

    s1_scores=append(s1_scores,s1)
    s2_scores=append(s2_scores,s2)
    s3_scores=append(s3_scores,s3)

    clusters_list=append(clusters_list,unique_clusters)

    # #Compute positives and negatives
    # chosen=(s1>results$meta$s1_threshold&s2>results$meta$s2_threshold&s3>results$meta$s3_threshold)
    # should_be_chosen=is.element(clusters,truth)
    #
    # tp=sum(chosen&should_be_chosen&s3!=0)
    # fp=sum(chosen&!should_be_chosen&s3!=0)
    # tn=sum(!chosen&!should_be_chosen&s3!=0)
    # fn=sum(!chosen&should_be_chosen&s3!=0)
    #
    # classification=rep("np",length(global_clusters))
    # global_index=match(clusters,global_clusters)
    # classification[global_index[which(chosen&should_be_chosen&s3!=0)]]="tp"
    # classification[global_index[which(chosen&!should_be_chosen&s3!=0)]]="fp"
    # classification[global_index[which(!chosen&!should_be_chosen&s3!=0)]]="tn"
    # classification[global_index[which(!chosen&should_be_chosen&s3!=0)]]="fn"
    # classification_list=rbind(classification_list,classification)
    #
    # fp_rate=fp/(fp+tn)
    # tp_rate=tp/(tp+fn)
    #
    # tp_scores=append(tp_scores,tp)
    # fp_scores=append(fp_scores,fp)
    # tn_scores=append(tn_scores,tn)
    # fn_scores=append(fn_scores,fn)
    #
    # fp_rate_scores=append(fp_rate_scores,fp_rate)
    # tp_rate_scores=append(tp_rate_scores,tp_rate)

    i=i+1

  }
  # clusters=unique(results$data[[1]]$patterns_out$cluster)
  # s1_roc_fp_rate=NULL
  # s1_roc_tp_rate=NULL
  # s2_roc_fp_rate=NULL
  # s2_roc_tp_rate=NULL
  # s3_roc_fp_rate=NULL
  # s3_roc_tp_rate=NULL
  # all_roc_fp_rate=NULL
  # all_roc_tp_rate=NULL
  # for(t in seq(1,-0.001,length=1000))
  # {
  #   #S1 ROC
  #   chosen=s1_scores>t
  #   should_be_chosen=is.element(clusters_list,truth)
  #
  #   tp=sum(chosen&should_be_chosen)
  #   fp=sum(chosen&!should_be_chosen)
  #   tn=sum(!chosen&!should_be_chosen)
  #   fn=sum(!chosen&should_be_chosen)
  #
  #   s1_roc_fp_rate=append(s1_roc_fp_rate,fp/(fp+tn))
  #   s1_roc_tp_rate=append(s1_roc_tp_rate,tp/(tp+fn))
  #
  #   #S2 ROC
  #   chosen=s2_scores>t
  #   should_be_chosen=is.element(clusters_list,truth)
  #
  #   tp=sum(chosen&should_be_chosen)
  #   fp=sum(chosen&!should_be_chosen)
  #   tn=sum(!chosen&!should_be_chosen)
  #   fn=sum(!chosen&should_be_chosen)
  #
  #   s2_roc_fp_rate=append(s2_roc_fp_rate,fp/(fp+tn))
  #   s2_roc_tp_rate=append(s2_roc_tp_rate,tp/(tp+fn))
  #
  #   #S3 ROC
  #   chosen=s3_scores>t
  #   should_be_chosen=is.element(clusters_list,truth)
  #
  #   tp=sum(chosen&should_be_chosen)
  #   fp=sum(chosen&!should_be_chosen)
  #   tn=sum(!chosen&!should_be_chosen)
  #   fn=sum(!chosen&should_be_chosen)
  #
  #   s3_roc_fp_rate=append(s3_roc_fp_rate,fp/(fp+tn))
  #   s3_roc_tp_rate=append(s3_roc_tp_rate,tp/(tp+fn))
  #
  #   #ALL ROC
  #   chosen=s1_scores>t&s2_scores>t&s3_scores>t
  #   should_be_chosen=is.element(clusters_list,truth)
  #
  #   tp=sum(chosen&should_be_chosen)
  #   fp=sum(chosen&!should_be_chosen)
  #   tn=sum(!chosen&!should_be_chosen)
  #   fn=sum(!chosen&should_be_chosen)
  #
  #   all_roc_fp_rate=append(all_roc_fp_rate,fp/(fp+tn))
  #   all_roc_tp_rate=append(all_roc_tp_rate,tp/(tp+fn))
  # }
  #
  # #Compute F1 score
  # tp=sum(tp_scores)
  # fp=sum(fp_scores)
  # fn=sum(fn_scores)
  # p=tp/(tp+fp) #precision
  # r=tp/(tp+fn) #recall
  # f1= 2*(r*p)/(r+p) #F1 score
  # print(f1)

  #FIG 2A Scattter S1 vs S2
  tiff(paste(output_dir,"Fig2A_Scatter_S1_S2.tiff",sep=""), width = 100, height = 100, units = 'mm', res = 300)

  df <- data.frame(x=s1_scores,y=s2_scores,col=as.factor(colors))
  p_a<-ggplot(df,aes(x,y,color=col))+ geom_point() +theme_bw()+xlab("S1")+ylab("S2")+
    scale_color_discrete("Ground truth", c("Negative class", "Positive class"),breaks=c(0,1))+
    scale_x_continuous(breaks = seq(0, 1, by = 0.2))+scale_y_continuous(breaks = seq(0, 1, by = 0.2))+
    theme(legend.justification = c(1, 0), legend.position = c(1, 0),legend.background=element_rect(size=0.2,colour=1),panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0))
  dev.off()

  #FIG 2B Precision Recall Curves
  tiff(paste(output_dir,"Fig2B_Precision_Recall_Curves.tiff",sep=""), width = 100, height = 100, units = 'mm', res = 300)

  prc1=PRROC::pr.curve(scores.class0 =s1_scores ,weights.class0 =is.element(clusters_list,truth) ,curve = T)
  prc2=PRROC::pr.curve(scores.class0 =s2_scores ,weights.class0 =is.element(clusters_list,truth) ,curve = T)
  prc_product=PRROC::pr.curve(scores.class0 =s1_scores*s2_scores*1 ,weights.class0 =is.element(clusters_list,truth) ,curve = T)

  df1<-data.frame(prc1$curve)
  df2<-data.frame(prc2$curve)
  df_product<-data.frame(prc_product$curve)
  df1$Score=paste("S1 (",round(prc1$auc.integral,2),"AUC )")
  df2$Score=paste("S2 (",round(prc2$auc.integral,2),"AUC )")
  df_product$Score=paste("S1·S2 (",round(prc_product$auc.integral,2),"AUC )")

  df<-rbind(df1,df2,df_product)

  p_b <- ggplot(df) + geom_line(aes(X1, X2, colour=Score))+
    theme_bw()+xlab("Recall")+ylab("Precision")+
    #scale_color_discrete("", c(paste("S1 (",round(prc1$auc.integral,2),"AUC )"), paste("S2 (",round(prc2$auc.integral,2),"AUC )"),paste("S1·S2 (",round(prc_product$auc.integral,2),"AUC )")))+
    scale_x_continuous(breaks = seq(0, 1, by = 0.2))+scale_y_continuous(breaks = seq(0, 1, by = 0.2))+
    theme(legend.justification = c(0, 0), legend.position = c(0, 0),legend.background=element_rect(size=0.2,colour=1),panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0))
  print(p_b)

  cols=c("#F8766D","#00BA38","#619CFF")

  # df <- data.frame(x=s1_scores,y=s2_scores,col=as.factor(colors))
  # p<-ggplot(df,aes(x,y,color=col))+ geom_point() +theme_bw()+xlab("S1")+ylab("S2")+
  #   #scale_color_discrete("Ground truth", c("Negative class", "Positive class"),breaks=c(0,1))+
  #   scale_x_continuous(breaks = seq(0, 1, by = 0.2))+scale_y_continuous(breaks = seq(0, 1, by = 0.2))+
  #   theme(legend.justification = c(1, 0), legend.position = c(1, 0))
  # print(p)

  # plot(prc1$curve,col=cols[1],type="l",xlim=c(0,1),ylim=c(0,1),xlab="Recall",ylab="Precision",asp=1)
  # lines(prc2$curve, add=T, col=cols[2])
  # lines(prc_product$curve, add=T, col=cols[3])
  #
  # legend=paste("S1 (",round(prc1$auc.integral,2),"AUC )")
  # legend=append(legend,paste("S2 (",round(prc2$auc.integral,2),"AUC )"))
  # legend=append(legend,paste("S1 * S2 (",round(prc_product$auc.integral,2),"AUC )"))
  # legend("bottomleft",legend=legend,col=cols,lty=1,bg="white")

  dev.off()


  #Fig 2C Scatter Datasets vs S1 * S2
  tiff(paste(output_dir,"Fig2C_Scatter_S1S2_Clusters_1.tiff",sep=""), width = 200, height = 100, units = 'mm', res = 300)

  unique_clusters=unique(clusters_list)
  scores=s1_scores*s2_scores
  mean_scores=NULL
  for(clus in unique_clusters)
  {
    mean_scores=append(mean_scores,mean(scores[which(clusters_list==clus)]))
  }
  sorted_mean_scores=sort(mean_scores,decreasing = T,index.return=T)

  unique_ids=1:length(unique_clusters)
  names(unique_ids)=unique_clusters[sorted_mean_scores$ix]
  labels=c("Ag1","Ag3","Ag6","Ag10","Ag1C28H58O1")
  labels_expresion=c("Ag","Ag[3]","Ag[6]","Ag[10]","paste(Ag,C[28],H[58],O)")
  #gsub("([0-9]+)", "[\\1] ", labels)
  indices=which(is.element(clusters_list,labels))
  lab=rep("",length(clusters_list))
  lab[indices]=clusters_list[indices]

  df <- data.frame(x=unique_ids[clusters_list],y=scores,z=as.factor(colors), lab=rep("  \n   ",length(clusters_list)))
  p_c<-ggplot(df,aes(label=lab))  + geom_vline(xintercept = unique_ids[labels],linetype="dashed",alpha=0.5)+
    annotate("text_repel", x = unique_ids[labels], y=0.6, label = labels_expresion,box.padding=3,size=4,parse=TRUE)+
    geom_point(aes(x, y, color=z))+
    theme_bw()+xlab("Cluster number")+ylab("S1·S2")+
    scale_color_discrete("Ground truth", c("Negative class", "Positive class"),breaks=c(0,1))+
    scale_x_continuous(breaks = seq(0, 85, by = 5))+scale_y_continuous(breaks = seq(0, 1, by = 0.2))+
    theme(legend.justification = c(1, 1), legend.position = c(1, 1),legend.background=element_rect(size=0.2,colour=1),panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0))
  print(p_c)
  #plot(unique_ids[clusters_list],scores,col=colors,xlab="Cluster",ylab="S1·S2",pch=16,cex=0.5)

  dev.off()

  tiff(paste(output_dir,"Fig2C_Scatter_S1S2_Clusters_2.tiff",sep=""), width = 200, height = 100, units = 'mm', res = 300)

  labels=names(unique_ids)
  ticks=seq(from=1,to = length(unique_clusters),by=5)
  labels=labels[ticks]

  plot(unique_ids[clusters_list],scores,col=colors,xlab="Cluster",ylab="S1·S2",pch=16,cex=0.5,xaxt='n')
  axis(1, at=ticks, labels=FALSE)
  text(ticks, par("usr")[3]-0.05-0.01*nchar(labels), labels = labels, srt = 45, pos = 1, adj=0, xpd = TRUE, cex=0.7)
  dev.off()

  tiff(paste(output_dir,"Fig2C_Scatter_S1S2_Clusters_3.tiff",sep=""), width = 200, height = 100, units = 'mm', res = 300)

  labels=c("Ag1","Ag3","Ag6","Ag10","Ag1C28H58O1")
  ticks=unique_ids[labels]

  plot(unique_ids[clusters_list],scores,col=colors,xlab="Cluster",ylab="S1·S2",pch=16,cex=0.5,xaxt='n')
  axis(1, at=ticks, labels=FALSE)
  text(ticks, par("usr")[3]-0.05-0.01*nchar(labels), labels = labels, srt = 45, pos = 1, adj=0, xpd = TRUE, cex=0.7)
  dev.off()

  #Fig 2D Scatter Clusters vs S1 * S2
  tiff(paste(output_dir,"Fig2D_Scatter_S1S2_Datasets.tiff",sep=""), width = 200, height = 100, units = 'mm', res = 300)

  ordered_datasets=c(11,12,1,2,3,4,5,6,7,8,9,10,13,14)

  plot(ordered_datasets[dataset_list],scores,col=colors,xlab="Dataset",ylab="S1·S2",pch=16,cex=0.5)

  dev.off()

  #Fig 2 Composition
  tiff(paste(output_dir,"Fig2.tiff",sep=""), width = 200, height = 200, units = 'mm', res = 300)
  print(grid.arrange(p_a,p_b,p_c,layout_matrix=matrix(c(1,3,2,3),nrow=2)))
  dev.off()


  #S1 vs. S2
  #plot(s1_scores,s2_scores,col=colors,cex.lab=2, cex.axis=2,cex=1,xlab="S1",ylab="S2")
  #text(s1_scores, s2_scores, labels=labels, cex= 0.7, pos=3,col=colors)
  #legend("topleft",legend=c("Not in peak matrix","Should not be present","Should be matrix-related"),col=1:3,pch = 1)
  #abline(v = results$meta$s1_threshold)
  #abline(h = results$meta$s2_threshold)

  #plot(s1_scores,s2_scores,xlim = c(results$meta$s1_threshold,1),ylim=c(results$meta$s2_threshold,1),col=colors)
  #text(s1_scores, s2_scores, labels=labels, cex= 0.7, pos=3,col=colors)
  #legend("topleft",legend=c("Not in peak matrix","Should not be present","Should be matrix-related"),col=1:3,pch = 1)
  #abline(v = results$meta$s1_threshold)
  #abline(h = results$meta$s2_threshold)

  #S1 vs. S3
  #plot(s1_scores,s3_scores,col=colors,cex.lab=2, cex.axis=2,cex=1,xlab="S1",ylab="S3")
  #text(s1_scores, s3_scores, labels=labels, cex= 0.7, pos=3,col=colors)
  #legend("topleft",legend=c("Not in peak matrix","Should not be present","Should be matrix-related"),col=1:3,pch = 1)
  #abline(v = results$meta$s1_threshold)
  #abline(h = results$meta$s3_threshold)

  # plot(s1_scores,s3_scores,xlim = c(results$meta$s1_threshold,1),ylim=c(results$meta$s3_threshold,1),col=colors)
  # text(s1_scores, s3_scores, labels=labels, cex= 0.7, pos=3,col=colors)
  # legend("topleft",legend=c("Not in peak matrix","Should not be present","Should be matrix-related"),col=1:3,pch = 1)
  # abline(v = results$meta$s1_threshold)
  # abline(h = results$meta$s3_threshold)

  #Common false negatives
  # cluster_indices=which(apply(classification_list=="fn",2,sum)>0)
  # if(length(cluster_indices))
  # {
  #   fn_matrix=(classification_list[,cluster_indices]=="fn")*2+(classification_list[,cluster_indices]=="np")*1-1
  #   colnames(fn_matrix)<-global_clusters[cluster_indices]
  #   rownames(fn_matrix)<-1:dim(fn_matrix)[1]
  #   print(levelplot(t(fn_matrix),xlab="Clusters",ylab="Dataset",main="False Negatives"))
  # }
  # #Common false positives
  # cluster_indices=which(apply(classification_list=="fp",2,sum)>0)
  # if(length(cluster_indices))
  # {
  #   fn_matrix=(classification_list[,cluster_indices]=="fp")*2+(classification_list[,cluster_indices]=="np")*1-1
  #   colnames(fn_matrix)<-global_clusters[cluster_indices]
  #   rownames(fn_matrix)<-1:dim(fn_matrix)[1]
  #   print(levelplot(t(fn_matrix),xlab="Clusters",ylab="Dataset",main="False Positives"))
  # }
  #Close pdf file
  #dev.off()
}
