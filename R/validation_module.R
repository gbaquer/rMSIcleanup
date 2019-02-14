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

#' Generate gold standard
#'
#' Given the matrix formula and the matrix peaks present in the image it returns the ground truth defined as the list of peaks present in the image that correspond to the matrix.
#'
#' @param matrix_formula String giving the chemical formula of the matrix in the enviPat notation
#' @param pks Peak Matrix Image
#' @param matching_method "loose": All peaks within a given tolerance are considered matrix peaks; "strict": Only peaks within a given tolerance that mantain the abundance ratios for each cluster are considered matrix peaks
#' @param cor_threshold Correlation between the theoretical and the real spectral pattern above which a given cluster is considered to be present and thus included in the gold truth (gt)
#' @param generate_pdf Boolean indicating whether to generate a pdf or not
#' @return Ground Truth: List of masses available in the image that correspond to the matrix.
#'
#' @export
generate_gt <- function (matrix_formula,pks,matching_method="strict",cor_threshold=0.85,generate_pdf=F,clus_num=1) {

  #SECTION 0 :: Preprocessing
  #Load semi-raw

  # Select first image if there are multiple
  if(length(pks$numPixels)>1)
  {
    rows=1:pks$numPixels[1]
    #rows=(pks$numPixels[1]+1):(pks$numPixels[1]+pks$numPixels[2])

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

  #Load semiraw
  semi_raw <- rMSI::LoadMsiData("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/comparativa_Ag_Au/postcode3/20160801-BrainCPF-Ag-CPF-proc.tar")


  #[CONTINUE HERE]
  # i=which.min(abs(semi_raw$mass-mass))
  # cols=(i-tol_scans):(i+tol_scans)
  # mass_range=semi_raw$mass[cols]
  # match(c(3,5),b)

  # Normalize to TIC
  pks$intensity=pks$intensity/pks$normalizations$TIC

  #Select one cluster
  # clus <- kmeans(pks$intensity, centers = 3)
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

  #0. Load data
  isotopes=NULL
  data("isotopes", package = "enviPat", envir = environment())
  adducts=NULL
  data("adducts", package = "enviPat", envir=environment())

  #1. Determine adducts depending on matrix_formula
  adducts_list=c("")
  sub_list=NULL
  #if(matrix_formula=="Ag1")
  #adducts_list=c("","Cl1","N1O3")
  # adducts_list=c("","H1","K1","Na1")
  # sub_list=c("H2O1")
  #For more complex matrix formulas adducts could be loaded from the library
  #adducts_formula=adducts$Formula_add[1:4]

  #2. [Not implemented yet] Assess which mass range is needed

  #Check max mass of the highest isotope and highest adduct. Divide and conquer.
  # for(elem in union())
  # tail(grep(elem,isotopes$element),1)

  #3.Generate list of possible chemical formulas
  base_forms=NULL
  for(i in 1:1)
  {
    base_forms=append(base_forms,paste(multiform(matrix_formula,i),adducts_list,sep=""))
  }
  if(!is.null(sub_list))
    base_forms=append(base_forms,unlist(lapply(sub_list,function(x) enviPat::subform(base_forms,x))))


  #bIPASS
  # base_forms=c("Ag9")
  #4. Open pdf report
  if(generate_pdf)
  {
    #Generate pdf name [pks$name_001]
    pdf_file=paste("output/",pks$names[1],"_000.pdf",sep="")
    i=0
    while(file.exists(pdf_file))
    {
      i=i+1
      pdf_file=paste("output/",pks$names[1],"_",str_pad(i, 3, pad = "0"),".pdf",sep="")
    }
    #open pdf file
    a4_width=8.27
    a4_height=11.69
    pdf(pdf_file,paper='a4',width=a4_width,height=a4_height)
    #First page of metadata [File name, mean image, matrix formula, adducts list, base forms, cor_threshold]
    text=NULL
    text=append(text,as.character(Sys.time()))
    text=append(text,"###################################################################")
    text=append(text,"IMAGE INFORMATION")
    text=append(text,paste(strwrap(paste("- File_name:",pks$names[1])),collapse = "\n"))
    text=append(text,paste(strwrap(paste("- Number of peaks:",length(pks$mass))),collapse = "\n"))
    text=append(text,paste(strwrap(paste("- Number of pixels:",pks$numPixels[1])),collapse = "\n"))
    text=append(text,paste(strwrap(paste("- Mass Range: [",min(pks$mass),", ", max(pks$mass),"]")),collapse = "\n"))
    #text=append(text,paste(strwrap(paste("- Masses:",paste(pks$mass,collapse="; "))),collapse = "\n"))
    text=append(text,"###################################################################")
    text=append(text,"MATRIX INFORMATION")
    text=append(text,paste(strwrap(paste("- Matrix formula:",matrix_formula)),collapse = "\n"))
    text=append(text,paste(strwrap(paste("- Adducts list:",paste(adducts_list,collapse="; "))),collapse = "\n"))
    text=append(text,paste(strwrap(paste("- Base forms:",paste(base_forms,collapse="; "))),collapse = "\n"))
    text=append(text,paste(strwrap(paste("- Correlation threshold:",cor_threshold)),collapse = "\n"))
    text=append(text,"###################################################################")

    text=paste(text,collapse = "\n")
    print(plot_text(text))
  }

  #5.Generate pattern list with enviPat
  page_layout=rbind(c(1,1,2,2),
                    c(1,1,2,2),
                    c(3,4,5,6),
                    c(7,8,9,10),
                    c(11,12,13,14),
                    c(15,16,17,18))
  patterns_out=NULL
  clus_num=1
  max_mass=rep(-1,length(base_forms))
  MALDI_resolution=cbind(c(1040.189125,1295.508274,1342.789598,1607.565012,2089.834515,2468.085106,3148.93617,4548.463357),c(26012.14575,34514.17004,36437.24696,41497.97571,44939.27126,44534.41296,42510.12146,37044.53441))
  dimnames(MALDI_resolution)[[2]]=c("m/z","R")
  while(min(max_mass)<max(pks$mass))
  {
    forms=multiform(base_forms,clus_num)
    #patterns=isopattern(isotopes,forms) [Improvement: The ppm threshold doesn't work]
    checked=check_chemform(isotopes,forms)
    checked=checked[which(checked$monoisotopic_mass<max(pks$mass)),]
    print(checked)
    if(nrow(checked)==0)
      break;
    patterns=isowrap(isotopes,checked,resmass = FALSE,resolution = 20000) #Pere said between 20000 and 25000 #Used rule of thumb FWHM_res=ppm*32
    #Append to final list
    for(i in 1:length(patterns))
    {
      max_mass[i]=max(patterns[[i]][,1])
      if(!is.element(attributes(patterns)$names[i],patterns_out$cluster))
      {
        patterns_out$mass=append(patterns_out$mass,patterns[[i]][,1])
        patterns_out$intensity=append(patterns_out$intensity,patterns[[i]][,2])
        patterns_out$cluster=append(patterns_out$cluster,rep(attributes(patterns)$names[i],length(patterns[[i]][,1])))
      }
    }
    clus_num=clus_num+1
  }
  #5. Determine matches with masses
  tol=500e-6
  gt=NULL
  if(matching_method=="loose")
  {
    #All peaks within a given tolerance are considered matrix peaks
    gt=pks$mass[which(apply(abs(outer(patterns_out[,1],pks$mass,'-')),2,min)/pks$mass<tol)]
  }
  if(matching_method=="strict")
  {
    #Only peaks within a given tolerance that mantain the abundance ratios for each cluster are considered matrix peaks
    mean_image=apply(pks$intensity,2,mean)
    clusters=unique(patterns_out$cluster)
    correlations=NULL
    for(c in clusters)
    {
      cluster_index=which(patterns_out$cluster==c)
      cluster_mass=patterns_out$mass[cluster_index]
      cluster_intensity=patterns_out$intensity[cluster_index]

      image_index=apply(abs(outer(cluster_mass,pks$mass,'-')),1,function(x) sort(x,index.return=TRUE)$ix[1])#which(apply(abs(outer(cluster_mass,pks$mass,'-')),2,min)/pks$mass<tol)
      image_mass=pks$mass[image_index]
      image_intensity=mean_image[image_index]

      final_image=NULL



      rel_error=abs(image_mass-cluster_mass)/cluster_mass
      # [Most likely wrong due to normalization]

      index_pksMat=which(rel_error<=tol)
      index_spec=which(rel_error>tol)
      for(i in 1:length(cluster_index))
      {
        if(rel_error>tol)
        {
          mass=cluster_mass[i]
          data=rMSI::loadImageSliceFromMass(semi_raw, mass,mass*tol/10)$data/pks$normalizations$TIC
          max_col=which.max(apply(data,2,mean))
          data_picture=data[,max_col]
          image_intensity[i]=mean(data_picture)

        }
        else
        {
          data_picture=pks$intensity[,image_index[i]]
        }
        final_image=cbind(final_image,data_picture)
      }
      colnames(final_image)<-NULL
      #image_intensity[which(rel_error>tol)]=0

      correl=cor(cluster_intensity,image_intensity)
      if(length(cluster_index)==1)
        correl=1

      #image_intensity[which(rel_error>tol)]=NA

      correlations=append(correlations,correl)


      # cluster_images=pks$intensity[,image_index]
      # cluster_images[,which(rel_error>tol)]=NA

      image_correl_matrix=cor(final_image)
      image_correl_matrix[which(is.na(image_correl_matrix))]=0
      image_correl = apply(image_correl_matrix,2,weighted.mean,w=cluster_intensity)
      image_correl_labels = append(rep("",length(image_correl)),round(image_correl,2))

      chosen=rep(F,length(image_mass))
      if(!is.na(correl))
      {
        #[Include a metric that checks for correlation of the image]
        if(correl>cor_threshold)
        {
          # chosen=apply(abs(outer(cluster_mass,image_mass,'-')),2,min)/image_mass<tol
          # new_gt=image_mass[which(chosen)]
          chosen[index_pksMat]=T
          new_gt=image_mass[index_pksMat]
          gt=append(gt,new_gt)
        }
      }
      #Print
      if(pkg_opt()$verbose_level<=-1 || generate_pdf)
      {
        intensities=image_intensity[which(!is.na(image_intensity))]
        if(length(intensities)!=0 && max(intensities)!=0)
          norm_image_intensity=image_intensity/max(intensities)
        if(max(cluster_intensity)!=0)
          cluster_intensity=cluster_intensity/max(cluster_intensity)

        print(paste(c,correl))
        df=data.frame(mass=cluster_mass,cluster_intensity=cluster_intensity,image_intensity=norm_image_intensity)
        melted_df=melt(df, id.vars="mass",measure.vars=c("cluster_intensity","image_intensity"))
        value=NULL
        variable=NULL
        mass=NULL

        s1=correl
        image_correl[which(is.na(image_correl))]=0
        s2=weighted.mean(image_correl,cluster_intensity)

        plt_correl= ggplot(melted_df, aes(mass,value,color=variable,ymin=0,ymax=value) ) + geom_linerange() + geom_point() + geom_text(aes(label=image_correl_labels,vjust=0))
        plt_correl= plt_correl + ggtitle(paste(c,"; S1:",round(s1,2),"; S2:",round(s2,2),"; MAX:", round(max(image_intensity),2))) + xlab("m/z") + ylab("Normalized Intensity") +scale_colour_discrete(name="",breaks=c("cluster_intensity","image_intensity"),labels=c("Calculated m/z","Experimental m/z"))+ theme(legend.position="bottom")

        text=paste(apply(round(image_correl_matrix,2),1,paste,collapse=" "),collapse="\n")
        plt_text=plot_text(text)
        plt_text=levelplot(image_correl_matrix,at=seq(-1,1,length.out = 100))
        plts=list(plt_correl,plt_text)

        page_layout=rbind(c(1,1,2,2),
                          c(1,1,2,2),
                          c(3,4,5,6),
                          c(7,8,9,10),
                          c(11,12,13,14),
                          c(15,16,17,18))
        offset=match(3,c(t(page_layout)))-1
        for(i in 1:length(image_index))
        {
          # print(pks$intensity[,image_index[i]])
          # print(paste("m/z",round(cluster_mass[i],4)," +/-",round(rel_error[i]*1e6),"ppm"))
          # print(is.na(image_intensity[i]))
          # print(chosen[i])
          title=paste("m/z",round(cluster_mass[i],4))
          if(is.element(i,index_pksMat))
            title=paste(title," +/-",round(rel_error[i]*1e6),"ppm")
          else
            title=paste(title,"(NOT IN PEAK MATRIX)")
          plts=c(plts,list(ggplot_peak_image(pks,final_image[,i],title,is.na(image_intensity[i]),chosen[i],is.element(i,index_pksMat))))
          if((i+offset)%%length(page_layout)==0||i==length(image_index))
          {
            grid.arrange(grobs=plts,layout_matrix=page_layout)
            plts=list()
            offset=0
            page_layout=matrix(1:length(page_layout), ncol = 4)
          }
        }

        #Generate images
        #for(i in 1:length)
        #grid.arrange(grobs=plts,layout_matrix=page_layout)

        #grid.arrange(plot(cluster_mass,cluster_intensity),plot(image_mass,image_intensity))
        print("---")
      }
    }
  }

  gt=unique(gt)
  not_gt=setdiff(pks$mass,gt)

  #Generate putative clusters

  # Close file
  if(generate_pdf)
  {
    if(length(gt)>0)
    {
      plts=list()
      gt_index=match(gt,pks$mass)
      for(i in 1:length(gt))
      {
        plts=c(plts,list(ggplot_peak_image(pks,pks$intensity[,gt_index[i]],paste("m/z",round(gt[i],4)),chosen=T)))
        if(i%%length(page_layout)==0||i==length(gt))
        {
          grid.arrange(grobs=plts,nrow=nrow(page_layout),ncol=ncol(page_layout))
          plts=list()
        }
      }
    }

    if(length(not_gt)>0)
    {
      plts=list()
      not_gt_index=match(not_gt,pks$mass)
      for(i in 1:length(not_gt))
      {
        plts=c(plts,list(ggplot_peak_image(pks,pks$intensity[,not_gt_index[i]],paste("m/z",round(not_gt[i],4)),chosen=F)))
        if(i%%length(page_layout)==0||i==length(not_gt))
        {
          grid.arrange(grobs=plts,nrow=nrow(page_layout),ncol=ncol(page_layout))
          plts=list()
        }
      }
    }

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
  tp_list=intersect(gt,pos)#pos[which(apply(abs(outer(gt,pos,'-')),2,min)/pos<tol)]
  fp_list=setdiff(pos,tp_list)#which(!is.element(pos,tp_list))
  fn_list=intersect(gt,neg)#neg[which(apply(abs(outer(gt,neg,'-')),2,min)/neg<tol)]
  tn_list=setdiff(neg,fn_list)#which(!is.element(neg,fn_list))

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
