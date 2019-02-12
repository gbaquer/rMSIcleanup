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
generate_gt <- function (matrix_formula,pks,matching_method="strict",cor_threshold=0.85,generate_pdf=F) {

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
  pks$intensity=pks$intensity/pks$normalizations$TIC

  #0. Load data
  isotopes=NULL
  data("isotopes", package = "enviPat", envir = environment())
  adducts=NULL
  data("adducts", package = "enviPat", envir=environment())

  #1. Determine adducts depending on matrix_formula
  adducts_list=c("")
  #if(matrix_formula=="Ag1")
  #adducts_list=c("","Cl1","N1O3")
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

  #4. Open pdf report
  if(generate_pdf)
  {
    #Generate pdf name [pks$name_001]
    pdf_file=paste("output/",pks_Ag$names[1],"_000.pdf",sep="")
    i=0
    while(file.exists(pdf_file))
    {
      i=i+1
      pdf_file=paste("output/",pks_Ag$names[1],"_",str_pad(i, 3, pad = "0"),".pdf",sep="")
    }
    #open pdf file
    a4_width=8.27
    a4_height=11.69
    pdf(pdf_file,paper='a4',width=a4_width,height=a4_height)
    #First page of metadata [File name, mean image, matrix formula, adducts list, base forms, cor_threshold]
    text=NULL
    text=append(text,"###################################################################")
    text=append(text,"IMAGE INFORMATION")
    text=append(text,paste(strwrap(paste("- File_name:",pks$names[1])),collapse = "\n"))
    text=append(text,paste(strwrap(paste("- Number of peaks:",length(pks$mass))),collapse = "\n"))
    text=append(text,paste(strwrap(paste("- Number of pixels:",pks$numPixels[1])),collapse = "\n"))
    text=append(text,paste(strwrap(paste("- Mass Range: [",min(pks$mass),", ", max(pks$mass),"]")),collapse = "\n"))
    text=append(text,paste(strwrap(paste("- Masses:",paste(pks$mass,collapse="; "))),collapse = "\n"))
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
    if(nrow(checked)==0)
      break;
    patterns=isowrap(isotopes,checked,resmass = FALSE,resolution = 20000) #Pere said between 20000 and 25000 #Used rule of thumb FWHM_res=ppm*32
    #Append to final list
    for(i in 1:length(patterns))
    {
      max_mass[i]=max(patterns[[i]][,1])
      patterns_out$mass=append(patterns_out$mass,patterns[[i]][,1])
      patterns_out$intensity=append(patterns_out$intensity,patterns[[i]][,2])
      patterns_out$cluster=append(patterns_out$cluster,rep(attributes(patterns)$names[i],length(patterns[[i]][,1])))
    }
    clus_num=clus_num+1
  }
  #5. Determine matches with masses
  tol=200e-6
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

      rel_error=abs(image_mass-cluster_mass)/image_mass
      image_intensity[which(rel_error>tol)]=0

      correl=cor(cluster_intensity,image_intensity)
      if(length(cluster_index)==1)
        correl=1

      correlations=append(correlations,correl)


      image_correl_matrix=cor(pks$intensity[,image_index])
      image_correl = round(apply(image_correl_matrix,2,mean),digits = 2)
      image_correl_labels = append(rep("",length(image_correl)),image_correl)


      if(!is.na(correl))
      {
        #[Include a metric that checks for correlation of the image]
        if(correl>cor_threshold)
        {
          new_gt=image_mass[which(apply(abs(outer(cluster_mass,image_mass,'-')),2,min)/image_mass<tol)]
          print(new_gt)
          gt=append(gt,new_gt)
        }
      }
      #Print
      if(pkg_opt()$verbose_level<=-1 || generate_pdf)
      {
        if(max(image_intensity)!=0)
          image_intensity=image_intensity/max(image_intensity)
        if(max(cluster_intensity)!=0)
          cluster_intensity=cluster_intensity/max(cluster_intensity)

        print(paste(c,correl))
        print(mean_image[image_index])
        df=data.frame(mass=cluster_mass,cluster_intensity=cluster_intensity,image_intensity=image_intensity)
        melted_df=melt(df, id.vars="mass",measure.vars=c("cluster_intensity","image_intensity"))
        value=NULL
        variable=NULL
        mass=NULL
        # plot1= ggplot(melted_df, aes(mass,value,color=variable,ymin=0,ymax=value) ) + geom_linerange() + geom_point() + geom_text(aes(label=image_correl))
        # plot1=plot1 + ggtitle(paste(c,correl)) + xlab("m/z") + ylab("Rel Intensity")
        # print(grid.arrange(plot1,plot1,plot1,plot1,layout_matrix=page_layout))

        plt_correl= ggplot(melted_df, aes(mass,value,color=variable,ymin=0,ymax=value) ) + geom_linerange() + geom_point() + geom_text(aes(label=image_correl_labels))
        plt_correl=plt_correl + ggtitle(paste(c,round(correl,2),round(mean(image_correl),2))) + xlab("m/z") + ylab("Rel Intensity")
        plt_example=plt_correl
        text=paste(apply(round(image_correl_matrix,2),1,paste,collapse=" "),collapse="\n")
        plt_text=plot_text(text)
        plt_text=levelplot(image_correl_matrix)
        plts=list(plt_correl,plt_text)
        for(i in image_index)
          plts=c(plts,list(ggplot_peak_image(pks,pks$intensity[,i],paste("m/z",pks$mass[i]))))


        #Generate images
        #for(i in 1:length)
        grid.arrange(grobs=plts,layout_matrix=page_layout)

        #grid.arrange(plot(cluster_mass,cluster_intensity),plot(image_mass,image_intensity))
        print("---")
      }
    }
  }


  #Each page

  # Close file
  if(generate_pdf)
  {
    dev.off()
  }

  return(unique(gt))
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
