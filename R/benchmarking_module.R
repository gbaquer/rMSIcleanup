#########################################################################
#
#     BENCHMARKING MODULE
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

#' Run experiment
#'
#' Run a given experiment
#'
#' @return None
#' @param base_dir Base directory where the images are stored and the output reports should be stored
#' @param dataset_indices Vector containing the indices of the datasets to be used
#' @inheritParams generate_gt
#'
#'
#' @export
run_experiment <- function (matrix_formula, base_dirs=c("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/Ag Software Test 1","/home/gbaquer/msidata/Ag Software Test 1"),
                            s1_threshold=0.80,s2_threshold=0.80, s3_threshold=0.7, similarity_method="euclidean", correlation_method="pearson", experiment_name="output",
                            MALDI_resolution=20000, tol_mode="scans",tol_ppm=200e-6,tol_scans=4,
                            mag_of_interest="intensity",normalization="None",
                            max_multi=10, add_list=NULL, sub_list=NULL, isobaric_detection=T,
                            save_results=T,generate_pdf=T,default_page_layout=NULL,include_summary=F,dataset_indices=NULL) {
  #0. Prepare directories
  base_dir=base_dirs[which(dir.exists(base_dirs))][1]

  images_dir=paste(base_dir,"/images",sep="")
  output_dir=paste(base_dir,"/output",sep="")
  experiment_dir=paste(generate_file_name(paste("/",experiment_name,sep=""),extension = "",folder = output_dir),"/",sep="")
  dir.create(experiment_dir)

  #1. Generate experiment metadata file
  #[PENDING]

  #2. Run experiment
  full_spectrum_name_list=list.files(images_dir,pattern = "*-proc.tar",recursive = T)
  pks_name_list=list.files(images_dir,pattern = "*mergeddata-peaks*",recursive = T,full.names = T)
  subfolder_list=unlist(lapply(strsplit(pks_name_list,'/'),function(x) paste(x[1:length(x)-1],collapse="/")))
  num_files=length(pks_name_list)

  results=list(meta=list(s1_threshold=s1_threshold,s2_threshold=s2_threshold,s3_threshold=s3_threshold,file_names=NULL),data=list())

  j=1
  for(i in 1:num_files)
  {
    if(is.null(dataset_indices) || is.element(i,dataset_indices))
    {
      pks_name=pks_name_list[i]
      subfolder=subfolder_list[i]
      print("Start peak matrix:")
      print(pks_name)
      pks=rMSIproc::LoadPeakMatrix(pks_name)

      pks_i=1
      for(name in pks$names)
      {
        print("Start full_spectrum:")
        print(name)
        full_spectrum_name=paste(subfolder,"/",unlist(strsplit(name,".",fixed=T))[1],"-proc.tar",sep="")
        print(full_spectrum_name)

        if(file.exists(full_spectrum_name))
        {
          full_spectrum=rMSI::LoadMsiData(full_spectrum_name)

          #[Potential improvement: Use ... instead]
          results$data[[j]]= generate_gt(matrix_formula=matrix_formula,pks=pks,full_spectrum=full_spectrum,folder=experiment_dir,
                      s1_threshold=s1_threshold,s2_threshold=s2_threshold, s3_threshold=s3_threshold, similarity_method=similarity_method,correlation_method=correlation_method,
                      MALDI_resolution=MALDI_resolution, tol_mode=tol_mode,tol_ppm=tol_ppm,tol_scans=tol_scans,
                      mag_of_interest=mag_of_interest,normalization=normalization,
                      max_multi=max_multi, add_list=add_list, sub_list=sub_list, isobaric_detection=isobaric_detection,
                      generate_pdf=generate_pdf,default_page_layout=default_page_layout,include_summary=include_summary,pks_i = pks_i)
          results$meta$file_names=append(results$meta$file_names,full_spectrum_name)
          j=j+1
        }
        pks_i=pks_i+1
      }
    }
  }

  #Store results
  if(save_results)
  {
    results_file=generate_file_name("global_results",folder = experiment_dir,extension = ".rds")
    saveRDS(results, results_file)
  }

  #Print results
  generate_pdf(results,experiment_dir)

}

#' Generate PDF report
#'
#' Generate PDF report
#'
#'
#' @return None
#'
#'
#' @export
generate_pdf_from_file <- function () {
  base_dirs=c("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/Ag Software Test 1","/home/gbaquer/msidata/Ag Software Test 1")
  base_dir=base_dirs[which(dir.exists(base_dirs))][1]
  output_dir=paste(base_dir,"/output/",sep="")
  #experiment_dir=paste(generate_file_name("/output",extension = "",folder = output_dir),"",sep="")
  #Load Results
  results_path = file.choose()
  results = readRDS(results_path)
  #Print Results
  generate_pdf(results,output_dir,results_path)
}

#' Generate PDF report
#'
#' Generate PDF report
#'
#' @param results Results matrix
#'
#' @return None
#'
#'
#' @export
generate_pdf <- function (results,experiment_dir,info=experiment_dir) {
  #Open pdf
  #Generate pdf name [pks$name_001]
  pdf_file=generate_file_name("global_results",folder = experiment_dir)
  #open pdf file
  a4_width=8.27
  a4_height=11.69
  pdf(pdf_file,width=a4_height,height=a4_height)

  info=paste(info,paste(results$meta$file_names,collapse="\n"),sep="\n\n")
  print(plot_text(info))

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
  truth=c("Ag1","Ag2","Ag3","Ag4","Ag5","Ag6","Ag7","Ag8","Ag9","Ag10","Ag1Na1")
  global_clusters=unique(unlist(lapply(results$data,function(x)unique(x$patterns_out$cluster))))
  dataset_names=NULL
  for(r in results$data)
  {
    #Print calculated clusters
    clusters=unique(r$patterns_out$cluster)
    cluster_indices=match(clusters,r$patterns_out$cluster)

    s1=r$patterns_out$s1_scores[cluster_indices]
    s2=r$patterns_out$s2_scores[cluster_indices]
    s3=r$patterns_out$s3_scores[cluster_indices]

    labels=append(labels,paste(i,"_",clusters,sep="")[which(s3!=0)])

    is_na=which(is.na(s1)|is.na(s2)|is.na(s3))
    col=is.element(clusters,truth)*1+2
    col[is_na]=1
    colors=append(colors,col[which(s3!=0)])

    s1[is.na(s1)]=0
    s2[is.na(s2)]=0
    s3[is.na(s3)]=0

    s1_scores=append(s1_scores,s1[which(s3!=0)])
    s2_scores=append(s2_scores,s2[which(s3!=0)])
    s3_scores=append(s3_scores,s3[which(s3!=0)])

    clusters_list=append(clusters_list,clusters[which(s3!=0)])

    #Compute positives and negatives
    chosen=(s1>results$meta$s1_threshold&s2>results$meta$s2_threshold&s3>results$meta$s3_threshold)
    should_be_chosen=is.element(clusters,truth)

    tp=sum(chosen&should_be_chosen&s3!=0)
    fp=sum(chosen&!should_be_chosen&s3!=0)
    tn=sum(!chosen&!should_be_chosen&s3!=0)
    fn=sum(!chosen&should_be_chosen&s3!=0)

    classification=rep("np",length(global_clusters))
    global_index=match(clusters,global_clusters)
    classification[global_index[which(chosen&should_be_chosen&s3!=0)]]="tp"
    classification[global_index[which(chosen&!should_be_chosen&s3!=0)]]="fp"
    classification[global_index[which(!chosen&!should_be_chosen&s3!=0)]]="tn"
    classification[global_index[which(!chosen&should_be_chosen&s3!=0)]]="fn"
    classification_list=rbind(classification_list,classification)

    fp_rate=fp/(fp+tn)
    tp_rate=tp/(tp+fn)

    tp_scores=append(tp_scores,tp)
    fp_scores=append(fp_scores,fp)
    tn_scores=append(tn_scores,tn)
    fn_scores=append(fn_scores,fn)

    fp_rate_scores=append(fp_rate_scores,fp_rate)
    tp_rate_scores=append(tp_rate_scores,tp_rate)

    i=i+1

  }
  clusters=unique(results$data[[1]]$patterns_out$cluster)
  s1_roc_fp_rate=NULL
  s1_roc_tp_rate=NULL
  s2_roc_fp_rate=NULL
  s2_roc_tp_rate=NULL
  s3_roc_fp_rate=NULL
  s3_roc_tp_rate=NULL
  all_roc_fp_rate=NULL
  all_roc_tp_rate=NULL
  for(t in seq(1,-0.001,length=1000))
  {
    #S1 ROC
    chosen=s1_scores>t
    should_be_chosen=is.element(clusters_list,truth)

    tp=sum(chosen&should_be_chosen)
    fp=sum(chosen&!should_be_chosen)
    tn=sum(!chosen&!should_be_chosen)
    fn=sum(!chosen&should_be_chosen)

    s1_roc_fp_rate=append(s1_roc_fp_rate,fp/(fp+tn))
    s1_roc_tp_rate=append(s1_roc_tp_rate,tp/(tp+fn))

    #S2 ROC
    chosen=s2_scores>t
    should_be_chosen=is.element(clusters_list,truth)

    tp=sum(chosen&should_be_chosen)
    fp=sum(chosen&!should_be_chosen)
    tn=sum(!chosen&!should_be_chosen)
    fn=sum(!chosen&should_be_chosen)

    s2_roc_fp_rate=append(s2_roc_fp_rate,fp/(fp+tn))
    s2_roc_tp_rate=append(s2_roc_tp_rate,tp/(tp+fn))

    #S3 ROC
    chosen=s3_scores>t
    should_be_chosen=is.element(clusters_list,truth)

    tp=sum(chosen&should_be_chosen)
    fp=sum(chosen&!should_be_chosen)
    tn=sum(!chosen&!should_be_chosen)
    fn=sum(!chosen&should_be_chosen)

    s3_roc_fp_rate=append(s3_roc_fp_rate,fp/(fp+tn))
    s3_roc_tp_rate=append(s3_roc_tp_rate,tp/(tp+fn))

    #ALL ROC
    chosen=s1_scores>t&s2_scores>t&s3_scores>t
    should_be_chosen=is.element(clusters_list,truth)

    tp=sum(chosen&should_be_chosen)
    fp=sum(chosen&!should_be_chosen)
    tn=sum(!chosen&!should_be_chosen)
    fn=sum(!chosen&should_be_chosen)

    all_roc_fp_rate=append(all_roc_fp_rate,fp/(fp+tn))
    all_roc_tp_rate=append(all_roc_tp_rate,tp/(tp+fn))
  }

  #Compute F1 score
  tp=sum(tp_scores)
  fp=sum(fp_scores)
  fn=sum(fn_scores)
  p=tp/(tp+fp) #precision
  r=tp/(tp+fn) #recall
  f1= 2*(r*p)/(r+p) #F1 score
  print(f1)

  #ROC AUC using pROC
  roc1=pROC::roc(response=is.element(clusters_list,truth),predictor=s1_scores)
  roc2=pROC::roc(response=is.element(clusters_list,truth),predictor=s2_scores)
  roc3=pROC::roc(response=is.element(clusters_list,truth),predictor=s3_scores)
  roc_product=pROC::roc(response=is.element(clusters_list,truth),predictor=s1_scores*s2_scores*1)
  roc_sum=pROC::roc(response=is.element(clusters_list,truth),predictor=s1_scores+s2_scores+0)
  roc_euclidean=pROC::roc(response=is.element(clusters_list,truth),predictor=sqrt(((1-s1_scores)^2)+((1-s2_scores)^2)+((1-1)^2)))
  roc_thresholds=pROC::roc(response=is.element(clusters_list,truth),predictor=s1_scores*(s2_scores>0.7)*(0.8>0.7))

  plot(1-roc1$specificities,roc1$sensitivities,col=2,type="l",xlim=c(0,1),ylim=c(0,1))
  lines(1-roc2$specificities,roc2$sensitivities, add=T, col=3)
  lines(1-roc3$specificities,roc3$sensitivities, add=T, col=4)
  lines(1-roc_product$specificities,roc_product$sensitivities, add=T, col=5)
  lines(1-roc_sum$specificities,roc_sum$sensitivities, add=T, col=6)
  lines(1-roc_euclidean$specificities,roc_euclidean$sensitivities, add=T, col=7)
  lines(1-roc_thresholds$specificities,roc_thresholds$sensitivities, add=T, col=7)
  lines(0:1,0:1)
  legend=paste("S1 (",round(roc1$auc,3),"AUC )")
  legend=append(legend,paste("S2 (",round(roc2$auc,3),"AUC )"))
  legend=append(legend,paste("S3 (",round(roc3$auc,3),"AUC )"))
  legend=append(legend,paste("Product (",round(roc_product$auc,3),"AUC )"))
  legend=append(legend,paste("Sum (",round(roc_sum$auc,3),"AUC )"))
  legend=append(legend,paste("Euclidean (",round(roc_euclidean$auc,3),"AUC )"))
  legend=append(legend,paste("Thresholds (",round(roc_thresholds$auc,3),"AUC )"))
  legend("bottomright",legend=legend,col=2:7,lty=1)

  # #ROC AUC
  # plot(s1_roc_fp_rate,s1_roc_tp_rate,col=2,type="l",xlim=c(0,1),ylim=c(0,1))
  # lines(s2_roc_fp_rate,s2_roc_tp_rate,col=3)
  # lines(s3_roc_fp_rate,s3_roc_tp_rate,col=4)
  # lines(all_roc_fp_rate,all_roc_tp_rate,col=5)
  # lines(0:1,0:1)
  # legend=paste("S1 (",round(mean(s1_roc_tp_rate),2),"AUC )")
  # legend=append(legend,paste("S2 (",round(mean(s2_roc_tp_rate),2),"AUC )"))
  # legend=append(legend,paste("S3 (",round(mean(s3_roc_tp_rate),2),"AUC )"))
  # legend=append(legend,paste("All (",round(mean(all_roc_tp_rate),2),"AUC )"))
  # legend=append(legend,"Reference (0.5 AUC)")
  # legend("bottomright",legend=legend,col=c(2,3,4,5,6,1),lty=1)

  #S1 vs. S2
  plot(s1_scores,s2_scores,col=colors)
  text(s1_scores, s2_scores, labels=labels, cex= 0.7, pos=3,col=colors)
  legend("topleft",legend=c("Not in peak matrix","Should not be present","Should be matrix-related"),col=1:3,pch = 1)
  abline(v = results$meta$s1_threshold)
  abline(h = results$meta$s2_threshold)

  plot(s1_scores,s2_scores,xlim = c(results$meta$s1_threshold,1),ylim=c(results$meta$s2_threshold,1),col=colors)
  text(s1_scores, s2_scores, labels=labels, cex= 0.7, pos=3,col=colors)
  legend("topleft",legend=c("Not in peak matrix","Should not be present","Should be matrix-related"),col=1:3,pch = 1)
  abline(v = results$meta$s1_threshold)
  abline(h = results$meta$s2_threshold)

  #S1 vs. S3
  plot(s1_scores,s3_scores,col=colors)
  text(s1_scores, s3_scores, labels=labels, cex= 0.7, pos=3,col=colors)
  legend("topleft",legend=c("Not in peak matrix","Should not be present","Should be matrix-related"),col=1:3,pch = 1)
  abline(v = results$meta$s1_threshold)
  abline(h = results$meta$s3_threshold)

  plot(s1_scores,s3_scores,xlim = c(results$meta$s1_threshold,1),ylim=c(results$meta$s3_threshold,1),col=colors)
  text(s1_scores, s3_scores, labels=labels, cex= 0.7, pos=3,col=colors)
  legend("topleft",legend=c("Not in peak matrix","Should not be present","Should be matrix-related"),col=1:3,pch = 1)
  abline(v = results$meta$s1_threshold)
  abline(h = results$meta$s3_threshold)

  #Common false negatives
  cluster_indices=which(apply(classification_list=="fn",2,sum)>0)
  fn_matrix=(classification_list[,cluster_indices]=="fn")*2+(classification_list[,cluster_indices]=="np")*1-1
  colnames(fn_matrix)<-global_clusters[cluster_indices]
  rownames(fn_matrix)<-1:dim(fn_matrix)[1]
  print(levelplot(t(fn_matrix),xlab="Clusters",ylab="Dataset",main="False Negatives"))

  #Common false positives
  cluster_indices=which(apply(classification_list=="fp",2,sum)>0)
  fn_matrix=(classification_list[,cluster_indices]=="fp")*2+(classification_list[,cluster_indices]=="np")*1-1
  colnames(fn_matrix)<-global_clusters[cluster_indices]
  rownames(fn_matrix)<-1:dim(fn_matrix)[1]
  print(levelplot(t(fn_matrix),xlab="Clusters",ylab="Dataset",main="False Positives"))
  #Close pdf file
  dev.off()
}
#' Test Paraffin
#'
#' Test generate_gt with Paraffin
#'
#' @return None
#'
#'
#' @export
test_paraffin <- function () {
  #Load images
  pks_Paraffin <- rMSIproc::LoadPeakMatrix("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/Parafina Software Test 1/images/20160719_TOF_Au_TumorFrescVsDespar/mergeddata-peaks.zip")
  full_spectrum_Paraffin1 <- rMSI::LoadMsiData("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/Parafina Software Test 1/images/20160719_TOF_Au_TumorFrescVsDespar/20160719-TumorDespar-Au_CAL-proc.tar")
  #full_spectrum_Paraffin2 <- rMSI::LoadMsiData("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/Parafina Software Test 1/images/20160719_TOF_Au_TumorFrescVsDespar/20160719-TumorFresc-Au_CAL-proc.tar")

  #Generate matrix formula
  n=10:100
  m=2*n+2
  matrix_formula=paste("C",n,"H",m,sep="")

  #Run generate gt
  rMSIcleanup::generate_gt(matrix_formula,pks_Paraffin,full_spectrum_Paraffin1,generate_pdf = T,max_multi = 1, normalization="TIC")
}


#' Cross validation.
#'
#' Compare the results of all methods to assess consistency.
#'
#' @return None
#'
#'
#' @export
cross_validation <- function () {
  #LOAD DATA
  #pks_Ag <- rMSIproc::LoadPeakMatrix("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/comparativa_Ag_Au/postcode3/20160801-BrainCPF-Ag-mergeddata-peaks.zip")
  pks_Ag <- rMSIproc::LoadPeakMatrix("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/Ag Software Test 1/images/20160801_TOF_AuAg_BrainCPF/mergeddata-peaks.zip")
  pks_Norharmane <- rMSIproc::LoadPeakMatrix("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/comparativa_matrix_au/peak_matrix_norharmane/mergeddata-peaks.zip")
  pks_Au <- rMSIproc::LoadPeakMatrix("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/comparativa_matrix_au/peak_matrix_au/mergeddata-peaks.zip")
  pks_WO3 <- rMSIproc::LoadPeakMatrix("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/kidney_auagwo3/postcode1/mergeddata-peaks.zip")
  pks_HCCA <- rMSIproc::LoadPeakMatrix("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/Lluc Images/postcode2/region7_calibrated/mergeddata-peaks.zip")
  pks_Garlic <- rMSIproc::LoadPeakMatrix("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/Garlic Alex/mergeddata-peaks.zip")
  pks_Paraffin <- rMSIproc::LoadPeakMatrix("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/Parafina Software Test 1/images/20180125_TOF_Au_Ghrelina_Parafina/mergeddata-peaks.zip")
  #LOAD Full spectra
  #full_spectrum_Ag <- rMSI::LoadMsiData("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/comparativa_Ag_Au/postcode3/20160801-BrainCPF-Ag-CPF-proc.tar")
  full_spectrum_Ag <- rMSI::LoadMsiData("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/Ag Software Test 1/images/20160801_TOF_AuAg_BrainCPF/CPF-Ag-proc.tar")
  full_spectrum_Paraffin1 <- rMSI::LoadMsiData("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/Parafina Software Test 1/images/20180125_TOF_Au_Ghrelina_Parafina/20180125_TOF_Au_Ghrelina-01C-proc.tar")
  full_spectrum_Paraffin2 <- rMSI::LoadMsiData("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/Parafina Software Test 1/images/20180125_TOF_Au_Ghrelina_Parafina/20180125_TOF_Au_Ghrelina-02A-proc.tar")
  full_spectrum_HCCA <- rMSI::LoadMsiData("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/Lluc Images/postcode2/region7_calibrated/18-3101_P18-011_Priscila_OvBov06-Fol_CHCA_Pos-proc.tar")

  #LOAD Ag ANNOTATIONS
  annotations_Ag=read.table("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/comparativa_Ag_Au/Ag.ref")

  #METHOD 1
  m1_indices=removeMatrix_padded(pks_Norharmane)
  m1_masses=pks_Norharmane$mass[m1_indices]
  #METHOD 2
  m2_indices=removeMatrix_compareAu(pks_Norharmane,pks_Au)
  m2_masses=pks_Norharmane$mass[m2_indices]
  # #METHOD 3
  # removeMatrix_kMeansTranspose(pks_Norharmane)
  # #METHOD 3.1
  # removeMatrix_kMeansTransposeCor(pks_Norharmane)
  #COMPARE RESULTS

  #PRINT REPORT
  print("METHOD 1: PADDED")
  print(m1_indices)
  print(m1_masses)
  print("METHOD 2: COMPARE AU")
  print(m2_indices)
  print(m2_masses)
  print("SIMILARITIES M1-M2")
  print(length(intersect(m1_indices,m2_indices)))
  distance_matrix=abs(outer(m1_masses,m2_masses,'-'))
  print(sort(distance_matrix))
  print(distance_matrix)
}

#' Benchmark.
#'
#' Compare the results of all methods to assess consistency.
#'
#' @return None
#'
#'
#' @export
benchmark <- function () {
  #LOAD DATA
  pks_Ag <- rMSIproc::LoadPeakMatrix("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/comparativa_Ag_Au/postcode3/20160801-BrainCPF-Ag-mergeddata-peaks.zip")
  pks_Au <- rMSIproc::LoadPeakMatrix("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/comparativa_Ag_Au/postcode3/20160801-BrainCPF-Au-mergeddata-peaks.zip")

  #SELECT INPUT
  pks=pks_Au
  chemical_formula="Au1"

  #BENCHMARK METHOD
  cor_thresholds=seq(0.6,1,length=100)
  scores_list=list()
  gt=generate_gt(chemical_formula,pks)#annotations_Ag[[2]]
  for(ct in cor_thresholds)
  {
    results=removeMatrix_padded_new(pks,matrix_cor = ct, exo_method = 1, max_exo=5)
    pos=pks$mass[results$pos]
    neg=pks$mass[results$neg]
    scores=compute_scores(gt,pos,neg)
    for(a in attributes(scores)$names)
    {
      scores_list[[a]]=append(scores_list[[a]],scores[[a]])
    }
    scores_list$u=append(scores_list$u,length(results$unknown))
  }

  #PLOT
  scores_list[["cor_thresholds"]]=cor_thresholds
  scores_list=data.frame(scores_list)
  scores_block1=melt(scores_list, id.vars="cor_thresholds",measure.vars=c("f1","b"))
  scores_block2=melt(scores_list, id.vars="cor_thresholds",measure.vars=c("p","r"))
  scores_block3=melt(scores_list, id.vars="cor_thresholds",measure.vars=c("u","tp","fn","fp","tn"))

  value=NULL
  variable=NULL
  plot1= ggplot(scores_block1, aes(cor_thresholds,value,color=variable) ) + geom_line()
  plot2=ggplot(scores_block2, aes(cor_thresholds,value,color=variable) ) + geom_line()
  plot3=ggplot(scores_block3, aes(cor_thresholds,value,fill=variable) ) + geom_area(stat = "identity")

  print(plot2)
  print(grid.arrange(plot1,plot3))

}
