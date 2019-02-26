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
#' @inheritParams generate_gt
#'
#'
#' @export
run_experiment <- function (matrix_formula, base_dir="C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/Ag Software Test 1",
                            s1_threshold=0.80,s2_threshold=0.85, s3_threshold=0.7, similarity_method="euclidean",
                            MALDI_resolution=20000, tol_mode="ppm",tol_ppm=200e-6,tol_scans=4,
                            mag_of_interest="intensity",normalization="None",
                            max_multi=10, add_list=NULL, sub_list=NULL,
                            generate_pdf=T,default_page_layout=NULL) {
  #0. Prepare directories
  images_dir=paste(base_dir,"/images",sep="")
  output_dir=paste(base_dir,"/output",sep="")
  report_folder=generate_file_name("Experiment",extension = "",folder = output_dir)

  #1. Generate experiment metadata file
  #[PENDING]

  #2. Run experiment
  full_spectrum_name_list=list.files(images_dir,pattern = "*-proc.tar",recursive = T)
  pks_name_list=list.files(images_dir,pattern = "*mergeddata-peaks*",recursive = T,full.names = T)
  subfolder_list=unlist(lapply(strsplit(pks_name_list,'/'),function(x) paste(x[1:length(x)-1],collapse="/")))
  num_files=length(pks_name_list)

  for(i in 1:num_files)
  {
    pks_name=pks_name_list[i]
    subfolder=subfolder_list[i]

    print(pks_name)
    pks=rMSIproc::LoadPeakMatrix(pks_name)


    for(name in pks$names)
    {
      print(name)
      full_spectrum_name=paste(subfolder,"/",unlist(strsplit(name,".",fixed=T))[1],"-proc.tar",sep="")
      print(full_spectrum_name)

      full_spectrum=rMSI::LoadMsiData(full_spectrum_name)

      #[Potential improvement: Use ... instead]
      generate_gt(matrix_formula=matrix_formula,pks=pks,full_spectrum=full_spectrum,
                  s1_threshold=s1_threshold,s2_threshold=s2_threshold, s3_threshold=s3_threshold, similarity_method=similarity_method,
                  MALDI_resolution=MALDI_resolution, tol_mode=tol_mode,tol_ppm=tol_ppm,tol_scans=tol_scans,
                  mag_of_interest=mag_of_interest,normalization=normalization,
                  max_multi=max_multi, add_list=add_list, sub_list=sub_list,
                  generate_pdf=generate_pdf,default_page_layout=default_page_layout)
    }
  }
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
  pks_Ag <- rMSIproc::LoadPeakMatrix("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/comparativa_Ag_Au/postcode3/20160801-BrainCPF-Ag-mergeddata-peaks.zip")
  pks_Norharmane <- rMSIproc::LoadPeakMatrix("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/comparativa_matrix_au/peak_matrix_norharmane/mergeddata-peaks.zip")
  pks_Au <- rMSIproc::LoadPeakMatrix("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/comparativa_matrix_au/peak_matrix_au/mergeddata-peaks.zip")
  pks_WO3 <- rMSIproc::LoadPeakMatrix("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/kidney_auagwo3/postcode1/mergeddata-peaks.zip")
  pks_HCCA <- rMSIproc::LoadPeakMatrix("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/Lluc Images/postcode2/region7_calibrated/mergeddata-peaks.zip")
  pks_Garlic <- rMSIproc::LoadPeakMatrix("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/Garlic Alex/mergeddata-peaks.zip")

  #LOAD Full spectra
  full_spectrum_Ag <- rMSI::LoadMsiData("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/comparativa_Ag_Au/postcode3/20160801-BrainCPF-Ag-CPF-proc.tar")

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
