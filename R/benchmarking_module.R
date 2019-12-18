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

ratios_experiment <- function ()
{
  labels<-c("Ag1","Ag2","Ag3","Ag4","Ag5","Ag6","Ag7","Ag8","Ag9")
  masses<-c(106.905097,215.809849,322.714946,431.619698,538.524795,647.429547,754.334644,863.239396000001,970.144493)
  files<-c("/home/gbaquer/msidata/Ag Software Test 1/images/20150430_TOF_AuAg_Pancreas/mergeddata-peaks.zip",
          "/home/gbaquer/msidata/Ag Software Test 1/images/20150430_TOF_AuAg_Pancreas/mergeddata-peaks.zip",
          "/home/gbaquer/msidata/Ag Software Test 1/images/20150526_TOF_Ag_KIDNEY/mergeddata-peaks.zip",
          "/home/gbaquer/msidata/Ag Software Test 1/images/20151001-Brain2_CS7_Ag/mergeddata-peaks.zip",
          "/home/gbaquer/msidata/Ag Software Test 1/images/20151001-Brain3_THS_Ag/mergeddata-peaks.zip",
          "/home/gbaquer/msidata/Ag Software Test 1/images/20151001-Brain5_BIT_Ag/mergeddata-peaks.zip",
          "/home/gbaquer/msidata/Ag Software Test 1/images/20160601_TOF_AuAg_BrainCPF/mergeddata-peaks.zip",
          "/home/gbaquer/msidata/Ag Software Test 1/images/20160601_TOF_AuAg_BrainCPF/mergeddata-peaks.zip",
          "/home/gbaquer/msidata/Ag Software Test 1/images/20160801_TOF_AuAg_BrainCPF/mergeddata-peaks.zip",
          "/home/gbaquer/msidata/Ag Software Test 1/images/20160801_TOF_AuAg_BrainCPF/mergeddata-peaks.zip",
          "/home/gbaquer/msidata/Ag Software Test 1/images/141106_BookMouseBrain_Silver/mergeddata-peaks.zip",
          "/home/gbaquer/msidata/Ag Software Test 1/images/150209_Fingermark_75 um/mergeddata-peaks.zip",
          "/home/gbaquer/msidata/Ag Software Test 1/images/5s-Ag_B73-maize-root_sec-13_pos/mergeddata-peaks.zip",
          "/home/gbaquer/msidata/Ag Software Test 1/images/5s-Ag_B73-maize-root_sec-14_neg/mergeddata-peaks.zip")
  pks_i<-c(1,2,1,1,1,1,1,2,1,2,1,1,1,1)
  pks_list<-list()
  tol_ppm = 200e-6
  results = matrix(0,length(files),length(masses))
  indices = matrix(0,length(files),length(masses))
  for(i in 1:length(files)){
    pks<-rMSIproc::LoadPeakMatrix(files[i])
    pks<-get_one_peakMatrix(pks,pks_i[i])
    pks_list<-append(pks_list,list(pks))
    mean_spectrum <- apply(pks$intensity,2,mean)
    for(j in 1:length(masses)){
      rel_error=(abs((pks$mass-masses[j])/masses[j]))
      if(min(rel_error)<tol_ppm){
        indices[i,j]=which.min(rel_error)
        results[i,j]=mean_spectrum[which.min(rel_error)]
      }
    }
  }

  #Plotting
  #cr <- readRDS("~/msidata/Ag Software Test 1/cluster_ratios.rds")
  cr<-results
  cr <- cr[1:10,]
  normalized_cr<-t(apply(cr, 1, function(x)(x-min(x))/(max(x)-min(x))))
  output_dir="/home/gbaquer/msidata/Ag Software Test 1/output/RESULTS/"
  #Fig A Patterns Coloured by 1-3 4-6 7-10 (11-12 13-14)
  tiff(paste(output_dir,"FigSXA_cluster_spectra.tiff",sep=""), width = 100, height = 100, units = 'mm', res = 300)

  legend_labels <- c("Run 1 (Datasets 1-2)", "Run 2 (Dataset 3)", "Run 3 (Datasets 4-6)","Run 4 (Datasets 7-10)", "TOF 4 (Datasets 11-12)(external)", "Orbitrap (Datasets 13-14)(external)")
  legend_labels <- legend_labels[1:4]
  df <- data.frame(x=rep(round(masses,4),each=dim(cr)[1]),y=c(normalized_cr),z=factor(c("run1","run1","run2","run3","run3","run3","run4","run4","run4","run4","run5","run5","run6","run6")[1:10]))
  p_a<-ggplot(df)  + geom_vline(xintercept = df$x,linetype="dashed",alpha=0.5)+
    geom_point(aes(x, y, color=z))+
    theme_bw()+xlab("m/z")+ylab("Intensity")+
    scale_color_discrete("Experimental run", legend_labels,breaks=factor(c("run1","run2","run3","run4","run5","run6"))[1:4])+
    scale_x_continuous(breaks=round(masses,4),labels=labels)+
    scale_y_continuous(breaks = seq(0, 1, by = 0.2))+
    theme(legend.justification = c(1, 1), legend.position = c(1, 1),legend.background=element_rect(size=0.2,colour=1),panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0))+
    ggtitle("A")
  print(p_a)

  dev.off()
  #Fig B S1 between datasets
  cr <- results
  S1_matrix=matrix(0,14,14)
  for(i in 1:14)
    for(j in i:14){
      S1_matrix[i,j]=compute_s1(cr[i,],cr[j,],"euclidean")
      S1_matrix[j,i]=S1_matrix[i,j]
    }
  levelplot(S1_matrix, col.regions = viridis(100))
  S1_matrix<-S1_matrix[1:10,1:10]
  tiff(paste(output_dir,"FigSXB_cluster_S1.tiff",sep=""), width = 100, height = 100, units = 'mm', res = 300)
  w=dim(S1_matrix)[1]
  x <-rep(seq(1,w),w)
  y <-unlist(lapply(seq(1,w),function(x)rep(x,w)))
  z <-c(S1_matrix)
  df<-data.frame(x,y,z)
  p_b <- ggplot(df, aes(x, y, fill = z)) + geom_tile() +
    theme_bw() +
    scale_fill_gradientn("S1",limits = c(0,1), colours=viridis(100)) +
    scale_x_continuous(breaks = seq(1, w),expand=c(0,0))+scale_y_continuous(breaks = seq(1, w),expand=c(0,0))+
    theme(panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),panel.grid = element_blank(), plot.title = element_text(hjust = 0))+
    labs(title = "B", x = "Dataset", y = "Dataset")
  print(p_b)
  dev.off()
  #Fig C S2 of all images vs S2 of several random 10 peaks.
  mean_cor=rep(0,14)
  mean_random_cor=matrix(0,14,100)
  for(i in 1:14){
    cor_matrix=cor(pks_list[[i]]$intensity[,indices[i,]])
    mean_cor[i]=mean(cor_matrix)
    mean_random_cor[i,]=rep(0,100)
    for(j in 1:100){
      cor_matrix=cor(pks_list[[i]]$intensity[,sample(1:length(pks_list[[i]]$mass),9)])
      mean_random_cor[i,j]=mean(cor_matrix)
    }
  }

  tiff(paste(output_dir,"FigSXC_cluster_spectra.tiff",sep=""), width = 100, height = 100, units = 'mm', res = 300)
  values=c(rbind(t(mean_random_cor),mean_cor))
  #legend_labels <- c("Run 1 (Datasets 1-2)", "Run 2 (Dataset 3)", "Run 3 (Datasets 4-6)","Run 4 (Datasets 7-10)", "TOF 4 (Datasets 11-12)(external)", "Orbitrap (Datasets 13-14)(external)")
  #legend_labels <- legend_labels[1:4]
  df <- data.frame(x=factor(rep(1:14,each=101)),y=values,z=factor(c(rep("rand",100),"silver")))
  p_c<-ggplot(df,aes(x, y))  + #geom_vline(xintercept = df$x,linetype="dashed",alpha=0.5)+
    geom_point(aes(color=z))+
    geom_boxplot(outlier.shape = NA, aes(color="rand"))+
    theme_bw()+xlab("Dataset")+ylab("Mean Correlation")+
    scale_color_discrete("Experimental run", c("Ag peaks","Random peaks"),breaks=factor(c("silver","rand")))+
    #scale_x_continuous(breaks=round(masses,4),labels=labels)+
    scale_y_continuous(breaks = seq(0, 1, by = 0.2))+
    theme(legend.justification = c(0, 0), legend.position = c(0, 0),legend.background=element_rect(size=0.2,colour=1),panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0))+
    ggtitle("C")
  print(p_c)

  dev.off()

  #Fig D Example 9 images
  tiff(paste(output_dir,"FigSXD_example.tiff",sep=""), width = 100, height = 100, units = 'mm', res = 300)

  i=5
  plts=list()
  for(j in indices[i,])
  {
    plts=append(plts,list(rMSIcleanup::ggplot_peak_image(pks_list[[i]],pks_list[[i]]$intensity[,j],only_image = T)))
  }
  p_d<-grid.arrange(grobs=c(plts),nrow=3,padding = 2)
  print(p_d)
  dev.off()

  tiff(paste(output_dir,"FigSX.tiff",sep=""), width = 200, height = 200, units = 'mm', res = 300)
  p_final<-grid.arrange(grobs=c(list(p_a,p_b,p_c,p_d)),nrow=2)
  print(p_final)
  dev.off()
  return(results)
}
#' Batch experiment
#'
#' Run a batch of experiments
#'
#' @return None
#'
#'

batch_experiment <- function ()
{
  halogens=c("F1","Cl1","Br1","I1")
  synthetic=c("H1","H2","He1","N1O3","Th2","F2","B1F4")
  plant_origin=c("C27H56","C29H60","C31H64","C26H54O1","C28H58O1","C30H62O1","C26H52O2","C30H60O2")
  add_list=append(append(halogens,synthetic),plant_origin)

  run_experiment(add_list = add_list)
}
figures_compare_peaks <- function(){

  output_dir = "/home/gbaquer/msidata/LUMC/output/"

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
  #global_clusters=unique(unlist(lapply(results$data,function(x)unique(x$patterns_out$cluster))))
  dataset_names=NULL

  # results_files = list("/home/gbaquer/msidata/LUMC/output/AU_1/global_results_000.rds",
  #                      "/home/gbaquer/msidata/LUMC/output/AU_2/global_results_000.rds",
  #                      "/home/gbaquer/msidata/LUMC/output/AU_matrix_1/global_results_000.rds",
  #                      #"/home/gbaquer/msidata/LUMC/output/AU_matrix_2/global_results_000.rds",
  #                      "/home/gbaquer/msidata/LUMC/output/DETECT_DHB_1/global_results_000.rds",
  #                      "/home/gbaquer/msidata/LUMC/output/DETECT_DHB_2/global_results_000.rds",
  #                      #"/home/gbaquer/msidata/LUMC/output/DHB_matrix_1/global_results_000.rds",
  #                      "/home/gbaquer/msidata/LUMC/output/DETECT_DHB_matrix_2/global_results_000.rds")
  results_files = list("/home/gbaquer/msidata/LUMC/output/AU_matrix_1/global_results_000.rds",
                       "/home/gbaquer/msidata/LUMC/output/DETECT_DHB_matrix_2/global_results_000.rds")
  for(file in results_files)
  {
    r = readRDS(file)
    #Print calculated clusters
    if(i<=3||T)
      present=1:length(r$patterns_out$cluster)
    else
      present=r$patterns_out$present

    clusters=r$patterns_out$cluster[present]

    clusters=r$patterns_out$cluster
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
    col=rep(i>1,length(unique_clusters))
    # c<#ol[is_na]=1
    colors=append(colors,col)

    # s1[is.na(s1)]=0
    # s2[is.na(s2)]=0
    # s3[is.na(s3)]=0

    s1_scores=append(s1_scores,s1)
    s2_scores=append(s2_scores,s2)
    s3_scores=append(s3_scores,s3)

    clusters_list=append(clusters_list,unique_clusters)

    i=i+1

  }

  #Fig 2C Scatter Datasets vs S1 * S2
  # tiff(paste(output_dir,"Fig2C_Scatter_S1S2_Clusters_1.tiff",sep=""), width = 200, height = 100, units = 'mm', res = 300)

  unique_clusters=unique(clusters_list)
  scores=s1_scores*s2_scores
  scores=s1_scores

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
    scale_color_discrete("Matrix type", c("Au", "DHB"),breaks=c(F,T))+
    scale_x_continuous(breaks = seq(0, 85, by = 5))+scale_y_continuous(breaks = seq(0, 1, by = 0.2))+
    theme(legend.justification = c(1, 1), legend.position = c(1, 1),legend.background=element_rect(size=0.2,colour=1),panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0))+
    ggtitle("Scores of DHB clusters")
  print(p_c)
}
compare_peaks<- function()
{
  filenames<-list("/home/gbaquer/msidata/LUMC/images/PR_AU_Pos_1_tissue-proc.zip",
                  "/home/gbaquer/msidata/LUMC/images/PR_AU_Pos_2_tissue-proc.zip",
                  "/home/gbaquer/msidata/LUMC/images/PR_AU_Pos_1_matrix-proc.zip",
                  "/home/gbaquer/msidata/LUMC/images/PR_AU_Pos_2_matrix-proc.zip",
                  "/home/gbaquer/msidata/LUMC/images/PR_DHB_Pos_1-proc.zip",
                  "/home/gbaquer/msidata/LUMC/images/PR_DHB_Pos_2-proc.zip",
                  "/home/gbaquer/msidata/LUMC/images/PR_DHB_Pos_1_matrix-proc.zip",
                  "/home/gbaquer/msidata/LUMC/images/PR_DHB_Pos_2_matrix-proc.zip")
  filenames_full<-list("/home/gbaquer/msidata/LUMC/images/PR_AU_Pos_1_tissue-proc.tar",
                  "/home/gbaquer/msidata/LUMC/images/PR_AU_Pos_2_tissue-proc.tar",
                  "/home/gbaquer/msidata/LUMC/images/PR_AU_Pos_1_matrix-proc.tar",
                  "/home/gbaquer/msidata/LUMC/images/PR_AU_Pos_2_matrix-proc.tar",
                  "/home/gbaquer/msidata/LUMC/images/PR_DHB_Pos_1-proc.tar",
                  "/home/gbaquer/msidata/LUMC/images/PR_DHB_Pos_2-proc.tar",
                  "/home/gbaquer/msidata/LUMC/images/PR_DHB_Pos_1_matrix-proc.tar",
                  "/home/gbaquer/msidata/LUMC/images/PR_DHB_Pos_2_matrix-proc.tar")
  mass_list=list()
  for(filename in filenames){
    pks=rMSIproc::LoadPeakMatrix(filename)
    mass_list=append(mass_list,list(pks$mass))
  }
  full_mass_list=list()
  for(filename in filenames_full){
    full=rMSI::LoadMsiData(filename)
    full_mass_list=append(full_mass_list,list(full$mass))
  }
  tol_mode="scans"
  tol_ppm=5e-6
  tol_scans=4
  overlap=matrix(0,8,8)
  overlap_percentage=matrix(0,8,8)
  for(i in 1:length(filenames)){
    for(j in 1:i){
      if(i==j)
        overlap[i,j]=length(mass_list[[i]])
      else{
        if(tol_mode=="ppm")
          overlap[i,j]=sum(unlist(lapply(mass_list[[i]],function(x,v=mass_list[[j]])(min(abs(v-x))/x)<=tol_ppm)))
        else{
          for(x in mass_list[[i]]){
            if(closest_within_scans(x,vs=mass_list[j],masses = full_mass_list[[i]],tol_scans = tol_scans))
              overlap[i,j]=overlap[i,j]+1
          }
        }
      }
      overlap[j,i]=overlap[i,j]
      overlap_percentage[i,j]=overlap[i,j]/length(mass_list[[i]])
      overlap_percentage[j,i]=overlap[i,j]/length(mass_list[[j]])
    }
  }
  endogenous_masses=mass_list[[1]][which(unlist(lapply(mass_list[[1]],function(x)closest_within_scans(x,vs=mass_list[c(2,5,6)],masses = full_mass_list[[1]],tol_scans = tol_scans))))]
  endogenous_count=length(endogenous_masses)

  AU_masses=mass_list[[1]][which(unlist(lapply(mass_list[[1]],function(x)closest_within_scans(x,vs=mass_list[c(2,3,4)],masses = full_mass_list[[1]],tol_scans = tol_scans))))]
  AU_count=length(AU_masses)

  DHB_masses=mass_list[[5]][which(unlist(lapply(mass_list[[1]],function(x)closest_within_scans(x,vs=mass_list[c(6,7,8)],masses = full_mass_list[[5]],tol_scans = tol_scans))))]
  DHB_count=length(DHB_masses)

  DHB_matrix_masses=mass_list[[5]][which(unlist(lapply(mass_list[[5]],function(x)closest_within_scans(x,vs=mass_list[c(6,1,2)],masses = full_mass_list[[5]],intersect = c(T,F,F),tol_scans = tol_scans))))]
  DHB_matrix_count=length(DHB_matrix_masses)

  DHB_matrix_all_masses=mass_list[[5]][which(unlist(lapply(mass_list[[5]],function(x)closest_within_scans(x,vs=mass_list[c(6,8,1,2)],masses = full_mass_list[[5]],intersect = c(T,T,F,F),tol_scans = tol_scans))))]
  DHB_matrix_all_count=length(DHB_matrix_all_masses)

  #Correlation with peaks from clusters we are searching for
  results = readRDS("/home/gbaquer/msidata/LUMC/output/DETECT_DHB_1/global_results_000.rds")

  DHB_enviPat_masses=DHB_masses[which(unlist(lapply(DHB_masses,function(x)closest_within_scans(x,vs=list(results$patterns_out$mass),masses = full_mass_list[[5]],tol_scans = tol_scans))))]
  DHB_enviPat_count=length(DHB_enviPat_masses)

  #Image correlation plot
  coul <- viridis(100)
  a=log10(10*overlap_percentage[-c(4,7),-c(4,7)])
  w=dim(a)[1]
  x <-rep(seq(1,w),w)
  y <-unlist(lapply(seq(1,w),function(x)rep(x,w)))
  z <-c(a)
  df<-data.frame(x,y,z)
  plt_image_correl <- ggplot(df, aes(x, y, fill = z)) + geom_tile() +
    xlab("Dataset") + ylab("Dataset") + theme_bw() +
    scale_fill_gradientn("",limits = c(-1,1), colours=coul) +
    scale_x_continuous(breaks = seq(1, w),expand=c(0,0))+scale_y_continuous(breaks = seq(1, w),expand=c(0,0))+
    theme(panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),panel.grid = element_blank(), plot.title = element_text(hjust = 0))+
    ggtitle("Masses Overlap Percentage")
  plt_image_correl

}
closest_within_scans <- function(x,vs,masses,intersect=rep(T,length(vs)),tol_scans){

  j=which.min(abs(masses-x))
  result=T
  for(v in vs){
    i=which.min(abs(v-x))
    overlapped=((v[i]>=masses[max(1,j-tol_scans)])&&(v[i]<=masses[min(length(masses),j+tol_scans)]))
    if(!intersect)
      overlapped=!overlapped
    result=result&&overlapped
  }
  return(result)
}
sponges_experiment<- function ()
{
  matrix_formula = "Br1"
  add_list=NULL
  sub_list=NULL
  full_spectrum=rMSI::LoadMsiData("/home/gbaquer/msidata/180517_Marine_Sponges/180517_Marine_Sponges-proc.tar")
  pks=rMSIproc::LoadPeakMatrix("/home/gbaquer/msidata/180517_Marine_Sponges/mergeddata-peaks.zip")
  results=generate_gt(matrix_formula=matrix_formula,pks=pks,full_spectrum=full_spectrum,folder="/home/gbaquer/msidata/180517_Marine_Sponges/output/",add_list=add_list,sub_list=sub_list,max_multi=10,tol_scans = 4,isobaric_detection = T)
  results_file=generate_file_name("global_results_2",folder = "/home/gbaquer/msidata/180517_Marine_Sponges/output/",extension = ".rds")
  saveRDS(results, results_file)
}
LUMC_experiment <- function ()
{
  # formula="C7H6O4"
  # matrix_formula=NULL
  # for(m in 0:2){
  #   for(mh2 in 0:2){
  #     for(h in 0:2){
  #       for(k in 0:2){
  #         for(na in 0:2){
  #           if(!(m==0&&mh2==0)){
  #             if(m!=0)
  #               current_formula= enviPat::multiform(formula,m)
  #             else
  #               current_formula=""
  #             if(mh2!=0)
  #               current_formula= enviPat::mergeform(current_formula,enviPat::multiform(enviPat::subform(formula,"H2O1"),mh2))
  #             if(h!=0)
  #               current_formula= enviPat::subform(current_formula,enviPat::multiform("H1",h))
  #             if(k!=0)
  #               current_formula= enviPat::mergeform(current_formula,enviPat::multiform("K1",k))
  #             if(na!=0)
  #               current_formula= enviPat::mergeform(current_formula,enviPat::multiform("Na1",na))
  #
  #             matrix_formula=append(matrix_formula,current_formula)
  #           }
  #         }
  #       }
  #     }
  #   }
  # }
  contents = read.csv("/home/gbaquer/msidata/LUMC/DHB_clusters_BETTER.csv")
  matrix_formula=as.character(contents$Formula)
  # matrix_formula="C7H6O4"
  # add_list=c("H1","Na1","K1","K1Na1")
  # sub_list=c("H2O1","H1O1")
  add_list=NULL
  sub_list=NULL
  full_spectrum=rMSI::LoadMsiData("/home/gbaquer/msidata/LUMC/images/PR_DHB_Pos_2-proc.tar")
  pks=rMSIproc::LoadPeakMatrix("/home/gbaquer/msidata/LUMC/images/PR_DHB_Pos_2-proc.zip")
  results=generate_gt(matrix_formula=matrix_formula,pks=pks,full_spectrum=full_spectrum,folder="/home/gbaquer/msidata/LUMC/output/",add_list=add_list,sub_list=sub_list,max_multi=1,tol_scans = 4,isobaric_detection = T)
  results_file=generate_file_name("global_results_2",folder = "/home/gbaquer/msidata/LUMC/output/",extension = ".rds")
  saveRDS(results, results_file)
}
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
run_experiment <- function (matrix_formula="Ag1", base_dirs=c("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/Ag Software Test 1","/home/gbaquer/msidata/Ag Software Test 1"),
                            s1_threshold=0.80,s2_threshold=0.80, s3_threshold=0.7, similarity_method="euclidean", correlation_method="pearson", experiment_name="output",
                            MALDI_resolution=20000, tol_mode="scans",tol_ppm=100e-6,tol_scans=4,
                            mag_of_interest="intensity",normalization="None",
                            max_multi=10, add_list=NULL, sub_list=NULL, isobaric_detection=T,
                            save_results=T,generate_pdf=T,generate_figures=F,default_page_layout=NULL,include_summary=F,dataset_indices=NULL) {
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

        full_spectrum=NULL
        if(file.exists(full_spectrum_name))
        {
          full_spectrum=rMSI::LoadMsiData(full_spectrum_name)
        }
        pks_individual=get_one_peakMatrix(pks,pks_i)
        #[Potential improvement: Use ... instead]
        results$data[[j]]= generate_gt(matrix_formula=matrix_formula,pks=pks_individual,full_spectrum=full_spectrum,folder=experiment_dir,
                    s1_threshold=s1_threshold,s2_threshold=s2_threshold, s3_threshold=s3_threshold, similarity_method=similarity_method,correlation_method=correlation_method,
                    MALDI_resolution=MALDI_resolution, tol_mode=tol_mode,tol_ppm=tol_ppm,tol_scans=tol_scans,
                    mag_of_interest=mag_of_interest,normalization=normalization,
                    max_multi=max_multi, add_list=add_list, sub_list=sub_list, isobaric_detection=isobaric_detection,
                    generate_pdf=generate_pdf,generate_figures=generate_figures,default_page_layout=default_page_layout,include_summary=include_summary,pks_i = pks_i)
        results$meta$file_names=append(results$meta$file_names,full_spectrum_name)
        pks_i=pks_i+1
        j=j+1
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

#' Generate before and after files
#'
#' Generate before and after zip files
#'
#'
#' @return None
#'
#'
generate_before_after <- function () {
  #1. Open File
  base_dirs=c("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/Ag Software Test 1","/home/gbaquer/msidata/Ag Software Test 1")
  base_dir=base_dirs[which(dir.exists(base_dirs))][1]
  output_dir=paste(base_dir,"/output/",sep="")
  #experiment_dir=paste(generate_file_name("/output",extension = "",folder = output_dir),"",sep="")
  #Load Results
  #results_path = file.choose()
  #results = readRDS(results_path)
  results = readRDS("/home/gbaquer/msidata/Ag Software Test 1/output/output_000/global_results_000.rds")
  pks_name=tools::file_path_sans_ext(basename(results$meta$file_names))
  pks_folder=dirname(results$meta$file_names)
  pks_file=paste(pks_folder,"/mergeddata-peaks.zip",sep = "")
  #pks_i=unlist(lapply(seq_along(pks_folder),function(i) sum(pks_folder[1:i]==pks_folder[i])))
  i=1
  for(filename in unique(pks_file))
  {
    pks=rMSIproc::LoadPeakMatrix(filename)
    for(pks_i in 1:sum(pks_file==filename))
    {
      #2. Load before peak matrix
      pks_before=get_one_peakMatrix(pks,pks_i)
      #3. Filter out matrix-related peaks
      pks_after=get_columns_peakMatrix(pks_before,which(!results$data[[i]]$gt))
      #4. Store results
      rMSIproc::StorePeakMatrix(paste(output_dir,paste(pks_name[i],pks_i,"before.zip",sep="_"),sep=""),pks_before)
      rMSIproc::StorePeakMatrix(paste(output_dir,paste(pks_name[i],pks_i,"after.zip",sep="_"),sep=""),pks_after)
      i=i+1
    }
  }
}

#' Generate before and after files
#'
#' Generate before and after zip files
#'
#'
#' @return None
#'
#'

tsne_before_after <- function (a=1,max_iter=5000,initial_dims=100, pca=TRUE, perplexity=20, eta=100, theta=0.5) {
  #1. Open File
  base_dirs=c("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/Ag Software Test 1","/home/gbaquer/msidata/Ag Software Test 1")
  base_dir=base_dirs[which(dir.exists(base_dirs))][1]
  output_dir=paste(base_dir,"/output/",sep="")
  #experiment_dir=paste(generate_file_name("/output",extension = "",folder = output_dir),"",sep="")
  #Load Results
  #results_path = file.choose()
  #results = readRDS(results_path)
  results = readRDS("/home/gbaquer/msidata/Ag Software Test 1/output/output_000/global_results_000.rds")
  pks_name=tools::file_path_sans_ext(basename(results$meta$file_names))
  pks_folder=dirname(results$meta$file_names)
  pks_file=paste(pks_folder,"/mergeddata-peaks.zip",sep = "")
  #pks_i=unlist(lapply(seq_along(pks_folder),function(i) sum(pks_folder[1:i]==pks_folder[i])))
  i=1
  pks_file=pks_file[a]
  for(filename in unique(pks_file))
  {
    pks=rMSIproc::LoadPeakMatrix(filename)
    for(pks_i in 1:sum(pks_file==filename))
    {
      #2. Load before peak matrix
      pks_before=get_one_peakMatrix(pks,pks_i)
      #3. Filter out matrix-related peaks
      ions_after=which(!results$data[[i]]$gt)
      pks_after=get_columns_peakMatrix(pks_before,which(!results$data[[i]]$gt))

      #4. Plot results
      tsne_model_before = Rtsne::Rtsne(as.matrix(((pks_before$intensity))), check_duplicates=FALSE,max_iter=max_iter,initial_dims=initial_dims, pca=pca, perplexity=perplexity, eta=eta, theta=theta, dims=3,num_threads=10,verbose=T)

      ## getting the two dimension matrix
      d_tsne_before = as.data.frame(tsne_model_before$Y)

      color_before=results$data[[i]]$s1_scores+results$data[[i]]$s2_scores
      color_before[which(is.na(color_before))]=0

      ## plotting the results without clustering
      p_before <- ggplot(d_tsne_before, aes(x=V1, y=V2, colour=(1:length(V1))/length(V1),shape=results$data[[i]]$gt)) +
        geom_point(size=1) +
        guides(colour=guide_legend(override.aes=list(size=6))) +
        xlab("") + ylab("") +
        ggtitle(paste(i,"BEFORE")) +
        theme_light(base_size=20) +
        theme(axis.text.x=element_blank(),
              axis.text.y=element_blank()) +
        scale_color_gradientn(colours = rainbow(5)) +
        theme(legend.position="none")
        #scale_colour_gradient(low = "black", high = "red")

      #After
      tsne_model_after = Rtsne::Rtsne(as.matrix(((pks_after$intensity))), check_duplicates=FALSE,max_iter=max_iter,initial_dims=initial_dims, pca=pca, perplexity=perplexity, eta=eta, theta=theta, dims=3,num_threads=10,verbose=T)

      ## getting the two dimension matrix
      d_tsne_after = as.data.frame(tsne_model_after$Y)

      color_after=color_before[ions_after]

      ## plotting the results without clustering
      p_after<- ggplot(d_tsne_after, aes(x=V1, y=V2, colour=(1:length(V1))/length(V1))) +
        geom_point(size=1) +
        guides(colour=guide_legend(override.aes=list(size=6))) +
        xlab("") + ylab("") +
        ggtitle(paste(i,"AFTER")) +
        theme_light(base_size=20) +
        theme(axis.text.x=element_blank(),
              axis.text.y=element_blank()) +
        scale_color_gradientn(colours = rainbow(5)) +
        theme(legend.position="none")
        #scale_colour_gradient(low = "black", high = "red")

      grid.arrange(p_before, p_after,  ncol=2)
      i=i+1
    }
  }
}

#' Before after
#'
#' Generate before after pdf with plots
#'
#'
#' @return None
#'
#'

before_after <- function (a=1,max_iter=1000,initial_dims=100, pca=TRUE, perplexity=35, eta=100, theta=0.5, centers=3, normalization="None",only_pca=T) {
  #1. Open File
  base_dirs=c("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/Ag Software Test 1","/home/gbaquer/msidata/Ag Software Test 1")
  base_dir=base_dirs[which(dir.exists(base_dirs))][1]
  output_dir=paste(base_dir,"/output/",sep="")
  #experiment_dir=paste(generate_file_name("/output",extension = "",folder = output_dir),"",sep="")
  #Load Results
  #results_path = file.choose()
  #results = readRDS(results_path)
  results = readRDS("/home/gbaquer/msidata/Ag Software Test 1/output/RESULTS no isobaric detection/global_results_000.rds")
  pks_name=tools::file_path_sans_ext(basename(results$meta$file_names))
  pks_folder=dirname(results$meta$file_names)
  pks_file=paste(pks_folder,"/mergeddata-peaks.zip",sep = "")
  #pks_i=unlist(lapply(seq_along(pks_folder),function(i) sum(pks_folder[1:i]==pks_folder[i])))

  # pdf_file=generate_file_name("tsne_pca_kmeans",folder = output_dir)
  #open pdf file
  # a4_width=8.27*2*4
  # a4_height=11.69*2*4
  # pdf(pdf_file,width=a4_height,height=a4_width)
  ordered_datasets=c(11,12,1,2,3,4,5,6,7,8,9,10,13,14)

  i=1
  for(filename in unique(pks_file))
  {
    if(file.exists(filename))
    {
      pks=rMSIproc::LoadPeakMatrix(filename)
      pks$intensity=pks$intensity/pks$normalizations$RMS
      for(pks_i in 1:sum(pks_file==filename))
      {
        print(paste(i,filename,pks_i))
        #2. Load before peak matrix
        pks_before=get_one_peakMatrix(pks,pks_i)
        if(normalization!="None")
          pks_before$intensity=pks_before$intensity/pks_before$normalizations[[normalization]]
        #3. Filter out matrix-related peaks
        ions_after=which(!results$data[[i]]$gt)
        pks_after=get_columns_peakMatrix(pks_before,which(!results$data[[i]]$gt))

        #4. PCA
        #4.1. Before
        pca_model_before<-prcomp(pks_before$intensity)
        pca_before<-as.data.frame(pca_model_before$x)

        # pca_transpose_model_before<-prcomp(t(pks_before$intensity))
        # pca_transpose_before<-as.data.frame(pca_transpose_model_before$x)

        #4.2. After
        pca_model_after<-prcomp(pks_after$intensity)
        pca_after<-as.data.frame(pca_model_after$x)

        # pca_transpose_model_after<-prcomp(t(pks_after$intensity))
        # pca_transpose_after<-as.data.frame(pca_transpose_model_after$x)

        if(!only_pca)
        {
          #5. tSNE
          #5.1. Before
          tsne_model_before = Rtsne::Rtsne(as.matrix(((pks_before$intensity))), check_duplicates=FALSE,max_iter=max_iter,initial_dims=initial_dims, pca=pca, perplexity=perplexity, eta=eta, theta=theta, dims=3,num_threads=10,verbose=T)
          tsne_before = as.data.frame(tsne_model_before$Y)

          tsne_transpose_model_before = Rtsne::Rtsne(as.matrix((t(pks_before$intensity))), check_duplicates=FALSE,max_iter=max_iter,initial_dims=initial_dims, pca=pca, perplexity=perplexity, eta=eta, theta=theta, dims=3,num_threads=10,verbose=T)
          tsne_transpose_before = as.data.frame(tsne_transpose_model_before$Y)

          #5.2. After
          tsne_model_after = Rtsne::Rtsne(as.matrix(((pks_after$intensity))), check_duplicates=FALSE,max_iter=max_iter,initial_dims=initial_dims, pca=pca, perplexity=perplexity, eta=eta, theta=theta, dims=3,num_threads=10,verbose=T)
          tsne_after = as.data.frame(tsne_model_after$Y)

          tsne_transpose_model_after = Rtsne::Rtsne(as.matrix((t(pks_after$intensity))), check_duplicates=FALSE,max_iter=max_iter,initial_dims=initial_dims, pca=pca, perplexity=perplexity, eta=eta, theta=theta, dims=3,num_threads=10,verbose=T)
          tsne_transpose_after = as.data.frame(tsne_transpose_model_after$Y)

          #6. kmeans
          #6.1. Only kmeans
          kmeans_before<-kmeans(pks_before$intensity,centers)
          kmeans_after<-kmeans(pks_after$intensity,centers)
          #6.2. PCA + kmeans
          kmeans_pca_before<-kmeans(pca_before,centers)
          kmeans_pca_after<-kmeans(pca_after,centers)
          #6.3. tSNE + kmeans
          kmeans_tsne_before<-kmeans(tsne_before,centers)
          kmeans_tsne_after<-kmeans(tsne_after,centers)
        }

        #7. Plots
        file_path=generate_file_name(paste("pca_",ordered_datasets[i],sep=""),folder = output_dir,extension=".tiff")
        tiff(file_path, width = 100, height = 50, units = 'mm', res = 300)

        rgb_channels=list(c(1,0,0),
                          c(0,1,0),
                          c(0,0,1),
                          c(1,1,1))

        plots_pca_before=lapply(rgb_channels,ggplot_rgb,pos=pks_before$pos,rgb_data=pca_before)
        plots_pca_after=lapply(rgb_channels,ggplot_rgb,pos=pks_after$pos,rgb_data=pca_after)
        plots_pca = arrangeGrob(grobs=c(plots_pca_before,plots_pca_after),nrow=2,top=paste("Dataset #",ordered_datasets[i],sep=""))

        if(only_pca)
        {
          grid.arrange(plots_pca)
        }
        else
        {
          plots_tsne_before=lapply(rgb_channels,ggplot_rgb,pos=pks_before$pos,rgb_data=tsne_before)
          plots_tsne_after=lapply(rgb_channels,ggplot_rgb,pos=pks_after$pos,rgb_data=tsne_after)
          plots_tsne = arrangeGrob(grobs=c(plots_tsne_before,plots_tsne_after),nrow=2,top="tSNE images")

          plots_kmeans_before=c(list(ggplot_clusters(pks_before$pos,kmeans_before$cluster)),list(ggplot_clusters(pks_before$pos,kmeans_pca_before$cluster)),list(ggplot_clusters(pks_before$pos,kmeans_tsne_before$cluster)))
          plots_kmeans_after=c(list(ggplot_clusters(pks_after$pos,kmeans_after$cluster)),list(ggplot_clusters(pks_after$pos,kmeans_pca_after$cluster)),list(ggplot_clusters(pks_after$pos,kmeans_tsne_after$cluster)))
          plots_kmeans = arrangeGrob(grobs=c(plots_kmeans_before,plots_kmeans_after),nrow=2,top="k-means")

          plots_images=arrangeGrob(plots_pca, plots_tsne, plots_kmeans,  nrow=3)

          plots_scatter=arrangeGrob(ggplot_scatter(pca_before),ggplot_scatter(pca_after),ggplot_scatter(pca_transpose_before),ggplot_scatter(pca_transpose_after),ggplot_scatter(tsne_before),ggplot_scatter(tsne_after),ggplot_scatter(tsne_transpose_before),ggplot_scatter(tsne_transpose_after),ncol=2)

          grid.arrange(plots_images)
          grid.arrange(plots_scatter)
        }
        dev.off()
        # [Pending PCA'S and TSNE's]

        # ## plotting the results without clustering
        # p_before <- ggplot(d_tsne_before, aes(x=V1, y=V2, colour=(1:length(V1))/length(V1),shape=results$data[[i]]$gt)) +
        #   geom_point(size=1) +
        #   guides(colour=guide_legend(override.aes=list(size=6))) +
        #   xlab("") + ylab("") +
        #   ggtitle(paste(i,"BEFORE")) +
        #   theme_light(base_size=20) +
        #   theme(axis.text.x=element_blank(),
        #         axis.text.y=element_blank()) +
        #   scale_color_gradientn(colours = rainbow(5)) +
        #   theme(legend.position="none")
        # #scale_colour_gradient(low = "black", high = "red")
        #
        # ## plotting the results without clustering
        # p_after<- ggplot(d_tsne_after, aes(x=V1, y=V2, colour=(1:length(V1))/length(V1))) +
        #   geom_point(size=1) +
        #   guides(colour=guide_legend(override.aes=list(size=6))) +
        #   xlab("") + ylab("") +
        #   ggtitle(paste(i,"AFTER")) +
        #   theme_light(base_size=20) +
        #   theme(axis.text.x=element_blank(),
        #         axis.text.y=element_blank()) +
        #   scale_color_gradientn(colours = rainbow(5)) +
        #   theme(legend.position="none")
        # #scale_colour_gradient(low = "black", high = "red")
        #
        # grid.arrange(p_before, p_after,  ncol=2)
        i=i+1
      }
    }
    else{
      i=i+1
    }
  }
}
normalize <- function(x){
  return((x-min(x))/(max(x)-min(x)))
}
ggplot_rgb <- function(pos,rgb_data,v=c(1,1,1)){
  names(rgb_data) <- paste("V",1:length(names(rgb_data)),sep="")
  return(ggplot(data=rgb_data, aes(x=pos[,1], y=pos[,2], fill=rgb(normalize(V1)*v[1],normalize(V2)*v[2],normalize(V3)*v[3]))) +
            geom_tile() + scale_fill_identity()+ coord_fixed() + scale_y_reverse()+theme_void()+
            theme(panel.background = element_rect(fill = "black")))
}
ggplot_clusters <- function(pos,clus){
  tmp=as.data.frame(clus)
  return(ggplot(data=tmp, aes(x=pos[,1], y=pos[,2], fill=clus+1)) +
           geom_tile() + scale_fill_identity()+ coord_fixed() + scale_y_reverse()+theme_void()+
           theme(panel.background = element_rect(fill = "black")))
}
ggplot_scatter <- function(data){
  names(data) <- paste("V",1:length(names(data)),sep="")
  return(ggplot(data=data, aes(x=V1, y=V2, colour=normalize(1:length(V1)))) +
           geom_point() + coord_fixed() + scale_y_reverse()+theme_void()+scale_color_gradientn(colours = rainbow(5)))
}
#' Generate PDF report
#'
#' Generate PDF report
#'
#'
#' @return None
#'
#'

generate_pdf_from_file <- function (folder="RESULTS") {
  base_dirs=c("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/Ag Software Test 1","/home/gbaquer/msidata/Ag Software Test 1")
  base_dir=base_dirs[which(dir.exists(base_dirs))][1]
  output_dir=paste(base_dir,"/output/",sep="")
  #experiment_dir=paste(generate_file_name("/output",extension = "",folder = output_dir),"",sep="")
  #Load Results
  #results_path = file.choose()
  results_path = paste("/home/gbaquer/msidata/Ag Software Test 1/output/",folder,"/global_results_000.rds",sep="")
  output_dir = paste("/home/gbaquer/msidata/Ag Software Test 1/output/",folder,"/",sep="")
  results = readRDS(results_path)
  #Print Results
  generate_pdf(results,output_dir,results_path)
}
generate_single_cluster<-function(folder="/home/gbaquer/msidata/Ag Software Test 1/images/20160601_TOF_AuAg_BrainCPF",file="2016_06_01_Brain_CPF-Ag-proc.tar",pks_i=2,cluster="Ag6")
{
  full_spectrum_name=paste(folder,file,sep="/")
  pks_name=list.files(folder,pattern = "*mergeddata-peaks*",recursive = T,full.names = T)[1]
  pks <- rMSIproc::LoadPeakMatrix(pks_name)
  full_spectrum <- rMSI::LoadMsiData(full_spectrum_name)
  pks_individual=get_one_peakMatrix(pks,pks_i)
  generate_gt(cluster,pks_individual,full_spectrum,max_multi=1,generate_pdf = T,folder="/home/gbaquer/msidata/Ag Software Test 1/output/")
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

    labels=append(labels,paste(i,"_",clusters,sep="")[which(present)])
    # labels=append(labels,paste(i,"_",clusters,sep=""))

    is_na=which(is.na(s1)|is.na(s2)|is.na(s3))
    col=is.element(clusters,truth)
    col[is_na]=1
    colors=append(colors,col[which(present)])
    # colors=append(colors,col)

    s1[is.na(s1)]=0
    s2[is.na(s2)]=0
    s3[is.na(s3)]=0

    s1_scores=append(s1_scores,s1[which(present)])
    s2_scores=append(s2_scores,s2[which(present)])
    s3_scores=append(s3_scores,s3[which(present)])

    clusters_list=append(clusters_list,clusters[which(present)])

    # s1_scores=append(s1_scores,s1)
    # s2_scores=append(s2_scores,s2)
    # s3_scores=append(s3_scores,s3)
    #
    # clusters_list=append(clusters_list,clusters)

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

  #Precision Recall AUC using PRROC
  roc1=PRROC::pr.curve(scores.class0 =s1_scores ,weights.class0 =is.element(clusters_list,truth) ,curve = T)
  roc2=PRROC::pr.curve(scores.class0 =s2_scores ,weights.class0 =is.element(clusters_list,truth) ,curve = T)
  roc3=PRROC::pr.curve(scores.class0 =s3_scores ,weights.class0 =is.element(clusters_list,truth) ,curve = T)
  roc_product=PRROC::pr.curve(scores.class0 =s1_scores*s2_scores*1 ,weights.class0 =is.element(clusters_list,truth) ,curve = T)
  roc_product_improved=PRROC::pr.curve(scores.class0 =s1_scores*s2_scores+((s1_scores>0.75)*(s2_scores>0.75)*s3_scores),weights.class0 =is.element(clusters_list,truth) ,curve = T)

  roc_sum=PRROC::pr.curve(scores.class0 =s1_scores+s2_scores+0 ,weights.class0 =is.element(clusters_list,truth) ,curve = T)
  roc_euclidean=PRROC::pr.curve(scores.class0 =sqrt(((s1_scores)^2)+((s2_scores)^2)+((0)^2)) ,weights.class0 =is.element(clusters_list,truth) ,curve = T)
  #roc_thresholds=pROC::roc(response=is.element(clusters_list,truth),predictor=s1_scores*(s2_scores>0.7)*(0.8>0.7))

  plot(roc1$curve,col=2,type="l",xlim=c(0,1),ylim=c(0,1),lwd=3,xlab="Recall",ylab="Precision",cex.lab=1.6)
  lines(roc2$curve, add=T, col=3,lwd=3)
  lines(roc3$curve, add=T, col=4,lwd=3)
  lines(roc_product$curve, add=T, col=5,lwd=3)
  #lines(roc_product_improved$curve, add=T, col=5,lwd=5)
  #lines(roc_sum$curve, add=T, col=6)
  #lines(roc_euclidean$curve, add=T, col=7)
  #lines(1-roc_thresholds$specificities,roc_thresholds$sensitivities, add=T, col=7)
  #lines(0:1,0:1)
  legend=paste("S1 (",round(roc1$auc.integral,2),"AUC )")
  legend=append(legend,paste("S2 (",round(roc2$auc.integral,2),"AUC )"))
  legend=append(legend,paste("S3 (",round(roc3$auc.integral,2),"AUC )"))
  legend=append(legend,paste("S1 * S2 (",round(roc_product$auc.integral,2),"AUC )"))
  #legend=append(legend,paste("S1 * S2 (",round(roc_product_improved$auc.integral,2),"AUC )"))
  #legend=append(legend,paste("Sum (",round(roc_sum$auc.integral,2),"AUC )"))
  #legend=append(legend,paste("Euclidean (",round(roc_euclidean$auc.integral,2),"AUC )"))
  #legend=append(legend,paste("Thresholds (",round(roc_thresholds$auc.integral,3),"AUC )"))
  legend("bottomleft",legend=legend,col=2:5,lty=1,cex = 1.2,bg="white",lwd=3)


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

  #S1 * S2
  scores=
  sorted_scores=sort(s1_scores*s2_scores,decreasing = T,index.return=T)
  scores=sorted_scores$x

  unique_clusters=unique(clusters_list)
  unique_ids=1:length(unique_clusters)
  names(unique_ids)=unique_clusters

  plot(unique_ids[clusters_list],s1_scores*s2_scores,col=colors,cex.lab=2, cex.axis=2,cex=1,xlab="S1",ylab="S2")
  #S1 vs. S2
  plot(s1_scores,s2_scores,col=colors,cex.lab=2, cex.axis=2,cex=1,xlab="S1",ylab="S2")
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
  plot(s1_scores,s3_scores,col=colors,cex.lab=2, cex.axis=2,cex=1,xlab="S1",ylab="S3")
  #text(s1_scores, s3_scores, labels=labels, cex= 0.7, pos=3,col=colors)
  #legend("topleft",legend=c("Not in peak matrix","Should not be present","Should be matrix-related"),col=1:3,pch = 1)
  #abline(v = results$meta$s1_threshold)
  #abline(h = results$meta$s3_threshold)

  plot(s1_scores,s3_scores,xlim = c(results$meta$s1_threshold,1),ylim=c(results$meta$s3_threshold,1),col=colors)
  text(s1_scores, s3_scores, labels=labels, cex= 0.7, pos=3,col=colors)
  legend("topleft",legend=c("Not in peak matrix","Should not be present","Should be matrix-related"),col=1:3,pch = 1)
  abline(v = results$meta$s1_threshold)
  abline(h = results$meta$s3_threshold)

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
  dev.off()
}
test_single_file <- function(){
  pks_path = file.choose()
  full_spectrum_path = file.choose()

  start_time <- Sys.time()

  pks <- rMSIproc::LoadPeakMatrix(pks_path)

  full_spectrum <- rMSI::LoadMsiData(full_spectrum_path)
  load_time <- Sys.time()
  generate_gt("Ag1",pks,full_spectrum,add_list=c("H1","Na1","Cl1","K1","N1O3"),generate_pdf = T,folder = "C:/Users/Gerard/Documents/1. Uni/1.5. PHD/rMSIcleanup/output")
  end_time <- Sys.time()
  print(load_time-start_time)
  print(end_time-load_time)

}
#' Generate PDF report
#'
#' Generate PDF report
#'
#'
#' @return None
#'
#'

generate_pdf_from_file <- function (folder="RESULTS") {
  base_dirs=c("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/Ag Software Test 1","/home/gbaquer/msidata/Ag Software Test 1")
  base_dir=base_dirs[which(dir.exists(base_dirs))][1]
  output_dir=paste(base_dir,"/output/",sep="")
  #experiment_dir=paste(generate_file_name("/output",extension = "",folder = output_dir),"",sep="")
  #Load Results
  #results_path = file.choose()
  results_path = paste("/home/gbaquer/msidata/Ag Software Test 1/output/",folder,"/global_results_000.rds",sep="")
  output_dir = paste("/home/gbaquer/msidata/Ag Software Test 1/output/",folder,"/",sep="")
  results = readRDS(results_path)
  #Print Results
  generate_pdf(results,output_dir,results_path)
}

generate_cluster_ratio_figure<-function()
{

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

generate_figures <- function (folder="RESULTS") {
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
              panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0))+
        ggtitle("A")
  print(p_a)
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
          panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0))+
    ggtitle("B")
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
              panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0))+
         ggtitle("C")
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
test_single_file <- function(){
  pks_path = file.choose()
  full_spectrum_path = file.choose()

  start_time <- Sys.time()

  pks <- rMSIproc::LoadPeakMatrix(pks_path)

  full_spectrum <- rMSI::LoadMsiData(full_spectrum_path)
  load_time <- Sys.time()
  generate_gt("Ag1",pks,full_spectrum,add_list=c("H1","Na1","Cl1","K1","N1O3"),generate_pdf = T,folder = "C:/Users/Gerard/Documents/1. Uni/1.5. PHD/rMSIcleanup/output")
  end_time <- Sys.time()
  print(load_time-start_time)
  print(end_time-load_time)

}

#' Test Paraffin
#'
#' Test generate_gt with Paraffin
#'
#' @return None
#'
#'

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
