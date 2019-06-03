print_experiment_names <- function (matrix_formula, base_dirs=c("C:/Users/Gerard/Documents/1. Uni/1.5. PHD/images/Ag Software Test 1","/home/gbaquer/msidata/Ag Software Test 1"),
                            s1_threshold=0.80,s2_threshold=0.80, s3_threshold=0.7, similarity_method="euclidean", correlation_method="pearson", experiment_name="output",
                            MALDI_resolution=20000, tol_mode="scans",tol_ppm=200e-6,tol_scans=4,
                            mag_of_interest="intensity",normalization="None",
                            max_multi=10, add_list=NULL, sub_list=NULL, isobaric_detection=T,
                            save_results=T,generate_pdf=T,default_page_layout=NULL,include_summary=F,dataset_indices=NULL) {
  #0. Prepare directories
  base_dir=base_dirs[which(dir.exists(base_dirs))][1]

  images_dir=paste(base_dir,"/images",sep="")
  output_dir=paste(base_dir,"/output",sep="")

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
      pks=rMSIproc::LoadPeakMatrix(pks_name)

      pks_i=1
      for(name in pks$names)
      {
        full_spectrum_name=paste(subfolder,"/",unlist(strsplit(name,".",fixed=T))[1],"-proc.tar",sep="")

        if(file.exists(full_spectrum_name))
        {
          cat(paste(j,":",full_spectrum_name,"\n"))
          j=j+1
        }
        pks_i=pks_i+1
      }
    }
  }

}
