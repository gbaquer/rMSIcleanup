#' Test Centroid Processing
#'
#' PTest whether the peaks are correct or not
#'
#'
#' @export
test_centroid_processing <- function()
{
  #Test with pks List
  pks_list<-rMSIproc::import_imzMLpeakList("/home/gbaquer/msidata/Ag Software Test 1/raw images/5s-Ag_B73-maize-root_sec-13_pos.imzML")
  pks_mass_list<-lapply(pks_list$peakList,function(x) x$mass)
  pks_mass<-unlist(lapply(pks_list$peakList,function(x) x$mass))
  pks_mass<-unique(pks_mass)

  b<-unlist(lapply(pks_list$peakList,function(x) length(x$mass)))
  plot(b)

  name_Ag<-c("Ag1_100","Ag1_93","Ag2_54","Ag2_100","Ag2_47","Ag3_36","Ag3_100","Ag3_93","Ag3_29","Ag4_19","Ag4_72","Ag4_100","Ag4_62","Ag4_14","Ag5_12","Ag5_54","Ag5_100","Ag5_93","Ag5_43","Ag6_35","Ag6_81","Ag6_100","Ag6_70","Ag6_26","Ag7_23","Ag7_65","Ag7_100","Ag7_93","Ag7_52","Ag7_16","Ag8_14","Ag8_46","Ag8_86","Ag8_100","Ag8_74","Ag8_35","Ag9_33","Ag9_72","Ag9_100","Ag9_93","Ag9_58","Ag9_23","Ag10_22","Ag10_55","Ag10_90","Ag10_100","Ag10_77","Ag10_41","Ag10_14","Ag11_15","Ag11_41","Ag11_77","Ag11_100","Ag11_93","Ag11_62","Ag11_29","Ag12_30","Ag12_62","Ag12_92","Ag12_100","Ag12_80","Ag12_46","Ag12_19")
  mass_Ag<-c(106.9051,108.90475,213.81019,215.80985,217.80951,320.71529,322.71494,324.7146,326.71426,427.62038,429.62004,431.6197,433.61936,435.61902,534.52548,536.52513,538.52479,540.52445,542.52411,643.43023,645.42989,647.42955,649.42921,651.42887,750.33532,752.33498,754.33464,756.3343,758.33396,760.33362,857.24042,859.24008,861.23974,863.2394,865.23906,867.23871,966.14517,968.14483,970.14449,972.14415,974.14381,976.14347,1073.05027,1075.04993,1077.04959,1079.04925,1081.0489,1083.04856,1085.04822,1179.95536,1181.95502,1183.95468,1185.95434,1187.954,1189.95366,1191.95332,1288.86012,1290.85978,1292.85944,1294.85909,1296.85875,1298.85841,1300.85807)

  tol=5e-6
  presence_list<-matrix(nrow=length(pks_mass_list),ncol=length(name_Ag))
  i=1
  for(pixel_mass in pks_mass_list)
  {
    rel_error=lapply(mass_Ag,function(m) min(abs((pixel_mass-m)/pixel_mass)))
    presence_list[i,]<-rel_error<tol
    i=i+1
  }
  print(apply(presence_list,2,sum))

  #Test with pks Matrix
  pks_mat<-rMSIproc::LoadPeakMatrix("/home/gbaquer/msidata/Ag Software Test 1/raw images/5s-Ag_B73-maize-root_sec-13_pos.zip")
  pks_mass_2<-pks_mat$mass

  rel_error=lapply(mass_Ag,function(m) min(abs((pks_mass_2-m)/pks_mass_2)))
  tol=25e-6
  print(name_Ag[which(rel_error<tol)])
}
