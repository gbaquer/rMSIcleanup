# rMSIcleanup

rMSIcleanup is an open-source R package to annotate matrix-related signals in MSI data. The algorithm takes into account the chemical formula and the spatial distribution to determine which ions are matrix-related. The algorithm incorporates an overlapping peak detection feature to prevent misclassification of overlapped or isobaric ions. Additionally, the package generates a visual report to transparently justify each annotation.


### 1. Installation

The simplest way to install rMSIcleanup and keep it updated is using devtools package. Install devtools from CRAN into your R session:
```R
>  install.packages("devtools")
```
Later, install rMSI and rMSIproc
```R
>  devtools::install_github("prafols/rMSI", ref = "0.8")
>  devtools::install_github("prafols/rMSIproc", ref = "0.2")
```
Then simply tell devtools to install rMSIcleanup from github latest release:
```R
> devtools::install_github("gbaquer/rMSIcleanup", ref = "0.1")
```
This will install rMSIcleanup package and all of its dependencies in your R environment. Then you can access its functions by loading the rMSIcleanup package or through the `::` operator.

### 2. Basic Usage
```R
## 2.0. Load Data
pks<-rMSIproc::LoadPeakMatrix("/home/gbaquer/msidata/Ag Software Test 1/images/20150526_TOF_Ag_KIDNEY/mergeddata-peaks.zip")
full<-rMSI::LoadMsiData("/home/gbaquer/msidata/Ag Software Test 1/images/20150526_TOF_Ag_KIDNEY/2015-05-26-KIDNEY-Ag-proc.tar")

## 2.1. Annotate Matrix
results<-rMSIcleanup::annotate_matrix(pks,"Ag1",full)

## 2.2. Generate Report
rMSIcleanup::generate_pdf_report(results,pks,full,"test",folder="/home/gbaquer/")

## 2.3. Remove Matrix
pks_clean<-rMSIcleanup::remove_matrix(pks,results)

# 2.4. Store Results
rMSIproc::StorePeakMatrix("before.zip",pks)
rMSIproc::StorePeakMatrix("after.zip",pks_clean)
```
