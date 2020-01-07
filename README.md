# rMSIcleanup

rMSIcleanup is an open-source R package to annotate matrix-related signals in MSI data. The algorithm takes into account the chemical formula and the spatial distribution to determine which ions are matrix-related. The algorithm incorporates an overlapping peak detection feature to prevent misclassification of overlapped or isobaric ions. Additionally, the package generates a visual report to transparently justify each annotation.

More information can be found in the accompanying preprint https://www.biorxiv.org/content/10.1101/2019.12.20.884957v1.

Baquer, G., Semente, L., Garcia-Altares, M., Lee, Y. J., Chaurand, P., Correig, X., & Rafols, P. (2019). rMSIcleanup: An open-source tool for matrix-related peak annotation in mass spectrometry imaging and its application to silver-assisted laser desorption/ionization. bioRxiv.

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
## 2.1. Load Data
pks<-rMSIproc::LoadPeakMatrix("[Full Path to Peak Matrix (.zip)]")
full<-rMSI::LoadMsiData("[Full Path to Processed Data (.tar)]")

## 2.2. Annotate Matrix
results<-rMSIcleanup::annotate_matrix(pks,"Ag1",full)

## 2.3. Generate Report
rMSIcleanup::generate_pdf_report(results,pks,full,"test",folder="/home/gbaquer/")

## 2.4. Remove Matrix
pks_clean<-rMSIcleanup::remove_matrix(pks,results)

# 2.5. Store Results
rMSIproc::StorePeakMatrix("before.zip",pks)
rMSIproc::StorePeakMatrix("after.zip",pks_clean)
```
### 3. Downloading sample data
To easily try out the functionality of the package we provide a sample dataset available at http://dx.doi.org/10.17632/vsk68tjcgh.1

You need to download both the Peak Matrix ("Ag_Pancreas_TOF_2015_Dataset2.zip") and the Processed File ("Ag_Pancreas_TOF_2015_Dataset2_proc.tar").

Baquer, Gerard; Ràfols, Pere (2019), “rMSIcleanup Sample Data”, Mendeley Data, V2, doi: 10.17632/vsk68tjcgh.1

### 4. Processing your own imzML
rMSIcleanup uses data in the rMSIproc format. To annotate your own data you will have to process the imzML using the following command:

```R
> rMSIproc::ProcessWizard()
```
A window will appear where you can set up all the processing parameters, the input data and the output directory to store the results.

Refer to the repository of rMSIproc for further details (https://github.com/prafols/rMSIproc)
