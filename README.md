# rMSIcleanup
### This package is still under development and no stable version has been released yet

rMSIcleanup is an open-source R package to annotate matrix-related signals in MSI data. The algorithm takes into account the chemical formula and the spatial distribution to determine which ions are matrix-related. The algorithm incorporates an overlapping peak detection feature to prevent misclassification of overlapped or isobaric ions. Additionally, the package generates a visual report to transparently justify each annotation.


### 1. Installation
The simplest way to install rMSIcleanup and keep it updated is using devtools package. Install devtools from CRAN into your R session:
```R
>  install.packages("devtools")
```
Then simply tell devtools to install rMSIcleanup from github latest release:
```R
> devtools::install_github("gbaquer/rMSIcleanup")
```
This will install rMSIcleanup package and all of its dependencies in your R environment. Then you can access to its functions by loading the rMSIcleanup package or through the `::` operator.

### 2. Usage example

### 3. Overview of the package
