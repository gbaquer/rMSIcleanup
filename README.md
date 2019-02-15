# rMSIcleanup
Package for annotating and removing the matrix peaks in MALDI data. It is compatible with the imzML and rMSI .tar formats.

### Installation
The simplest way to install rMSI and keep it updated is using devtools package. Install devtools from CRAN into your R session:
```R
>  install.packages("devtools")
```
Then simply tell devtools to install rMSI from github latest release:
```R
> devtools::install_github("prafols/rMSI", ref = "0.7")
```
This will install rMSI package and all of its dependencies in your R environment. Then you can access to its functions by loading the rMSI package or through the `::` operator.

### Structure
1.A paragraph that describes the high-level purpose of the package.

2.An example that shows how to use the package to solve a simple problem.

3.Installation instructions, giving code that can be copied and pasted into R.

4.An overview that describes the main components of the package. For more complex packages, this will point to vignettes for more details.