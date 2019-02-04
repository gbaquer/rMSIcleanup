#########################################################################
#
#     EXPORT MODULE
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

#' Export to mmass
#'
#' Exports the results to a text file which can be read by mmass for easy interpretation and validation of the results.
#'
#' @param pks_Matrix Peak matrix
#' @param matrix_annotation Vector determining whether each peak in the spectrum corresponds to the matrix or the tissue
#' @param metadata List containing metadata on the matrix annotation process
#'
#' @return None
#'
#'
#' @export
export_mmass <- function (pks_Matrix,matrix_annotation=rep(1,length(pks_Matrix$mass)),metadata=FALSE) {
  #IMPORT PYTHON MODULES
  struct=reticulate::import("struct")
  zlib=reticulate::import("zlib")
  base64=reticulate::import("base64")

  #Cluster image
  centers=3
  clus <- kmeans((pks_Matrix$intensity/pks_Matrix$normalizations$TIC), centers = centers)

  #Plot cluster image
  dev.new()
  rMSIproc::plotClusterImage(pks_Matrix,clus$cluster)

  for(i in 0:centers )
  {
    if(i==0){
      pixels=1:pks_Matrix$numPixels
      file_suffix="whole_image"
    }
    else{
      pixels=which(clus$cluster==i)
      file_suffix=paste("cluster_",i,sep="")
    }

    #Calculate mean spectra
    mean_spectra=apply(pks_Matrix$intensity[pixels,],2,mean)
    masses=lapply(pks_Matrix$mass,function(x) x)
    #Create table
    export_table=data.frame(mass=pks_Matrix$mass,intensity=mean_spectra)
    #Write txt file
    write.table(export_table, file = paste("output/complete_spectrum_",file_suffix,".txt",sep=""), sep = "\t", row.names=FALSE, col.names=FALSE)

    #CONVERT TO BINARY WITH THE PROPER UTF ENCODING
    intArray=""
    mzArray=""
    removed=0
    matrix_peaks=which(matrix_annotation>0)
    for (i in matrix_peaks) {
      tryCatch({
        tmp1=iconv(struct$pack("f",as.single(mean_spectra[i])))
        tmp2=iconv(struct$pack("f",as.single(masses[i])))
      },
      error=function(cond) {
        print("Skipped one peak")
        print(cond)
        tmp1=0
        tmp2=0
        removed=removed+1
      },
      finally = {
        intArray=paste(intArray,tmp1,sep = "")
        mzArray=paste(mzArray,tmp2,sep = "")
      }
      )
    }

    # COMPRESSION [IT FAILS PROBABLY DUE TO ENCODING ISSUES AGAIN]
    # intArray=zlib$compress(intArray)
    # mzArray=zlib$compress(mzArray)

    #BASE 64 ENCODING
    mzArray = base64$b64encode(mzArray)
    intArray = base64$b64encode(intArray)

    #WRITE mSD FILE
    xml <- XML::xmlTree("mSD_file")
    # names(xml)
    xml$addNode("mSD", close=FALSE, attrs=c(version="2.2"))
    #Description
    xml$addNode("description", close=FALSE)
    xml$addNode("title","MATRIX_ANNOTATION_REPORT")
    xml$addNode("date", attrs=c(value=Sys.time()))
    xml$addNode("operator", attrs=c(value=""))
    xml$addNode("contact", attrs=c(value="Gerard Baquer Gomez (gerard.baquer@urv.cat)"))
    xml$addNode("institution", attrs=c(value="MIL@B (URV)"))
    xml$addNode("instrument", attrs=c(value=""))
    xml$addNode("notes","THIS SPACE IS RESERVED FOR NOTES")
    xml$closeTag()

    #Spectrum
    xml$addNode("spectrum",attrs=c(points=(length(mean_spectra)-removed)), close=FALSE)
    #xml$addNode("spectrum",attrs=c(points=1), close=FALSE)
    xml$addNode("mzArray", mzArray, attrs=c(precision="32", endian="little"))
    xml$addNode("intArray",intArray, attrs=c(precision="32", endian="little"))
    xml$closeTag()

    #Peaklist
    #[THE PEAK LIST WILL PROBABLY INCLUDE THE PEAK MATRIX WHEN THE FUNCTION IS UPGRADED TO INCLUDE THE FULL SPECTRUM]

    #Annotations
    xml$addNode("annotations", close=FALSE)
    for (i in matrix_peaks) {
      xml$addNode("annotation", matrix_annotation[i], attrs=c(peakMZ=masses[i], peakIntensity=mean_spectra[i]))
    }
    xml$closeTag()

    #Save document
    xml$closeTag()
    saveXML(xml,paste("output/matrix_peaks_",file_suffix,".msd",sep=""),prefix='<?xml version="1.0" encoding="utf-8" ?>\n')
  }
}
