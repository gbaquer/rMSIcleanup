#########################################################################
#
#     UTILS MODULE
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


#' Package Global Options
#'
#' Global options used by the package
#'
#' @section verbose_level:
#' Verbose level determining which level of information is printed or plotted
#' \itemize{
#'  \item  1: Silent
#'  \item  0: Application messages
#'  \item -1: Code Sections
#'  \item -2: Code subsections
#'  \item -3: Debug mode
#' }

#'
#' @param ... Options to change
#' @param RESET Reset the package options to the default
#' @param READ.ONLY Return only the READ.ONLY options
#' @param LOCAL Change to from global to local mode when TRUE and change from local to global mode when FALSE. In the local mode a copy of the options is created and the copy is modified.
#' @param ADD Add a new option on the fly

pkg_opt= set_opt(
  verbose_level=list(.value=1,
                     .read.only=FALSE,
                     .validate= function(x) is.numeric(x) && x%%1==0 && x>=-3 && x<=1,
                     .failed_msg = "Verbose should be an integer in the range [-3,1]",
                     .description = "Verbose level determining which level of information is printed or plotted"
  )
)

#' Cosine similarity two vectors
#'
#' Returns the cosine similarity between two vectors
#'
#' @param ix Input vectors
#'
cos_sim_vectors <- function(ix)
{
  A = X[ix[1],]
  B = X[ix[2],]
  return( sum(A*B)/sqrt(sum(A^2)*sum(B^2)) )
}

#' Cosine similarity two vectors
#'
#' Returns the cosine similarity between two vectors
#'
#' @param X Input matrix
#'
cos_sim <- function(X)
{
  n <- nrow(X)
  cmb <- expand.grid(i=1:n, j=1:n)
  C <- matrix(apply(cmb,1,cos_sim_vectors),n,n)
  return(C)
}
#' Plot text
#'
#' Plots text with white background
#'
#' @param text Text to be plotted
#' @param size Size of the text
#'
# plot_text <- function(text, size=5)
# {
#   default <- par()
#   par(mar = rep(1, 4))
#   par(oma = rep(0, 4))
#   plot(NA, xlim=c(0,1), ylim=c(0,1), bty='n',
#        xaxt='n', yaxt='n', xlab='', ylab='')
#   text(0,1,text, pos=4)
#   par(default)
# }
plot_text <- function(text,size=3)
{
  p=ggplot() + geom_label()+ annotate("text", x = 0:1, y = c(1,0), size=size, label = c(text,NA), vjust=1,hjust=0) +
    coord_fixed(sqrt(2),expand = F)+
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())
  return(p)
}

ggplot_peak_image <- function(pks,values=NA,title="",isNA=F,chosen=T,in_pksMat=T)
{
  # zplots<-matrix(NA, nrow=max(pks$pos[,1]), ncol=max(pks$pos[,2]))
  # for( i in 1:nrow(pks$pos))
  # {
  #   zplots[pks$pos[ i , 1 ], pks$pos[ i , 2 ]] <- values[i]
  # }
  #levelplot(zplots)
  if(isNA)
  {
    p=ggplot() + geom_label()+ annotate("text", x = 0, y = 0, size=5, label = "NA") +
      ggtitle(title) +
      theme(axis.line=element_blank(),axis.text.x=element_blank(),
            axis.text.y=element_blank(),axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),legend.position="none",panel.border=element_blank(),panel.grid.major=element_blank(), plot.title = element_text(size=8),
            panel.grid.minor=element_blank(),plot.background=element_blank())
  }
  else
  {
    if(in_pksMat)
      if(chosen)
        background=element_rect(fill = "#45F442FF")
      else
        background=element_rect(fill = "#F46242FF")
    else
      background=element_rect(fill = "#918A88FF")

    df=data.frame(x=pks$pos[,2],y=pks$pos[,1],z=values)
    p=ggplot(df, aes(x, y, fill = z)) + geom_raster() +
      coord_fixed(1,expand = F) +
      ggtitle(title) +
      scale_fill_gradientn(colours=append("#000000FF",rev(rainbow(200,start=5.5/6,end=4.2/6)))) +
      theme(axis.line=element_blank(),axis.text.x=element_blank(),
            axis.text.y=element_blank(),axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),legend.position="none",
            panel.background=element_rect(fill = "#000000FF"),panel.border=element_blank(),panel.grid.major=element_blank(), plot.title = element_text(size=8),
            panel.grid.minor=element_blank(),plot.background=background)
  }
  return(p)
}

#' Specify decimal
#'
#' Exports the results to a text file which can be read by mmass for easy interpretation and validation of the results.
#'
#' @param x Input double
#' @param k Precision
#'
#' @return x with k precision
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

#RANDOM CHUNCKS OF CODE
# plot_text <- function(text,size=5)
# {
#   p=ggplot() + geom_label()+ annotate("text", x = 0:1, y = c(1,0), size=size, label = c(text,NA), vjust=1,hjust=0) + scale_x_continuous(expand=c(0,10))+ scale_y_continuous(expand=c(0,sqrt(2)*10)) #+ coord_fixed(sqrt(2))#+
#     # scale_x_continuous(expand=c(0,0)) +
#     # scale_y_continuous(expand=c(0,0)) #+ theme(axis.line=element_blank(),axis.text.x=element_blank(),
#   #                                             axis.text.y=element_blank(),axis.ticks=element_blank(),
#   #                                             axis.title.x=element_blank(),
#   #                                             axis.title.y=element_blank(),legend.position="none",
#   #                                             panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
#   #                                             panel.grid.minor=element_blank(),plot.background=element_blank())
#   return(p)
# }
