#' Prepare P_IZ and T_IZ data for use in COSERO
#' @author Simon Frey
#' @import TigR
#' @import tools
#' @export
#' @param x an xts object, normally created by reading a T_IZ or P_IZ file using \code{\link{read.xts}}
#' @param filename character string for the output file(s)
#' @param fill logical. Should missing values be filled (see \code{\link{fill.missing}} for details)
#' @param divide NULL or "years". If "years", x will be divided into hydrological years
#' @param ... further arguments passed to \code{\link{write.xts}}
#' @description Prepare a P_IZ or T_IZ file to be processed by COSERO. Any missing values can be filled using \code{\link{nalocf}}. The output can be divided into hydrological years.
#' @details If the output is divided more than one file is written. The filename then is extended by the hydrological year (e.g. file_2019.txt for the period from 2018-10-01 to 2019-09-30)
#' @return an xts object (if divide == NULL) or a list with xts objects
#'
prep.data <- function(x, filename, fill = TRUE, divide = NULL, ...){

  library(TigR)
  library(tools)

  if(class(x)[1] != "xts"){
    stop("x must be an xts object")
  }

  if(fill){
    x <- fill.missing(x)
  }

  if(divide == "years"){

    index.x.is.hydrological <- month(first(index(x))) %in% c(9:10) & month(last(index(x))) %in% c(9:10)

    ext <- paste(".",file_ext(filename),sep="")
    filename <- gsub(ext,"",filename)

    if(index.x.is.hydrological){
      nHY <- length(split(x[.indexmon(x) %in% 0:8], f="years"))
    } else {
      nHY <- length(split(x[.indexmon(x) %in% 0:8], f="years"))-1
    }


    myList <- list()

    for(k in 1:nHY){

      oct2dec <- split(x[.indexmon(x) %in% 9:11], f="years")[[k]]
      if(index.x.is.hydrological){
        jan2sep <- split(x[.indexmon(x) %in% 0:8], f="years")[[k]]
      } else {
        jan2sep <- split(x[.indexmon(x) %in% 0:8], f="years")[[k+1]]
      }

      myList[[k]] <- rbind(oct2dec, jan2sep)

      myYear <- as.character(format(last(index(myList[[k]])), format = "%Y"))

      # write file
      write.xts(myList[[k]], file = paste(filename, "_", myYear, ext, sep=""),
                format = "%Y  %m  %d  %H  %M", fmt = "%6.1f", col.names = FALSE, ...)
    }

    x <- myList

  } else {
    write.xts(x, file = filename, format = "%Y  %m  %d  %H  %M", fmt = "%6.1f", col.names = FALSE, ... )
  }


  return(x)
}
