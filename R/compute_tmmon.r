#' Compute tmmon from a T_IZ file for COSERO
#' @author Simon Frey
#' @description Compute the parameter tmmon (mean monthly temperature for use in the Thornthwaite method) from a T_IZ file
#' @param x an xts object
#' @param filename NULL or a character string representing the filename of the outputfile
#' @param ... further arguments passed on to \code{\link{apply.monthly}}
#' @import xts
#' @export
#' @return a vector or nothing if filename != NULL
#'
compute_tmmon <- function(x, filename = NULL, ...){

  if(class(x)[1] != "xts"){
    stop("x must be an xts object")
  }

  x <- apply.monthly(x, FUN = mean, ...)

  mon <- month(index(x))

  xx <- NULL

  for(k in 1:12){
    xx <- c(xx, as.numeric(colMeans(x[which(mon == k),])))
  }

  if(!is.null(filename)){
    write.table(xx, file = filename, col.names = F, row.names = F, quote = F, sep = "\t")
  } else {
    return(xx)
  }

}
