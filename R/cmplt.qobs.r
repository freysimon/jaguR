#' Complete QOBS for use in COSERO
#' @description Stick together an xts time series with start and end dates that might lay outside the original time series
#' @author Simon Frey
#' @export
#' @import TigR
#' @param x an xts object
#' @param start a POSIXct object or a character string that can be coerced to such. Gives the desired start time of the time series.
#' @param end a POSIXct object or a character string that can be coerced to such. Gives the desired end time of the time series.
#' @param TZ time zone specification to be used for the conversion, if one is required. System-specific (see time zones), but "" is the current time zone, and "GMT" is UTC (Universal Time, Coordinated). Invalid values are most commonly treated as UTC, on some platforms with a warning.
#' @param file name of the outputfile or NULL. If the latter, no file will be written.
#' @param ... additional parameters passed on to \code{\link{as.POSIXct}}
#' @return an xts object
#' @details An xts object is extended at the beginning and end of the dates. Duplicated values (as they might exist because of time changes) are removed. NA values are substitued with -0.01 which is COSEROS standard for missing.
#'
#'    If file != NULL a txt file will be written using \code{\link{write.xts}} from the TigR package

cmplt.qobs <- function(x, start=NULL, end=NULL, TZ = "utc", file=NULL, ...){
  if(class(x)[1] != "xts"){
    stop("x must be an xts object")
  }

  # estimate periodicity of x
  p <- periodicity(x)$frequency

  if(is.null(start)){
    start <- first(index(x))
  }
  if(is.null(end)){
    end <- last(index(x))
  }

  if(is.character(start)){
    start <- as.POSIXct(start, tz = TZ, ...)
  }
  if(is.character(end)){
    end <- as.POSIXct(end, tz = TZ, ...)
  }

  if(class(start)[1] != "POSIXct"){
    stop("start must be a POSIXct object or a string that can be coerced to such")
  }
  if(class(end)[1] != "POSIXct"){
    stop("end must be a POSIXct object or a string that can be coerced to such")
  }

  myTS <- seq.POSIXt(from = start, to = end, by = p)
  myTS <- xts(rep(NA, length(myTS)), order.by = myTS)

  # merge x with myTS
  x <- cbind(myTS, x)
  x <- x[,-1]

  # clear duplicated entries
  x <- clear.duplicated(x)

  x[is.na(x)] <- -0.01

  if(!is.null(file)){
    write.xts(x, file = file, col.names = TRUE, format = "%Y  %m  %d  %H  %M", fmt = "%8.2f")
  }

  return(x)

}
