#' @author Simon Frey
#' @title Calculate NVAR
#' @description Calculate COSERO's parameter NVAR from roughness properties derived from a DEM
#' @param rough SpatRaster object containing roughness information. Might be NA. See details.
#' @param DEM SpatRaster object containing a digital elevation model where roughness information can be derived from. Might be NA. See details.
#' @param x SpatVector object with points geometry. At these locations, NVAR will be extracted.
#' @param resamp SpatRaster object or FALSE. If resamp is a SpatRaster object, rough is resampled to match the spatial resolution of resamp.
#' @imports terra
#' @imports MASS
#' @export
#' @return A two column matrix containing the ID of the points in x and the computed values of NVAR.
#' @details
#' This function replaced \\link{compute.nvar}, which is deprecated. It calculates NVAR at certain points given by the parameter x. Eigher DEM or rough must be given, the other one then need's to be NA. If the first, roughness will be calculated.
#'

calc.nvar <- function(rough=NA, DEM=NA, x, resamp = FALSE){

  # check if the parameters match

  if(is.logical(rough) && is.logical(DEM)){
    stop("Either rough or DEM must be provided.")
  }
  if(class(rough) != "SpatRaster" && is.logical(DEM)){
    stop("rough must be a SpatRaster object")
  }
  if(class(DEM) != "SpatRaster" && is.logical(rough)){
    stop("DEM must be a SpatRaster object")
  }
  if(class(rough) == "SpatRaster" && class(DEM) == "SpatRaster"){
    warning("Both rough and DEM were provided. Using rough and ignoring DEM.")
    DEM <- NA
  }


  use.resample = TRUE
  if(is.logical(resamp)){
    if(resamp){
      stop("If the roughness should be resamples, resamp must be a SpatRaster object.")
    }
    else {
      use.resample = FALSE
    }
  }

  if(class(x) != "SpatVector"){
    stop("x must be a SpatVector object.")
  }

  # computing roughness if DEM is given
  if(class(rough) != "SpatRaster"){
    rough <- terra::terrain(DEM, v = "rougness")
    DEM <- NA
  }

  # calculating quantiles of roughness to fit a lognormal distribuion to
  rough1 <- (rough <= stats::quantile(values(rough), probs=0.2, na.rm = TRUE)) * rough
  rough2 <- (rough <= stats::quantile(values(rough), probs=0.4, na.rm = TRUE) & rough > stats::quantile(values(rough), probs=0.2, na.rm = TRUE)) * rough
  rough3 <- (rough <= stats::quantile(values(rough), probs=0.6, na.rm = TRUE) & rough > stats::quantile(values(rough), probs=0.4, na.rm = TRUE)) * rough
  rough4 <- (rough <= stats::quantile(values(rough), probs=0.8, na.rm = TRUE) & rough > stats::quantile(values(rough), probs=0.6, na.rm = TRUE)) * rough
  rough5 <- (rough <= stats::quantile(values(rough), probs=1, na.rm = TRUE) & rough > stats::quantile(values(rough), probs=0.8, na.rm = TRUE)) * rough

  # resample rough1 to rough5 rasters if resamp != FALSE
  if(isTRUE(use.resample)){

    if(prod(res(DEM)) >= prod(res(resamp))){
      warning("Resolution of resamp is higher than the resolution of rough. No resampling possible. Skipping resample.")
    } else {
      rough1 <- terra::resample(rough1, resamp, method = "bilinear")
      rough2 <- terra::resample(rough2, resamp, method = "bilinear")
      rough3 <- terra::resample(rough3, resamp, method = "bilinear")
      rough4 <- terra::resample(rough4, resamp, method = "bilinear")
      rough5 <- terra::resample(rough5, resamp, method = "bilinear")
    }
  }

  # extract points from rasters
  r1 <- terra::extract(rough1, x)
  r2 <- terra::extract(rough2, x)
  r3 <- terra::extract(rough3, x)
  r4 <- terra::extract(rough4, x)
  r5 <- terra::extract(rough5, x)

  rx <- cbind(r1[,2],r2[,2],r3[,2],r4[,2],r5[,2])

  rx[rx <= 0] <- 0.001 # Because log(0) = NA
  rx[rx == Inf] <- 500
  rx[rx == -Inf] <- 0.001
  rx[is.na(rx)] <- 0.001

  # define fitting function
  fitfunc <- function(x){
    y <- MASS::fitdistr(x,"lognormal")
    return(y$estimate[2])
  }

  fit <- apply(rx, MARGIN = 1, FUN = fitfunc)

  val <- cbind(r1[,1], fit)

  colnames(val) <- c("x","NVAR")

  return(val)

}
