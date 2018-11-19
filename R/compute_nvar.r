
#' @author Simon Frey
#' @description Derive the parameter nvar for the hydrological model COSERO from roughness parameters of a DEM
#' @title Compute the nvar parameter from a DEM
#' @param SP a spatialPoints object at which points the values of the raster should be extracted or a character string pointing at a ESRI shapefile
#' @param DEM a raster object containing the DEM information or a character string pointing at such a raster
#' @param NZ a raster object containing the information of the NZ values of the COSERO model or a character string pointing at such a raster
#'
#' @return Returns a matrix of nvar values at specific points
#' @export
#' @import raster
#' @import rgdal
#' @import MASS
#' @import TigR

compute_nvar <- function(SP, DEM, NZ){

  library(raster)
  library(rgdal)
  library(MASS)
  library(TigR)

  if(class(DEM) == "character"){
    DEM <- raster::raster(DEM)
  } else if(class(DEM) != "RasterLayer"){
    stop("DEM must be a raster object or a character string pointing towards a raster")
  }

  if(class(NZ) == "character"){
    NZ <- raster::raster(NZ)
  } else if(class(NZ) != "RasterLayer"){
    stop("NZ must be a raster object or a character string pointing towards a raster")
  }


  if(class(SP) == "character"){
    pathparts <- strsplit(SP, "/", fixed = TRUE)[[1]]
    dsn <- paste(pathparts[1:(length(pathparts)-1)], collapse = "/")
    shp <- substr(pathparts[length(pathparts)], 1, nchar(pathparts[length(pathparts)])-4)
    SP <- readOGR(dsn = dsn, layer = shp)
  } else if (class(SP) != "SpatialPointsDataFrame"){
    stop("SP must be a SpatialPointsDataFrame object or a character string pointing towards a shapefile")
  }

  rough <- terrain(DEM,opt="roughness")

  rough1 <- (rough <= quantile(rough, 0.2)) * rough
  rough2 <- (rough <= quantile(rough, 0.4) & rough > quantile(rough, 0.2)) * rough
  rough3 <- (rough <= quantile(rough, 0.6) & rough > quantile(rough, 0.4)) * rough
  rough4 <- (rough <= quantile(rough, 0.8) & rough > quantile(rough, 0.6)) * rough
  rough5 <- (rough <= quantile(rough, 1) & rough > quantile(rough, 0.8)) * rough

  # resampling
  rough1 <- resample(rough1, NZ)
  rough2 <- resample(rough2, NZ)
  rough3 <- resample(rough3, NZ)
  rough4 <- resample(rough4, NZ)
  rough5 <- resample(rough5, NZ)

  NZMAT <- as.matrix(NZ)
  RMAT <- list(
    as.matrix(rough1),
    as.matrix(rough2),
    as.matrix(rough3),
    as.matrix(rough4),
    as.matrix(rough5)
  )

  # Über das Raster (in Form der Matrizen) loopen und die Zellen auswerten, die NZ!=NA,
  # das heißt im Gebiet liegen.
  NVAR <- matrix(data=NA,nrow=nrow(NZMAT),ncol=ncol(NZMAT))
  LLIKE <- NVAR # Log-Likelihood
  for(i in 1:nrow(NVAR)){
    for(k in 1:ncol(NVAR)){
      if(!is.na(NZMAT[i,k])){
        x <- c(RMAT[[1]][i,k],RMAT[[2]][i,k],RMAT[[3]][i,k],RMAT[[4]][i,k],RMAT[[5]][i,k])
        x = 1/x
        if(any(x %in% c(0, Inf, -Inf))){
          for(j in 1:5){
            if(x[j] %in% c(0, -Inf)) x[j] <- 0.001
            if(x[j] == Inf){
              if(j > 1){
                x[j] <- (x[j-1] + 1)^2
              } else {
                x[j] <- 0.001
              }
            }
          }
        }
        x[x <= 0] <- 0.01 # Da kein Log von 0 gezogen werden kann
        x[x == Inf] <- 500
        x[x == -Inf] <- 0.01
        fit <- MASS::fitdistr(x,"lognormal")
        NVAR[i,k] <- fit$estimate[2]
        LLIKE[i,k] <- fit$loglik
      }
    }
  }

  # Umwandeln der Matrizen in vectoren
  NZvec <- as.vector(NZMAT)
  NVARvec <- as.vector(NVAR)
  LLIKEvec <- as.vector(LLIKE)

  NVARvec <- (NVARvec-max(NVARvec,na.rm=TRUE))*-1

  output <- cbind(NZvec,NVARvec,LLIKEvec)

  # Normieren der Werte (0...1)
  max_out <- max(output[,2],na.rm=TRUE)
  min_out <- min(output[,2],na.rm=TRUE)

  normout <- (max_out-output[,2])/(max_out - min_out)
  test <- exp(normout)-1
  output[,2] <- normout*2.5
  output[,2] <- test
  output <- output[complete.cases(output),]


  return(output)

}

