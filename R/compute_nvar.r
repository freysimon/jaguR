
#' @author Simon Frey
#' @description Derive the parameter nvar for the hydrological model COSERO from roughness parameters of a DEM
#' @title Compute the nvar parameter from a DEM
#' @param SP a spatialPoints object at which points the values of the raster should be extracted or a character string pointing at a ESRI shapefile
#' @param DEM a raster object containing the DEM information or a character string pointing at such a raster
#' @param NZ a raster object containing the information of the NZ values of the COSERO model or a character string pointing at such a raster
#' @param showplot logical. Should the results be plotted (rastermap)
#' @param dev device for optional plotting of the results. See \code{\link{dev.new.file}}
#' @param file filename for optional plotting of the results.
#' @param ... additional parameters passed from other methods
#'
#' @return Returns a matrix of nvar values at specific points
#' @export

compute_nvar <- function(SP, DEM, NZ, showplot = TRUE, dev="dev", file = "RPlot%03d", ...){

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
  #output <- output[complete.cases(output),]

  # Normieren der Werte (0...1)
  max_out <- max(output[,2],na.rm=TRUE)
  min_out <- min(output[,2],na.rm=TRUE)

  normout <- (max_out-output[,2])/(max_out - min_out)
  test <- exp(normout)-1
  output[,2] <- normout*2.5
  output[,2] <- test
  output <- output[complete.cases(output),]

  if(!showplot){
    return(output)
  } else {
    er <- rasterize(SP, rough1, "NVAR")
    clrmp <- colorRampPalette(c("GreenYellow","yellow","gold","red"))
    dev.new.file(device = dev, file = file, width = 14, height = 7, units = "in")
    par(mar=c(1,1,4,6),oma=c(0,0,0,0))
    image.plot(er,col=clrmp(255),zlim=c(0,exp(1)-1),xaxt="n",yaxt="n",ylab="",xlab="",main="Verteilung von nvar")
    if(dev != "dev"){
      dev.off()
    }
    return(output)
  }

}

