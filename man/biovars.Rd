\name{biovars}

\docType{methods}

\alias{bcvars}
\alias{bcvars,matrix,matrix,matrix-method}
\alias{bcvars,SpatRaster,SpatRaster,SpatRaster-method}
\alias{bcvars,numeric,numeric,numeric-method}

\title{bioclimatic variables}

\description{
Function to create 'bioclimatic variables' from monthly climate data. 
}

\usage{
\S4method{bcvars}{SpatRaster,SpatRaster,SpatRaster}(prec, tmin, tmax, filename="", ...)

\S4method{bcvars}{numeric,numeric,numeric}(prec, tmin, tmax)

\S4method{bcvars}{matrix,matrix,matrix}(prec, tmin, tmax)
}

\arguments{
  \item{prec}{numeric vector (12 values), matrix (12 columns), or SpatRaster with monthly (12 layers) precipitation data}
  \item{tmin}{same as \code{prec}}
  \item{tmax}{same as \code{prec}}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link[terra]{writeRaster}}}
}

\details{
Input data is normally monthly. I.e. there should be 12 values (layers) for each variable, but the function should also work for e.g. weekly data (with some changes in the meaning of the output variables. E.g. #8 would then not be for a quarter (3 months), but for a 3 week period). 
}

\value{
Same class as input, but 19 values/variables

bio1 = Mean annual temperature

bio2 = Mean diurnal range (mean of max temp - min temp)

bio3 = Isothermality (bio2/bio7) (* 100)

bio4 = Temperature seasonality (standard deviation *100)

bio5 = Max temperature of warmest month

bio6 = Min temperature of coldest month

bio7 = Temperature annual range (bio5-bio6)

bio8 = Mean temperature of the wettest quarter

bio9 = Mean temperature of driest quarter 

bio10 = Mean temperature of warmest quarter

bio11 = Mean temperature of coldest quarter

bio12 = Total (annual) precipitation

bio13 = Precipitation of wettest month

bio14 = Precipitation of driest month

bio15 = Precipitation seasonality (coefficient of variation)

bio16 = Precipitation of wettest quarter

bio17 = Precipitation of driest quarter

bio18 = Precipitation of warmest quarter
}


\examples{
tmin <- c(10,12,14,16,18,20,22,21,19,17,15,12)
tmax <- tmin + 5
prec <- c(0,2,10,30,80,160,80,20,40,60,20,0)
bcvars(prec, tmin, tmax)

tmn <- tmx <- prc <- rast(nrow=1, ncol=1, nlyr=12)
values(tmn) <- t(matrix(c(10,12,14,16,18,20,22,21,19,17,15,12)))
tmx <- tmn + 5
values(prc) <- t(matrix(c(0,2,10,30,80,160,80,20,40,60,20,0)))
b <- bcvars(prc, tmn, tmx)
as.matrix(b)
}


\keyword{spatial}

