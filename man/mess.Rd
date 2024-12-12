\name{mess}

\alias{mess}
\alias{mess,SpatRaster-method}
\alias{mess,data.frame-method}

\title{Multivariate environmental similarity surfaces (MESS)}

\description{Compute multivariate environmental similarity surfaces (MESS), as described by Elith et al., 2010}

\usage{
\S4method{mess}{SpatRaster}(x, v, full=FALSE, filename="", ...)

\S4method{mess}{data.frame}(x, v, full=FALSE)
}

\arguments{
  \item{x}{SpatRaster or data.frame}
  \item{v}{matrix or data.frame containing the reference values; each column should correspond to one layer of the SpatRaster object. If \code{x} is a SpatRaster, it can also be a SpatVector with reference locations (points)}
  \item{full}{logical. If \code{FALSE} a SpatRaster with the MESS values is returned. If \code{TRUE}, a SpatRaster is returned with \code{n} layers corresponding to the layers of the input SpatRaster and an additional layer with the MESS values}
  \item{filename}{character. Output filename (optional)}
  \item{...}{additional arguments as for \code{\link[terra]{writeRaster}}}
}

\details{
\code{v} can be obtained for a set of points using \code{\link[terra]{extract}} .
}

\value{
SpatRaster (or data.frame) with layers (columns) corresponding to the input layers and an additional layer with the mess values (if \code{full=TRUE} and \code{nlyr(x) > 1}) or a SpatRaster (data.frame) with the MESS values (if \code{full=FALSE}).
}

\references{
Elith J., M. Kearney M., and S. Phillips, 2010. The art of modelling range-shifting species. Methods in Ecology and Evolution 1:330-342. \doi{10.1111/j.2041-210X.2010.00036.x}
}

\author{
Jean-Pierre Rossi, Robert Hijmans, Paulo van Breugel
}


\examples{

set.seed(9)
r <- rast(ncol=10, nrow=10)
r1 <- setValues(r, (1:ncell(r))/10 + rnorm(ncell(r)))
r2 <- setValues(r, (1:ncell(r))/10 + rnorm(ncell(r)))
r3 <- setValues(r, (1:ncell(r))/10 + rnorm(ncell(r)))
s <- c(r1,r2,r3)
names(s) <- c('a', 'b', 'c')
xy <- cbind(rep(c(10,30,50), 3), rep(c(10,30,50), each=3))
refpt <- extract(s, xy)

ms <- mess(s, refpt, full=TRUE)
plot(ms)

\dontrun{
filename <- paste0(system.file(package="predicts"), "/ex/bradypus.csv")
bradypus <- read.table(filename, header=TRUE, sep=',')
bradypus <- bradypus[,2:3]

predfile <- paste0(system.file(package="predicts"), "/ex/bio.tif")
predictors <- rast(predfile)
reference_points <- extract(predictors, bradypus, ID=FALSE)
mss <- mess(x=predictors, v=reference_points, full=TRUE)

breaks <- c(-500, -50, -25, -5, 0, 5, 25, 50, 100)
fcol <- colorRampPalette(c("blue", "beige", "red"))
plot(mss[[10]], breaks=breaks, col=fcol(9), plg=list(x="bottomleft"))
}

}

