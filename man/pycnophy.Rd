\name{pycnophy}

\alias{pycnophy}

\title{
Pycnophylactic interpolation.
}

\description{
Given a \code{SpatVector} of polygons and population data for each polygon, compute a population density estimate based on Tobler's pycnophylactic interpolation algorithm.
}

\usage{
pycnophy(x, v, pop, r = 0.2, converge = 3, verbose=FALSE)
}

\arguments{
  \item{x}{SpatRaster to interpolate to}
  \item{v}{SpatVector of polygons}
  \item{pop}{Either a character (name in \code{v}) or a numeric vector of length \code{nrow(v)}}
  \item{r}{ A relaxation parameter for the iterative step in the pycnophylactic algorithm.  Prevents over-compensation in the smoothing step.  In practice the default value works well}
  \item{converge}{ A convergence parameter, informing the decision on when iterative improvements on the smooth surface have converged sufficiently - see details}
  \item{verbose}{If \code{TRUE} the function report the maximum change in any grid cell value for each iterative step}
}

\details{
This method uses an iterative approach, and for each iteration notes the maximum change in a pixel.  When this value falls below a certain level (10^(-converge) times the largest initial grid cell value) the iteration stops.  
}

\value{
SpatRaster 
}

\references{
Tobler, W.R. (1979) \emph{Smooth Pycnophylactic Interpolation for Geographical Regions}. Journal of the American Statistical Association, v74(367) pp. 519-530.
}

\author{
Chris Brunsdon (adapted for terra objects by Robert Hijmans)
}

\note{
Pycnophylatic interpolation has the property that the sum of the estimated values associated with all of the pixels in any polygon equals the supplied population for that polygon.  A further property is that all pixel values are greater than or equal to zero.  The method is generally used to obtain pixel-based population estimates when total populations for a set of irregular polygons (eg. counties) are known.
}

\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- vect(f)
r <- rast(v, resolution = 0.01)
p <- pycnophy(r, v, "POP", converge=3, verbose=FALSE)
plot(p); lines(v)
}

\keyword{ smoothing }
\keyword{ spatial }
