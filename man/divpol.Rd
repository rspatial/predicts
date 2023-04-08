\name{divide_polygons}

\alias{divide_polygons}

\title{
Divide polygons into equal area parts
}

\description{
Divide a \code{SpatVector} of polygons into \code{n} compact and approximately equal area parts.
The results are not deterministic so you should use \code{\link{set.seed}} to be able to reproduce your results.

If you get a warning about non-convergence, you can increase the number of iterations used with additional argument \code{iter.max}
}

\usage{
divide_polygons(x, n, ...)
}

\arguments{
  \item{x}{SpatVector of polygons}
  \item{n}{positive integer. The number of parts requested}
  \item{...}{additional arguments such as \code{iter.max} passed on to \code{\link{kmeans}}}
}


\value{
SpatVector 
}

\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- aggregate(vect(f))
d1 <- divide_polygons(v, 10)
d2 <- divide_polygons(v, 100)
}

\keyword{ spatial }
