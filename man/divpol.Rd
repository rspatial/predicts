\name{divide_polygons}

\alias{divider}
\alias{stripper}

\title{
Divide polygons into equal area parts
}

\description{
\code{stripper} divides polygons into horizontal or vertical strips of a specified relative size.

\code{divider} divides a \code{SpatVector} of polygons into \code{n} compact and approximately equal area parts. The results are not deterministic so you should use \code{\link{set.seed}} to be able to reproduce your results. If you get a warning about non-convergence, you can increase the number of iterations used with additional argument \code{iter.max}
}

\usage{
divider(x, n, env=NULL, alpha=1, ...)
stripper(x, f=c(1/3, 2/3), vertical=TRUE)
}

\arguments{
  \item{x}{SpatVector of polygons}
  \item{n}{positive integer. The number of parts requested}
  \item{env}{SpatRaster with environmental data}
  \item{alpha}{numeric. One or two numbers that act as weights for the x and y coordinates}
  \item{...}{additional arguments such as \code{iter.max} passed on to \code{\link{kmeans}}}
  \item{f}{numeric vector of fractions. These must be > 0 and < 1, and in ascending order}
  \item{vertical}{logical. If \code{TRUE} the strips are vertical}
}


\value{
SpatVector 
}

\author{stripper was derived from a function by Barry Rowlingson}

\examples{
f <- system.file("ex/lux.shp", package="terra")
v <- aggregate(vect(f))
set.seed(33)
d1 <- divider(v, 10)
plot(d1)

d2 <- divider(v, 100)
boxplot(expanse(d2, "km"))

x <- stripper(v, seq(0.1, 0.9, 0.1))
round(expanse(x,"km"), 1)
plot(x, col=rainbow(12))

}

\keyword{ spatial }

