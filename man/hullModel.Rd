\name{hullModel}

\alias{HullModel-class}

\alias{hullModel}
\alias{hullModel,SpatVector-method}
\alias{hullModel,matrix-method}
\alias{hullModel,data.frame-method}

\alias{geometry}
\alias{geometry,HullModel-method}

\alias{plot,HullModel,missing-method}


\title{hull model}


\description{
The hull model predicts that a species is present at sites inside the a hull that contains the training points, and is absent outside that circle.

The hull can be "convex", "circle", or "rectangle"
}


\usage{
\S4method{hullModel}{SpatVector}(p, type="convex", n=1)

\S4method{hullModel}{data.frame}(p, type="convex", crs="", n=1)

\S4method{hullModel}{matrix}(p, type="convex", crs="", n=1)
}

\arguments{
  \item{p}{point locations (presence). Two column matrix, data.frame or SpatVector}
  \item{type}{character. The type of hull. One of "convex", "circle", or "rectangle"}
  \item{crs}{character. The coordinate reference system}
  \item{n}{positive integer. The number of hulls to make}
}


\value{
"HullModel"
}


\examples{
r <- rast(system.file("ex/logo.tif", package="terra"))   
#presence data
pts <- matrix(c(17, 42, 85, 70, 19, 53, 26, 84, 84, 46, 48, 85, 4,
    95, 48, 54, 66, 74, 50, 48, 28, 73, 38, 56, 43, 29, 63, 22, 46, 45,
    7, 60, 46, 34, 14, 51, 70, 31, 39, 26), ncol=2)
train <- pts[1:12, ]
test <- pts[13:20, ]
				 
ch <- hullModel(train, crs="+proj=longlat")
predict(ch, test)

plot(r)
plot(ch, border="red", lwd=2, add=TRUE)
points(train, col="red", pch=20, cex=2)
points(test, col="black", pch=20, cex=2)

pr <- predict(ch, r)
plot(pr)
points(test, col="black", pch=20, cex=2)
points(train, col="red", pch=20, cex=2)

# to get the polygons:
p <- geometry(ch)
p
}

\keyword{spatial}

