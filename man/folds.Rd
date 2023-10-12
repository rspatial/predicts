\name{folds}

\alias{folds}

\title{
Make folds for k-fold partitioning
}

\description{
k-fold partitioning of a data set for model testing purposes. Each record in a matrix (or similar data structure) is randomly assigned to a group. Group numbers are between 1 and \code{k}. The function assures that each fold has the same size (or as close to that as possible).
}


\usage{
folds(x, k=5, by)
}
 
\arguments{
  \item{x}{a vector, matrix, data.frame, or Spatial object}
  \item{k}{number of groups}
  \item{by}{Optional argument. A vector or factor with sub-groups (e.g. species). Its length should be the same as the number of records in x}
}

\value{
a vector with group assignments
}

\author{
Robert J. Hijmans
}

\examples{

library(disdat)
train <- disPo("NSW")
## a single species
srsp1 <- subset(train, spid=="nsw01")
folds(srsp1, k = 5)

## all species
k = folds(train, k=5, by=train$spid)

## each group has the same number of records 
##(except for adjustments if the number of records 
## divided by k is not an integer) 

table(k[train$spid=="nsw01"])
}

\keyword{spatial}

