\name{varImportance}

\alias{varImportance}

\title{
Get variable importance
}

\description{
Get variable importance.
}


\usage{
varImportance(model, y, x, n=10, ...)
}
 
\arguments{
  \item{model}{a model object}
  \item{y}{the response variable used to fit the \code{model}}
  \item{x}{data.frame with the predictor variables used to fit the \code{model}}
  \item{n}{positive integer. Number of simulations}
  \item{...}{model specific additional arguments passed to \code{predict}}
}

\value{
named numeric vector
}

\examples{

set.seed(1)
d <- data.frame(y=1:10, x1=runif(10), x2=runif(10))
m <- lm(y~., data=d)

predicts::varImportance(m, d[,1], d[,2:3])
}
