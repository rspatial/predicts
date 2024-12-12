\name{varImportance}

\alias{varImportance}

\title{
Get variable importance
}

\description{
Get variable importance. The importance is expressed as the deterioration of the evaluation statistic. The statistic is computed \code{n} times for model predictions after randomizing a predictor variable and subtracting the statistic for the non-randomized data. The larger the difference, the more important the variable is. 
}


\usage{
varImportance(model, y, x, n=10, stat, value="relative", ...)
}
 
\arguments{
  \item{model}{a model object}
  \item{y}{the response variable used to fit the \code{model}. If missing, it is attempted to extract it from \code{model}. If that fails, it is computed from \code{x}. In the latter case the model would be assumed to have no error}
  \item{x}{data.frame with the predictor variables used to fit the \code{model}. If missing, it is attemted to extract it from \code{model}}
  \item{n}{positive integer. Number of simulations}
  \item{stat}{character. For models with a continuous response variable this can be one of "RMSE" (the default), "AUC", or "cor". See \code{\link{RMSE}} or \code{\link{pa_evaluate}}. For models with a categorical response variable this can be one of "overall" (overall accuracy, the default) or "kappa", see \code{\link{cm_evaluate}}}
  \item{value}{character specifying how to express the output. One of , "relative" (), "difference" (), "absolute" (no adjustments)}
  \item{...}{model specific additional arguments passed to \code{predict}}
}

\value{
named numeric vector
}

\examples{

set.seed(1)
d <- data.frame(y=1:10, x1=runif(10), x2=runif(10))
m <- lm(y~., data=d)

varImportance(m, d[,1], d[,2:3])
}
