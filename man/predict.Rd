\name{predict}

\docType{methods}

\alias{predict}
\alias{predict,envelope_model-method}
\alias{predict,hull_model-method}
\alias{predict,MaxEnt_model-method}
\alias{predict,MaxEnt_model_replicates-method}


\title{Spatial model predictions}

\description{
Make predictions with models defined in the \code{predicts} package
}


\usage{
\S4method{predict}{MaxEnt_model}(object, x, ext=NULL, args="", filename="", ...)

\S4method{predict}{envelope_model}(object, x, tails=NULL, ext=NULL, filename="", ...)

\S4method{predict}{hull_model}(object, x, ext=NULL, mask=FALSE, filename="", ...)
}

\arguments{
  \item{object}{model defined in this package (e.g. "envelope_model" and "maxent_model")}
  \item{x}{data to predict to. Either a data.frame or a SpatRaster}
  \item{tails}{character. You can use this to ignore the left or right tail of the percentile distribution for a variable. If supplied, tails should be a character vector with a length equal to the number of variables used in the model. Valid values are "both" (the default), "low" and "high". For example, if you have a variable x with an observed distribution between 10 and 20 and you are predicting the bioclim value for a value of 25, the default result would be zero (outside of all observed values); but if you use tail='low', the high (right) tail is ignored and the value returned will be 1. }
  \item{args}{Pass *prediction* arguments (options) to the maxent software. See \code{\link[predicts]{maxent}}}
  \item{ext}{\code{NULL} or a \code{\link[terra]{SpatExtent}} to limit the prediction to a sub-region of \code{x}}
  \item{mask}{logical. If \code{TRUE} areas that are \code{NA} in \code{x} are set to \code{NA} in the output}
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link[terra]{writeRaster}}} 
}


\value{
SpatRaster or vector (if \code{x} is a data.frame). 
}

\seealso{
\code{\link[terra]{predict}} function in the "terra" package for spatial predictions with glm, randomForest, etc.
}


\keyword{methods}
\keyword{spatial}
