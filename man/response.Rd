\name{partialResponse}

\alias{partialResponse}
\alias{partialResponse2}

\title{
Get partial response data
}

\description{
Get partial response data.
}


\usage{
partialResponse(model, data, var=NULL, rng=NULL, nsteps=25, plot=TRUE, nr, nc, ...)
partialResponse2(model, data, var1, var2, var2levels, rng=NULL, nsteps=25, ...)
}
 
\arguments{
  \item{model}{a model object}
  \item{data}{data.frame with data for all model variables}
  \item{var}{character or positive integer to identify the variable(s) of interest in \code{data}. If this is \code{NULL}, the partial response is computed for all variables}
  \item{var1}{character or positive integer to identify the variable of interest in \code{data}}
  \item{var2}{character. A second variable of interest}
  \item{var2levels}{character. The levels of the second variable to consider}
  \item{rng}{optional vector of two numbers to set the range or the variable}
  \item{nsteps}{positive integer. Number of steps to consider for the variable}
  \item{plot}{logical. If \code{TRUE}, the responses are plotted}
  \item{nc}{positive integer. Optional. The number of columns to divide the plotting device in (when plotting multiple variables)}
  \item{nr}{positive integer. Optional. The number of rows to divide the plotting device in (when plotting multiple variables)}
  \item{...}{model specific additional arguments passed to \code{predict}}
}

\value{
list (invisible if \code{plot=TRUE})
}

\examples{
fsp <- system.file("/ex/bradypus.csv", package="predicts")
occ <- read.csv(fsp)[,-1]
f <- system.file("ex/bio.tif", package="predicts")
preds <- rast(f)[[c(1,7,9)]]
v <- extract(preds, occ, ID=FALSE)

bc <- envelope(v)

pr <- partialResponse(bc, data=v, var=c("bio1", "bio12"), nsteps=30)
str(pr)
}
