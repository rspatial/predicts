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
partialResponse(model, data, var=1, rng=NULL, nsteps=25)
partialResponse2(model, data, var, var2, var2levels, rng=NULL, nsteps=25)
}
 
\arguments{
  \item{model}{a model object}
  \item{data}{data.frame with data for all model variables}
  \item{var}{character or positive integer to identify the variable of interest in \code{data}}
  \item{var2}{character. A second variable of interest}
  \item{var2levels}{character. The levels of the second variable to consider}
  \item{rng}{optional vector of two numbers to set the range or the variable}
  \item{nsteps}{positive integer. Number of steps to consider for the variable}
}

\value{
data.frame
}

\examples{
fsp <- system.file("/ex/bradypus.csv", package="predicts")
occ <- read.csv(fsp)[,-1]
f <- system.file("ex/bio.tif", package="predicts")
preds <- rast(f)[[c(1,7,9)]]
v <- extract(preds, occ, ID=FALSE)

bc <- envelope(v)

pr <- partialResponse(bc, data=v, var="bio12", nsteps=30)
plot(pr, type="l", las=1)
}
