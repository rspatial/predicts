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
partialResponse(model, data, var, rng=NULL, nsteps=25)
partialResponse2(model, data, var, var2, var2levels, rng=NULL, nsteps=25)
}
 
\arguments{
  \item{model}{a model object}
  \item{data}{data.frame with data for all model variables}
  \item{var}{character. The variable of interest}
  \item{var2}{character. A second variable of interest}
  \item{var2levels}{character. The levels of the second variable to consider}
  \item{rng}{optional vector of two numbers to set the range or the variable}
  \item{nsteps}{postive integer. Number of steps to consider for the variable}
}

\value{
data.frame
}

