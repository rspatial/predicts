\name{varImportance}

\alias{varImportance}

\title{
Get partial response data
}

\description{
Get partial response data.
}


\usage{
varImportance(model, data, vars=colnames(data), n=10)
}
 
\arguments{
  \item{model}{a model object}
  \item{data}{data.frame with data for all model variables}
  \item{vars}{character. The variables of interest}
  \item{n}{poistive integer. Number of simulations}
}

\value{
data.frame
}

