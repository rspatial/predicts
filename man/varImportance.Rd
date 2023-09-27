\name{varImportance}

\alias{varImportance}

\title{
Get variable importance
}

\description{
Get variable importance.
}


\usage{
varImportance(model, data, vars=colnames(data), n=10)
}
 
\arguments{
  \item{model}{a model object}
  \item{data}{data.frame with data for all model variables}
  \item{vars}{character. The variables of interest}
  \item{n}{positive integer. Number of simulations}
}

\value{
named numeric vector
}

