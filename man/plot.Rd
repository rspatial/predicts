\name{plot}

\docType{methods}

\alias{plot,envelope_model,missing-method}
\alias{plot,MaxEnt_model,ANY-method}

\title{
Plot predictor values 
}

\description{
Plot predictor values for occurrence (presence and absence) data in a model object.
}

\usage{
\S4method{plot}{envelope_model,missing}(x, a = 1, b = 2, p = 0.9, ...)

\S4method{plot}{MaxEnt_model,ANY}(x, y, ...)
}

\arguments{
\item{x}{model object}
\item{a}{name or position of the variable to plot in the x axis}
\item{b}{name or position of the variable to plot in the y axis}
\item{p}{percentile for coloring the points in the plot and delimiting the envelope rectangle}
\item{y}{missing}
\item{...}{arguments passed to \link[base]{plot}}
}


\keyword{spatial}

