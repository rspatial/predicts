\name{threshold}

\alias{threshold}
\alias{threshold,paModelEvaluation-method}


\title{ Find a threshold }

\description{ 
Find a threshold (cut-off) to transform model predictions (probabilities, distances, or similar values) to a binary score (presence or absence). 
}

\usage{
\S4method{threshold}{paModelEvaluation}(x)
}
 
\arguments{
  \item{x}{paModelEvaluation object (see \code{\link{pa_evaluate}}}
}

\value{
data.frame with the following columns:

kappa: the threshold at which kappa is highest ("max kappa")

spec_sens: the threshold at which the sum of the sensitivity (true positive rate) and specificity (true negative rate) is highest

no_omission: the highest threshold at which there is no omission 

prevalence: modeled prevalence is closest to observed prevalence
		
equal_sens_spec: equal sensitivity and specificity
}

\author{
Robert J. Hijmans and Diego Nieto-Lugilde
}


\seealso{\code{\link{pa_evaluate}}}

\examples{
## See ?maxent for an example with real data.
# this is a contrived example:
# p has the predicted values for 50 known cases (locations)
# with presence of the phenomenon (species)
p <- rnorm(50, mean=0.7, sd=0.3)
# b has the predicted values for 50 background locations (or absence)
a <- rnorm(50, mean=0.4, sd=0.4)
e <- pa_evaluate(p=p, a=a)

threshold(e)
}

\keyword{spatial}

