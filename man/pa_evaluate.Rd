\name{pa_evaluate}

\alias{paModelEvaluation-class}

\alias{pa_evaluate}
\alias{plot,paModelEvaluation,ANY-method}


\title{Presence/absence Model evaluation }

\description{ 
Evaluation of models with presence/absence data. Given a vector of presence and a vector of absence values, confusion matrices are computed for a sequence of thresholds, and model evaluation statistics are computed for each confusion matrix / threshold. 
}

\usage{
pa_evaluate(p, a, model=NULL, x=NULL, tr, ...)
}
 
\arguments{
  \item{p}{either (1) predictions for presence points  (\code{model} is \code{NULL}); or (2) predictor values for presence points (\code{model} is not \code{NULL}, \code{x} is \code{NULL}; or locations for presence points (\code{model} and \code{x} are not \code{NULL})}   
  \item{a}{as above for absence or background points}
  \item{model}{A fitted model used to make predictions}
  \item{x}{SpatRaster used to extract predictor values from}
  \item{tr}{Optional. a vector of threshold values to use for computing the confusion matrices}
  \item{...}{Additional arguments passed on to \code{predict(model,...)}}
}

\value{
pa_ModelEvaluation object
}


\section{details}{

A pa_ModelEvaluation object has the the following slots 

\describe{
   \item{\code{presence}:}{presence values used }
   \item{\code{absence}:}{absence values used}
   \item{\code{confusion}:}{confusion matrix for each threshold}
   \item{\code{stats}:}{statistics that are not threshold dependent}
   \item{\code{tr_stats}:}{statistics that are threshold dependent}
   \item{\code{thresholds}:}{optimal thresholds to classify values into presence and absence}
}

\code{stats} has the following values
\describe{
    \item{\code{np}:}{number of presence points}
    \item{\code{na}:}{number of absence points}
    \item{\code{auc}:}{Area under the receiver operator (ROC) curve}
    \item{\code{pauc}:}{ p-value for the AUC (for the Wilcoxon test W statistic}
    \item{\code{cor}:}{ Correlation coefficient}
    \item{\code{pcor}:}{p-value for correlation coefficient }
    \item{\code{prevalence}:}{ Prevalence }
    \item{\code{ODP}:}{ Overall diagnostic power }
}

\code{tr_stats} has the following values
\describe{
    \item{\code{tresholds}:}{ vector of thresholds used to compute confusion matrices }
    \item{\code{CCR}:}{ Correct classification rate }
    \item{\code{TPR}:}{ True positive rate }
    \item{\code{TNR}:}{ True negative rate }
    \item{\code{FPR}:}{ False positive rate }
    \item{\code{FNR}:}{ False negative rate }
    \item{\code{PPP}:}{ Positive predictive power }
    \item{\code{NPP}:}{ Negative predictive power }
    \item{\code{MCR}:}{ Misclassification rate }
    \item{\code{OR}:}{ Odds-ratio }
    \item{\code{kappa}:}{ Cohen's kappa }
  }

\code{thresholds} has the following values
\describe{
	\item{\code{max_kappa}:}{the threshold at which kappa is highest}
	\item{\code{max_spec_sens}:}{the threshold at which the sum of the sensitivity (true positive rate) and specificity (true negative rate) is highest}
	\item{\code{no_omission}:}{the highest threshold at which there is no omission}
	\item{\code{prevalence}:}{modeled prevalence is closest to observed prevalence}
	\item{\code{equal_sens_spec}:}{equal sensitivity and specificity}
}
}


\references{
Fielding, A.H. and J.F. Bell, 1997. A review of methods for the assessment of prediction errors in conservation presence/absence models. Environmental Conservation 24:38-49

Liu, C., M. White & G. Newell, 2011. Measuring and comparing the accuracy of species distribution models with presence-absence data. Ecography 34: 232-243.
}


\seealso{\code{\link{cm_evaluate}}}


\examples{
set.seed(0)
# p has the predicted values for 50 known cases (locations) 
# with presence of the phenomenon (species)
p <- rnorm(50, mean=0.6, sd=0.3)
# a has the predicted values for 50 background locations (or absence)
a <- rnorm(50, mean=0.4, sd=0.4)

e <- pa_evaluate(p=p, a=a)
e

e@stats

plot(e, "ROC")
plot(e, "TPR")
plot(e, "boxplot")
plot(e, "density")

str(e)
}

\keyword{spatial}

