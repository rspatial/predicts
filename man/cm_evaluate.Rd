\name{cm_evaluate}

\alias{cm_evaluate}

\title{Model evaluation with a confusion matrix}

\description{
Get model evaluation statistics from a confusion matrix. This is useful when predicting (multiple) classes.
}

\usage{
cm_evaluate(cmat, stat="overall")
}

\arguments{
  \item{cmat}{confusion matrix. Normally created with table (see examples)}
  \item{stat}{character. Either "overall" (overall accuracy), "kappa", "class" for user and producer accuracy}
}


\value{
numeric
}


\examples{  
classes <- c("forest", "water", "urban", "agriculture")
set.seed(1)
observed <- sample(classes, 100, replace=TRUE)
predicted <- observed
i <- seq(1,100,2)
predicted[i] <- sample(classes, length(i), replace=TRUE)
conmat <- table(observed, predicted)
conmat

cm_evaluate(conmat, "kappa")
cm_evaluate(conmat, "class")
}

