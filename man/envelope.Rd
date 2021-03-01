\name{envelope}

\alias{envelope}
\alias{envelope,SpatRaster-method}
\alias{envelope,matrix-method}
\alias{envelope,data.frame-method}

\alias{envelope_model-class}


\title{Fit a (climate) envelope model and make predictions}


\description{
The envelope algorithm has been extensively used for species distribution modeling under the name "bioclim model". This is the classic 'climate-envelope-model' that started what was later called species distribution modeling and ecological niche modeling. Although it generally does not perform as good as some other methods (Elith et al. 2006) and is unsuited for predicting climate change effects (Hijmans and Graham, 2006). It may be useful in certain cases, among other reasons because the algorithm is easy to understand and thus useful in teaching species distribution modeling. 

The algorithm computes the similarity of a location by comparing the values of environmental variables at any location to a percentile distribution of the values at known locations of occurrence ('training sites'). The closer to the 50th percentile (the median), the more suitable the location is. The tails of the distribution are not distinguished, that is, 10 percentile is treated as equivalent to 90 percentile. 

In this R implementation, percentile scores are between 0 and 1, but predicted values larger than 0.5 are subtracted from 1. Then, the minimum percentile score across all the environmental variables is computed (i.e. this is like Liebig's law of the minimum, except that high values can also be limiting factors). The final value is subtracted from 1 and multiplied with 2 so that the results are between 0 and 1. The reason for this transformation is that the results become more like that of other distribution modeling methods and are thus easier to interpret. The value 1 will rarely be observed as it would require a location that has the median value of the training data for all the variables considered. The value 0 is very common as it is assigned to all cells with a value of an environmental variable that is outside the percentile distribution (the range of the training data) for at least one of the variables. 

When using the \link[predicts]{predict} function you can choose to ignore one of the tails of the distribution (for example, to make low rainfall a limiting factor, but not high rainfall).
}


\usage{
envelope(x, ...)
}
 

\arguments{
  \item{x}{matrix or data.frame where each column is an environmental variable and each row an occurrence}
  \item{...}{ Additional arguments }
}


\value{
An object of class 'envelope_model'
}


\references{
Nix, H.A., 1986. A biogeographic analysis of Australian elapid snakes. In: Atlas of Elapid Snakes of Australia. (Ed.) R. Longmore, pp. 4-15. Australian Flora and Fauna Series Number 7. Australian Government Publishing Service: Canberra.

Booth, T.H., H.A. Nix, J.R. Busby and M.F. Hutchinson, 2014. BIOCLIM: the first species distribution modelling package, its early applications and relevance to most current MAXENT studies. Diversity and Distributions 20: 1-9

Elith, J., C.H. Graham, R.P. Anderson, M. Dudik, S. Ferrier, A. Guisan, R.J. Hijmans, F. Huettmann, J. Leathwick, A. Lehmann, J. Li, L.G. Lohmann, B. Loiselle, 
G. Manion, C. Moritz, M. Nakamura, Y. Nakazawa, J. McC. Overton, A.T. Peterson, S. Phillips, K. Richardson, R. Scachetti-Pereira, R. Schapire, J. Soberon, 
S. Williams, M. Wisz and N. Zimmerman, 2006. 
Novel methods improve prediction of species' distributions from occurrence data. Ecography 29: 129-151. \url{https://dx.doi.org/10.1111/j.2006.0906-7590.04596.x}

Hijmans R.J., and C.H. Graham, 2006. Testing the ability of climate envelope models to predict the effect of climate change on species distributions. Global change biology 12: 2272-2281. \url{https://dx.doi.org/10.1111/j.1365-2486.2006.01256.x}
}



\author{
Robert J. Hijmans
}

\seealso{
\code{\link[predicts]{predict}, \link[predicts]{maxent}}
}


\examples{
logo <- rast(system.file("ex/logo.tif", package="terra"))
#presence data
pts <- matrix(c(48.24, 48.24, 47.99, 52.88, 49.53, 46.18, 54.17, 69.62, 83.79, 85.34, 74.26, 83.79, 95.13, 84.57, 66.28, 41.8, 25.83, 3.94, 18.88, 17.33, 7.05, 13.65, 26.09, 28.54, 39.1, 44.57, 51.17, 56.26, 46.27, 38.16, 30.62, 21.95, 34.39, 59.66, 69.84, 73.23, 63.24, 45.89, 43.25, 28.36), ncol=2)

v <- extract(logo, pts)
bc <- envelope(v)

d <- logo[2360:2370]
predict(bc, d)

p1 <- predict(bc, logo)
p2 <- predict(bc, logo, tails=c('both', 'low', 'high'))

#or
}

\keyword{spatial}

