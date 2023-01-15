\name{MaxEnt}

\alias{MaxEnt_model-class}
\alias{MaxEnt_model_replicates-class}

\alias{MaxEnt}

\alias{MaxEnt,missing,missing-method}
\alias{MaxEnt,SpatRaster,ANY-method}
\alias{MaxEnt,SpatRaster,SpatVector-method}
\alias{MaxEnt,data.frame,numeric-method}

\alias{MaxEnt_model-class}
\alias{MaxEnt_model_replicates-class}

\title{MaxEnt}

\description{
Build a "maxent" (Maximum Entropy) species distribution model (see references below). The function uses environmental data for locations of known presence and for a large number of 'background' locations. Environmental data can be extracted from raster files. The result is a model object that can be used to predict the suitability of other locations, for example, to predict the entire range of a species. 

Background points are sampled randomly from the cells that are not \code{NA} in the first predictor variable, unless background points are specified with argument \code{a}. 

This function uses the MaxEnt species distribution model software by Phillips, Dudik and Schapire.
}

\usage{
\S4method{MaxEnt}{SpatRaster,SpatVector}(x, p, a=NULL, removeDuplicates=TRUE, nbg=10000, ...)

\S4method{MaxEnt}{data.frame,numeric}(x, p, args=NULL, path, silent=FALSE, ...)

\S4method{MaxEnt}{missing,missing}(x, p, silent=FALSE, ...)
}

\arguments{
\item{x}{Predictors. Either a SpatRaster to extract values from for the locations in \code{y}; or a data.frame, in which case each column should be a predictor variable and each row a presence or background record}

\item{p}{If \code{x} is a SpatRaster: occurence data. This can be a data.frame, matrix, or SpatVector. If \code{p} is a data.frame or matrix it represents a set of point locations; and it must have two columns with the first being the x-coordinate (longitude) and the second the y-coordinate (latitude). 

If \code{x} is a data.frame, \code{p} should be a vector with a length equal to \code{nrow(x)} and contain 0 (background) and 1 (presence) values, to indicate which records (rows) in data.frame \code{x} are presence records, and which are background records}

\item{a}{Background points. Only used if \code{p} is not a vector and not missing}

\item{nbg}{Number of background points to use. These are sampled randomly from the cells that are not \code{NA} in the first predictor variable. Ignored if background points are specified with argument \code{a}} 


\item{args}{character. Additional argument that can be passed to MaxEnt. See the MaxEnt help for more information. The R MaxEnt function only uses the arguments relevant to model fitting. There is no point in using args='outputformat=raw' when *fitting* the model; but you can use arguments relevant for *prediction* when using the predict function. Some other arguments do not apply at all to the R implementation. An example is 'outputfiletype', because the 'predict' function has its own 'filename' argument for that}

\item{removeDuplicates}{Boolean. If \code{TRUE}, duplicate presence points (that fall in the same grid cell) are removed}

\item{path}{character. Optional argument to set where you want the MaxEnt output files to be stored. This allows you to permanently keep these files. If not supplied the MaxEnt files will be stored in a temporary file. These are the files that are shown in a browser when typing the model name or when you use "show(model)"}

\item{silent}{Boolean. If \code{TRUE} a message is printed}

\item{...}{Additional arguments}
}


\value{
An object of class 'MaxEnt_model'. Or a 'MaxEnt_model_replicates' object if you use 'replicates=' as part of the \code{args} argument. 

If the function is run without any arguments a boolean value is returned (\code{TRUE} if MaxEnt.jar was found).
}

\seealso{
\code{\link[predicts]{predict}}
}

\author{Steven Phillips and Robert J. Hijmans}

\references{

\url{https://biodiversityinformatics.amnh.org/open_source/maxent/}

Steven J. Phillips, Miroslav Dudik, Robert E. Schapire, 2004.  A maximum entropy approach to species distribution modeling. Proceedings of the Twenty-First International Conference on Machine Learning. p. 655-662.

Steven J. Phillips, Robert P. Anderson, Robert E. Schapire, 2006. Maximum entropy modeling of species geographic distributions. Ecological Modelling 190:231-259. 

Jane Elith, Steven J. Phillips, Trevor Hastie, Miroslav Dudik, Yung En Chee, Colin J. Yates, 2011. A statistical explanation of MaxEnt for ecologists. Diversity and Distributions 17:43-57. \doi{j.1472-4642.2010.00725.x}
}

\examples{

\donttest{

# test if you can use MaxEnt 
MaxEnt()

if (MaxEnt()) {

# get predictor variables
f <- system.file("ex/bio.tif", package="predicts")
preds <- rast(f)
plot(preds)

# file with presence points
occurence <- system.file("/ex/bradypus.csv", package="predicts")
occ <- read.csv(occurence)[,-1]

# witholding a 20% sample for testing 
fold <- folds(occ, k=5)
occtest <- occ[fold == 1, ]
occtrain <- occ[fold != 1, ]

# fit model
me <- MaxEnt(preds, occtrain)

# see the MaxEnt results in a browser:
me

# use "args"
me2 <- MaxEnt(preds, occtrain, factors='biome', args=c("-J", "-P"))

# plot showing importance of each variable
plot(me)

# response curves
# response(me)

# predict to entire dataset
r <- predict(me, preds) 

# with some options:
r <- predict(me, preds, args=c("outputformat=raw"))

plot(r)
points(occ)

#testing
# background sample
bg <- backgroundSample(preds, 1000)

#simplest way to use 'evaluate'
e1 <- pa_evaluate(me, p=occtest, a=bg, x=preds)

# alternative 1
# extract values
pvtest <- data.frame(extract(preds, occtest))
avtest <- data.frame(extract(preds, bg))

e2 <- pa_evaluate(me, p=pvtest, a=avtest)

# alternative 2 
# predict to testing points 
testp <- predict(me, pvtest) 
head(testp)
testa <- predict(me, avtest) 

e3 <- pa_evaluate(p=testp, a=testa)
e3
threshold(e3)

plot(e3, 'ROC')
}
}

}

\keyword{spatial}

