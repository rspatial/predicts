# Author: Robert J. Hijmans
# Date : December 2009
# Version 1.0
# Licence GPL v3


if (!isGeneric("envelope")) {
	setGeneric("envelope", function(x, ...)
		standardGeneric("envelope"))
}	


setClass("envelope_model",
	contains = "SDM",
	representation (
		min = "numeric",
		max = "numeric",
		factor = "list",
		names = "character"
	),	
	validity = function(object)	{
		return(TRUE)
	}
)


setMethod("envelope", signature(x="matrix"), 
	function(x, ...) {
		envelope(data.frame(x), ...)
	}
)

setMethod("envelope", signature(x="data.frame"), 
	function(x, ...) {
		n <- ncol(x)
		nv <- 0
		nms <- NULL
		bc <- new("envelope_model")
		datatype <- sapply(x, class)
		i <- (datatype == "factor" | datatype == "character")
		if (any(i)) {
			y <- x[,i,drop=FALSE]
			x <- x[,!i,drop=FALSE]
			y <- lapply(y, function(i) as.character(as.vector(stats::na.omit(i))))
			y <- y[sapply(y, length) > 0]
			if (length(y) > 0) { 
				y <- lapply(y, function(i) { 
					i <- as.data.frame(table(i))
					i[,2] <- i[,2]/max(i[,2])
					i
				})	
				nv <- length(y)
				nms <- names(y)
				bc@factor <- y
			}
		}
		haveNumeric = ncol(x) > 0
		if (haveNumeric) {
			x <- lapply(x, function(i) as.vector(stats::na.omit(i)))
			x <- x[sapply(x, length) > 1]
			nv <- nv + length(x)
			nms <- c(names(x), nms)
			bc@min <- sapply(x, min)
			bc@max <- sapply(x, max)
			bc@presence <- data.frame(x)
		}
		if (nv == 0) {
			stop("no variables with sufficient data")
		}
		if (nv < n) {
			warning("variables with insufficient data were removed")
		}
		bc@names <- nms
		bc
	}
)

setMethod("envelope", signature(x="SpatRaster"), 
	function(x, p, ...) {
		m <- extract(x, p, ID=FALSE)
		envelope(data.frame(m), ...)
	}
)

#bc <- envelope(d)

setMethod("show", signature(object="envelope_model"), 
	function(object) {
		utils::str(object)
	}
)






.percRank <- function(x, y, tail) {
	x <- sort(x)
	b <- apply(y, 1, FUN=function(z)sum(x<z))
	t <- apply(y, 1, FUN=function(z)sum(x==z))
	r <- (b + 0.5 * t)/length(x)
	
	if (tail=="both") {
		i <- which(r > 0.5)
		r[i] <- 1-r[i]
	} else if (tail == "high") {
		r[ r < 0.5 ] <- 0.5
		r <- 1-r
	} else { # if tail == "low"
		r[ r > 0.5 ] <- 0.5
	}
	r * 2
}


.envelope_predict <- function(object, tails, mincomp, maxcomp, xn, xf=NULL) {
	bc <- matrix(0, ncol=ncol(xn), nrow=nrow(xn))
	na <- apply(is.na(xn), 1, any)
	bc[na,] <- NA
	k <- apply(t(xn) >= mincomp, 2, all) & apply(t(xn) <= maxcomp, 2, all)
	k[na] <- FALSE
	nms <- colnames(xn)
	for (i in 1:ncol(bc)) {
		bc[k,i] <- .percRank(object@presence[[ nms[i] ]], xn[k, nms[i], drop=FALSE], tails[i])
	}
	return( apply(bc, 1, min) )
}


setMethod("predict", signature(object="envelope_model"), 
function(object, x, tails=NULL, ext=NULL, filename="", ...) {

	ln <- object@names

	if (is.null(tails) ) {
		tails <- rep("both", times=length(ln))
	} else {
		stopifnot(all(tails %in% c("low", "high", "both")))
		if (length(tails) == 1) {
			tails <- rep(tails, length(ln))
		} else if (length(tails) != length(ln)) {
			stop('length of "tails" is: ', length(tails), '.\nThis does not match the number of variables in the model which is: ', length(ln))
		}
	}
	mincomp <- object@min
	mincomp[tails=="high"] <- -Inf
	maxcomp <- object@max
	maxcomp[tails=="low"] <- Inf

	if (! all(ln %in% names(x)) ) {
		stop("missing variables in x")
	}
	
	if (! (inherits(x, "SpatRaster")) ) {
		x <- x[, ln ,drop=FALSE]
		return(.envelope_predict(object, tails, mincomp, maxcomp, x))

	} else {
		if ((length(ln) < length(x))) {
			x <- x[[ln]]
		}
		if (!is.null(ext)) {
			x <- terra::crop(x, ext)
		}
		out <- terra::rast(x, nlyr=1)
		names(out)  <- "envelope"
		ncols <- terra::ncol(out)
		terra::readStart(x)
		on.exit(terra::readStop(x))
		b <- terra::writeStart(out, filename, ...)
		for (i in 1:b$n) {
			rowvals <- terra::readValues(x, b$row[i], b$nrows[i], 1, ncol(x), TRUE, FALSE)
			p <- .envelope_predict(object, tails, mincomp, maxcomp, rowvals)
			terra::writeValues(out, p, b$row[i], b$nrows[i])
		}
		terra::writeStop(out)
		return(out)
	}
}
)


setMethod("plot", signature(x="envelope_model", y='missing'), 
	function(x, a=1, b=2, p=0.9, ocol="gray", icol="red", bcol="blue", cex=c(0.6, 0.6), ...) {
			
		myquantile <- function(x, p) {
			p <- min(1, max(0, p))
			x <- sort(as.vector(stats::na.omit(x)))
			if (p == 0) return(x[1])
			if (p == 1) return(x[length(x)])
			i = (length(x)-1) * p + 1
			ti <- trunc(i)
			below = x[ti]
			above = x[ti+1]
			below + (above-below)*(i-ti)  
		}
	
		p <- min(1,  max(0, p))
		if (p > 0.5) p <- 1 - p
		p <- p / 2
		d <- x@presence
		prd <- predict(x, d)
		i <- prd > p & prd < (1-p)
		if (is.character(a)) xlab <- a else xlab <- colnames(d)[a]
		if (is.character(b)) ylab <- b else ylab <- colnames(d)[b]
		plot(d[,a], d[,b], xlab=xlab, ylab=ylab, cex=0, type="n")
		type <- 6
		x1 <- quantile(d[,a], probs=p, type=type)	
		x2 <- quantile(d[,a], probs=1-p, type=type)	
		y1 <- quantile(d[,b], probs=p, type=type)	
		y2 <- quantile(d[,b], probs=1-p, type=type)	
#		x1 <- myquantile(x[,a], p)	
#		x2 <- myquantile(x[,a], 1-p)	
#		y1 <- myquantile(x[,b], p)	
#		y2 <- myquantile(x[,b], 1-p)
		cex <- rep_len(cex, 2)
		points(d[!i,a], d[!i,b], col=ocol, pch=3, cex=cex[1])
		points(d[i,a], d[i,b], xlab=colnames(x)[a], ylab=colnames(x)[b], col=icol, cex=cex[2])
		graphics::polygon(rbind(c(x1,y1), c(x1,y2), c(x2,y2), c(x2,y1), c(x1,y1)), border=bcol, lwd=2)
	}
)

