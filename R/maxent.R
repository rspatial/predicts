# Author: Robert J. Hijmans, r.hijmans@gmail.com
# Date: 2009-2021
# Version 0.1
# Licence GPL v3

setClass("maxent_model",
	representation (
		lambdas  = "vector",
		results = "matrix",
		path = "character",
		html = "character"
	),	
	prototype (	
		lambdas = as.vector(NA),
		results = as.matrix(NA),
		path = "",
		html = ""
	),
)



setClass("maxent_model_replicates",
	representation (
		models  = "list",
		results = "matrix",
		html = "character"
	),	
	prototype (	
		models = list(),
		results = as.matrix(NA),
		html = ""
	),
)


setMethod ("show" , "maxent_model_replicates", 
	function(object) {
		cat("class     :" , class(object), "\n")
		cat("replicates:", length(object@models), "\n")
		if (file.exists(object@html)) {
			utils::browseURL( paste("file:///", object@html, sep="") )
		} else {
			cat("model html no longer exists\n")
		}
	}
)	

		
		
		

setMethod ("show" , "maxent_model", 
	function(object) {
		cat("class    :" , class(object), "\n")
		cat("variables:", colnames(object@presence), "\n")
		# cat("lambdas\n")
		# print(object@lambdas)
#		pp <- nrow(object@presence)
#		cat("\npresence points:", pp, "\n")
#		if (pp < 5) { 
#			print(object@presence)
#		} else {
#			print(object@presence[1:5,])
#			cat("  (... ...  ...)\n")
#			cat("\n")
#		}
#		pp <- nrow(object@absence)
#		cat("\nabsence points:", pp, "\n")
#		if (pp < 5) {
#			print(object@absence)
#		} else {
#			print(object@absence[1:5,])
#			cat("  (... ...  ...)\n")
#			cat("\n")
#		}
#		cat("\nmodel fit\n")
#		print(object@results)
#		cat("\n")

		if (file.exists(object@html)) {
			utils::browseURL( paste("file:///", object@html, sep="") )
		} else {
			cat("output html file no longer exists\n")
		}
	}
)	


if (!isGeneric("maxentropy")) {
	setGeneric("maxentropy", function(x, p, ...)
		standardGeneric("maxentropy"))
}	


.getMeVersion <- function(...) {}

setMethod("maxentropy", signature(x="missing", p="missing"), 
	function(x, p, silent=FALSE, ...) {

		if (is.null(getOption("dismo_rJavaLoaded"))) {
			# to avoid trouble on macs
			Sys.setenv(NOAWT=TRUE)
			if ( requireNamespace("rJava") ) {
				rJava::.jpackage("dismo")
				options(dismo_rJavaLoaded=TRUE)
			} else {
				if (!silent) {
					cat("Cannot load rJava\n")			
				}
				return(FALSE)
			}
		}

		if (is.null(getOption("dismo_maxent"))) {
			mxe <- rJava::.jnew("meversion") 
			v <- try(rJava::.jcall(mxe, "S", "meversion"), silent=TRUE)
			if (class(v) == "try-error") {
				if (!silent) {
					cat("maxent_model is missing or incompatible with your version of Java\n")
				}
				return(FALSE)
			} else if (v == "3.3.3a") {
				if (!silent) {
					cat("This is not a compatible version of Maxent\n")
				}
				return(FALSE)
			}
			options(dismo_maxent=v)
		} else {
			v = getOption("dismo_maxent")
		}
		if (!silent) {
			cat("This is maxent_model version", v, "\n")
		}
		invisible(TRUE)
	}
)



.getMatrix <- function(x) {
	if (inherits(x, "SpatVector")) {
		x <- geom(x)[,c("x", "y")]
	} 
	if (inherits(x, "matrix")) {
		x <- data.frame(x)
	}
	if (! class(x) == "data.frame" ) {
		stop("data should be  a matrix, data.frame, or SpatVector")
	}
	if (dim(x)[2] != 2) {
		stop("presence or absence coordinates data should be a matrix or data.frame with 2 columns" ) 	
	}
	colnames(x) <- c("x", "y")
	return(x)
} 


setMethod("maxentropy", signature(x="SpatRaster", p="ANY"), 
	function(x, p, a=NULL, factors=NULL, removeDuplicates=TRUE, nbg=10000, ...) {

		p <- .getMatrix(p)
		if (removeDuplicates) {
			cells <- unique(cellFromXY(x, p))
			pv <- extract(x, cells)
		} else {
			pv <- extract(x, p)
		}
		if (is.null(dim(pv))) {
			pv <- matrix(pv, ncol=1)
			colnames(pv) <- names(x)
		}
		pv <- data.frame(pv)
		

		lpv <- nrow(pv)
		pv <- stats::na.omit(pv)
		nas <- lpv - nrow(pv)
		if (nas > 0) {
			if (nas >= 0.5 * lpv) {
				stop("more than half of the presence points have NA predictor values")
			} else {
				warning(nas, " (", round(100*nas/lpv,2), "%) of the presence points have NA predictor values")
			}
		} 
		
		if (! is.null(a) ) {
			a <- .getMatrix(a)
			av <- extract(x, a)
			if (is.null(dim(av))) {
				av <- matrix(av, ncol=1)
				colnames(av) <- names(x)
			}
			av <- data.frame(av)
			
			avr <- nrow(av)
			av <- stats::na.omit(av)
			nas <- length(as.vector(attr(av, "na.action")))
			if (nas > 0) {
				if (nas >= 0.5 * avr) {
					stop("more than half of the absence points have NA predictor values")
				} else {
					warning(nas, " (", round(100*nas/avr, 2), "%) of the presence points have NA predictor values")
				}
			}
		} else { 
		# random absence
			if (is.null(nbg)) {
				nbg <- 10000 
			} else {
				if (nbg < 100) {
					stop("number of background points is very low")
				} else if (nbg < 1000) {
					warning("number of background points is very low")
				}
			}

			if (nlyr(x) > 1) {
				xy <- background( rast(x,1), nbg, p, warn=0 )
			} else {
				xy <- background(x, nbg, p, warn=0 )			
			}
			av <- extract(x, xy)
			if (is.null(dim(av))) {
				av <- matrix(av, ncol=1)
				colnames(av) <- names(x)
			}
			av <- data.frame(av)
			
			
			av <- stats::na.omit(av)
			if (nrow(av) == 0) {
				stop("could not get valid background point values; is there a layer with only NA values?")
			}
			if (nrow(av) < 100) {
				stop("only got:", nrow(av), "random background point values; is there a layer with many NA values?")
			}
			if (nrow(av) < 1000) {
				warning("only got:", nrow(av), "random background point values; Small exent? Or is there a layer with many NA values?")
			}
		}
		
		# Signature = data.frame, missing

		x <- rbind(pv, av)
		
		if (!is.null(factors)) {
			for (f in factors) {
				x[,f] <- factor(x[,f])
			}
		}
		
		p <- c(rep(1, nrow(pv)), rep(0, nrow(av)))
		maxentropy(x, p, ...)	
	}
)


.getreps <- function(args) {
	if (is.null(args)) { return(1) } 
	args <- trim(args)
	i <- which(substr(args,1,10) == "replicates")
	if (! isTRUE(i > 0)) {
		return(1)
	} else {
		i <- args[i]
		i <- strsplit(i, "=")[[1]][[2]]
		return(as.integer(i))
	}
}



setMethod("maxentropy", signature(x="data.frame", p="vector"), 
	function(x, p, args=NULL, path, silent=FALSE, ...) {
	
		stopifnot(maxentropy())

		x <- cbind(p, x)
		x <- stats::na.omit(x)
		x[is.na(x)] <- -9999  # maxent flag for NA, unless changed with args(nodata= ), so we should check for that rather than use this fixed value.

		p <- x[,1]
		x <- x[, -1 ,drop=FALSE]

		factors <- NULL
		for (i in 1:ncol(x)) {
			if (class(x[,i]) == "factor") {
				factors <- c(factors, colnames(x)[i])
			}
		}
		
		if (!missing(path)) {
			path <- trim(path)
			dir.create(path, recursive=TRUE, showWarnings=FALSE)
			if (!file.exists(path)) {
				stop("cannot create output directory: ", path)
			}
			dirout <- path			
		} else {
			dirout <- .meTmpDir()
			f <- paste(round(stats::runif(10)*10), collapse="")
			dirout <- paste(dirout, "/", f, sep="")
			dir.create(dirout, recursive=TRUE, showWarnings=FALSE)
			if (! file.exists(dirout)) {
				stop("cannot create output directory: ", f)
			}
		}
		
		pv <- x[p==1, ,drop=FALSE]
		av <- x[p==0, ,drop=FALSE]
		me <- new("maxent_model")
		me@presence <- pv
		me@absence <- av
		me@hasabsence <- TRUE
		me@path <- dirout

		pv <- cbind(data.frame(species="species"), x=1:nrow(pv), y=1:nrow(pv), pv)
		av <- cbind(data.frame(species="background"), x=1:nrow(av), y=1:nrow(av), av)
		
		pfn <- paste(dirout, "/presence", sep="")
		afn <- paste(dirout, "/absence", sep="")
		utils::write.table(pv, file=pfn, sep=",", row.names=FALSE)
		utils::write.table(av, file=afn, sep=",", row.names=FALSE)

		mxe <- rJava::.jnew("mebridge")
		
		names(args) = NULL
		replicates <- .getreps(args) 
		args <- c("-z", args)

		if (is.null(factors)) {
			str <- rJava::.jcall(mxe, "S", "fit", c("autorun", "-e", afn, "-o", dirout, "-s", pfn, args)) 
		} else {
			str <- rJava::.jcall(mxe, "S", "fit", c("autorun", "-e", afn, "-o", dirout, "-s", pfn, args), rJava::.jarray(factors))
		}
		if (!is.null(str)) {
			stop("args not understood:\n", str)
		}

	
		if (replicates > 1) {
		
			mer <- new("MaxEntReplicates")
			d <- t(utils::read.csv(paste(dirout, "/maxentResults.csv", sep="") ))
			d1 <- d[1,]
			d <- d[-1, ,drop=FALSE]
			dd <- matrix(as.numeric(d), ncol=ncol(d))
			rownames(dd) <- rownames(d)
			colnames(dd) <- d1
			mer@results <- dd
			f <- paste(dirout, "/species.html", sep="")
			html <- readLines(f)
			html[1] <- "<title>Maxent model</title>"
			html[2] <- "<CENTER><H1>Maxent model</H1></CENTER>"
			html[3] <- sub("model for species", "model result", html[3])
			newtext <- paste("using 'predicts' version ", utils::packageDescription("predicts")$Version, "& Maxent version")
			html[3] <- sub("using Maxent version", newtext, html[3])
			f <- paste(dirout, "/maxent.html", sep="")
			writeLines(html, f)	
			mer@html <- f
			
			for (i in 0:(replicates-1)) {	
				mex <- me
				mex@lambdas <- unlist( readLines( paste(dirout, "/species_", i, ".lambdas", sep="") ) )
					
				f <- paste(mex@path, "/species_", i, ".html", sep="")
				html <- readLines(f)
				html[1] <- "<title>Maxent model</title>"
				html[2] <- "<CENTER><H1>Maxent model</H1></CENTER>"
				html[3] <- sub("model for species", "model result", html[3])
				newtext <- paste("using 'predicts' version ", utils::packageDescription("predicts")$Version, "& Maxent version")
				html[3] <- sub("using Maxent version", newtext, html[3])
				f <- paste(mex@path, "/maxent_", i, ".html", sep="")
				writeLines(html, f)
				mex@html <- f
				mer@models[[i+1]] <- mex
				mer@models[[i+1]]@results <- dd[, 1+1, drop=FALSE]				
			}
			
			return(mer)
			
		} else {
			
			me@lambdas <- unlist( readLines( paste(dirout, "/species.lambdas", sep="") ) )
			d <- t(utils::read.csv(paste(dirout, "/maxentResults.csv", sep="") ))
			d <- d[-1, ,drop=FALSE]
			dd <- matrix(as.numeric(d))
			rownames(dd) <- rownames(d)
			me@results <- dd
			
			f <- paste(me@path, "/species.html", sep="")
			html <- readLines(f)
			html[1] <- "<title>Maxent model</title>"
			html[2] <- "<CENTER><H1>Maxent model</H1></CENTER>"
			html[3] <- sub("model for species", "model result", html[3])
			newtext <- paste("using 'predicts' version ", utils::packageDescription("predicts")$Version, "& Maxent version")
			html[3] <- sub("using Maxent version", newtext, html[3])
			f <- paste(me@path, "/maxent.html", sep="")
			writeLines(html, f)	
			me@html <- f
		}
		
		me
	}
)


.meTmpDir <- function() {
	return( file.path(tempdir(), "maxent") )
}


.maxentRemoveTmpFiles <- function() {
	d <- .meTmpDir()
	if (file.exists(d)) {
		unlink(paste(d, "/*", sep=""), recursive = TRUE)
	}
}

setMethod("plot", signature(x="maxent_model", y="missing"), 
	function(x, sort=TRUE, main="Variable contribution", xlab="Percentage", ...) {
		r <- x@results
		rnames <- rownames(r)
		i <- grep(".contribution", rnames)
		r <- r[i, ]
		names(r) <- gsub(".contribution", "", names(r))
		if (sort) {
			r <- sort(r)
		}
		graphics::dotchart(r, main=main, xlab=xlab, ...)
		invisible(r)
	}
)



setMethod("predict", signature(object="maxent_model_replicates"), 
	function(object, x, ext=NULL, filename="", args="", ...) {
		stopifnot(maxentropy())

		n <- length(object@models)
		if (filename != "") {
			filename <- trim(filename)
			fxt <- tools::file_ext(filename)
			filename <- tools::file_path_sans_ext(filename)
			fname <- paste(filename, "_", 1:n, fxt, sep="")
		} else {
			fname <- rep("", n)
		}
		lst <- list()
		for (i in 1:n) {
			lst[[i]] <- predict(object@models[[i]], x, ext=ext, filename=fname[i], args=args, ...)
		}
		return(rast(lst))
	}
)




.maxent_predict <- function(object, mxe, args, x) {
	lambdas <- paste(object@lambdas, collapse='\n')
	variables <- colnames(object@presence)
	x <- x[,variables,drop=FALSE]
	if (inherits(x, "data.frame")) {
		for (i in 1:ncol(x)) {
			if (class(x[,i]) == "factor") {
				x[,i] <- as.numeric(as.character(x[,i]))
			} else if (class(x[,i]) == "character") {
				x[,i] <- as.numeric(x[,i])
			}
		}
	} else {
		x[] <- as.numeric(x)
	}
	
	out <- rep(NA, times=nrow(x))
	ok <- rowSums(is.na(x)) == 0
	if (sum(ok) > 0) {
		x <- as.matrix(x[ok, ,drop=FALSE])
		p <- rJava::.jcall(mxe, "[D", "predict", lambdas, rJava::.jarray(colnames(x)), rJava::.jarray(x, dispatch=TRUE), args) 
		p[p == -9999] <- NA
		out[ok] <- p
	}
	out
}



setMethod("predict", signature(object="maxent_model"), 
	function(object, x, ext=NULL, args="", filename="", ...) {

		stopifnot(maxentropy())

		args <- c(args, "")
		lambdas <- paste(object@lambdas, collapse="\n")
		variables <- colnames(object@presence)
			
		mxe <- rJava::.jnew("mebridge") 		
		args <- c("-z", args)
		tst <- rJava::.jcall(mxe, "S", "testPredictArgs", lambdas, args) 
		if (!is.null(tst)) {
			stop("args not understood:\n", tst)
		}

		if (!inherits(x, "SpatRaster")) {
			if (! all(variables %in% colnames(x))) {
				stop("missing predictor variables in x")
			}
			me <- .maxent_predict(object, mxe, args, x)
			return(me)
		}

		filename <- trimws(filename)
			
		if (! all(variables  %in%  names(x) )) {
			stop("missing layers (or wrong names)")
		}
		if (nlyr(x) > length(variables)) {
			x <- x[[variables]]
		}
		if (!is.null(ext)) {
			x <- terra::crop(x, ext)
		}
		out <- terra::rast(x, nlyr=1)
		names(out)  <- "maxent"
		ncols <- terra::ncol(out)
		if (!terra::readStart(x)) { stop(x@ptr$messages$getError()) }
		on.exit(terra::readStop(x))
		b <- terra::writeStart(out, filename, ...)
		for (i in 1:b$n) {
			rowvals <- terra::readValues(x, b$row[i], b$nrows[i], 1, ncol(x), TRUE, FALSE)
			p <- .maxent_predict(object, mxe, args, rowvals)
			terra::writeValues(out, res, b$row[i], b$nrows[i])
		}
		terra::writeStop(out)
		return(out)
	}
)




