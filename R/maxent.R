# Author: Robert J. Hijmans, r.hijmans@gmail.com
# Date: 2009-2021
# Version 0.1
# Licence GPL v3

setClass("MaxEnt_model",
	contains = "SDM",
	representation (
		lambdas  = "vector",
		results = "matrix",
		path = "character",
		html = "character",
		levels = "list"		
	),	
	prototype (	
		lambdas = as.vector(NA),
		results = as.matrix(NA),
		path = "",
		html = "",
		levels = list()
	),
)


setClass("MaxEnt_model_replicates",
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


.getMeVersion <- function(...) {}

setMethod("MaxEnt", signature(x="missing", p="missing"), 
	function(x, p, silent=FALSE, ...) {

		if (is.null(getOption("predicts_rJavaLoaded"))) {
			# to avoid trouble on macs
			Sys.setenv(NOAWT=TRUE)
			if ( requireNamespace("rJava") ) {
				rJava::.jpackage("predicts")
				options(predicts_rJavaLoaded=TRUE)
			} else {
				if (!silent) {
					message("Cannot load rJava")			
				}
				return(FALSE)
			}
		}

		if (is.null(getOption("predicts_maxent"))) {
			mxe <- rJava::.jnew("meversion") 
			v <- try(rJava::.jcall(mxe, "S", "meversion"), silent=TRUE)
			if (inherits(v, "try-error")) {
				if (!silent) {
					message("MaxEnt_model is missing or incompatible with your version of Java")
				}
				return(FALSE)
			} else if (v == "3.3.3a") {
				if (!silent) {
					message("This is not a compatible version of Maxent")
				}
				return(FALSE)
			}
			options(predicts_maxent=v)
		} else {
			v = getOption("dismo_maxent")
		}
		if (!silent) {
			message(paste("This is MaxEnt version", v))
		}
		invisible(TRUE)
	}
)



.getMatrix <- function(x) {
	if (inherits(x, "SpatVector")) {
		x <- geom(x)[,c("x", "y"), drop=FALSE]
	} 
	if (inherits(x, "data.frame")) {
		x <- as.matrix(x)
	}
	if (!inherits(x, "matrix" )) {
		stop("data should be  a matrix, data.frame, or SpatVector")
	}
	if (dim(x)[2] != 2) {
		stop("presence or background coordinates data should be a matrix or data.frame with 2 columns" ) 	
	}
	colnames(x) <- c("x", "y")
	return(x)
} 


.maxentFailedMsg <- function(dirout) {
	logf <- file.path(dirout, "maxent.log")
	if (!file.exists(logf)) {
		return(paste0("MaxEnt failed (no log file at '", logf, "')"))
	}
	lines <- tryCatch(readLines(logf, warn=FALSE), error=function(e) character())
	errlines <- grep("^(Error|Exception)", lines, value=TRUE)
	if (length(errlines) > 0) {
		return(paste0("MaxEnt failed:\n  ", paste(errlines, collapse="\n  "),
			"\n(full log: '", logf, "')"))
	}
	if (length(lines) == 0) {
		return(paste0("MaxEnt failed (empty log at '", logf, "')"))
	}
	n <- min(20, length(lines))
	paste0("MaxEnt failed. Tail of '", logf, "':\n  ",
		paste(utils::tail(lines, n), collapse="\n  "))
}


.biasBackground <- function(x, nbg, biasfile) {
	if (!inherits(biasfile, "SpatRaster")) {
		stop("'biasfile' must be a SpatRaster")
	}
	biasfile <- biasfile[[1]]
	compareGeom(biasfile, x[[1]], stopOnError=TRUE)
	bpts <- spatSample(biasfile, nbg, method="weights", na.rm=TRUE, as.points=TRUE, values=FALSE, warn=FALSE)
	av <- extract(x, bpts, ID=FALSE)
	stats::na.omit(av)
}


#factors=NULL, 
setMethod("MaxEnt", signature(x="SpatRaster", p="ANY"), 
	function(x, p, a=NULL, removeDuplicates=TRUE, nbg=10000, biasfile=NULL, ...) {

		dots <- list(...)
		if (!is.null(dots$args)) { #redundant with <data.frame> method check but nice to catch early
			bi <- grep("^biasfile=", trimws(dots$args), ignore.case=TRUE)
			if (length(bi) > 0) {
				stop("args='biasfile=..' cannot be used here. Either use the 'biasfile' argument directly\nor supply a bias-weighted background via 'a'.", call.=FALSE)
			}
		}
		p <- predicts:::.getMatrix(p)
		if (removeDuplicates) {
			cells <- unique(cellFromXY(x, p))
			pv <- extract(x, cells)
		} else {
			pv <- extract(x, p)
		}

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
		
		if (!is.null(a) ) {
			if (!is.null(biasfile)) {
				stop("'biasfile' cannot be used as background points 'a' were supplied (draw 'a' with the bias grid as sampling weights?)") 
			}
			a <- .getMatrix(a)
			av <- extract(x, a)
			avr <- nrow(av)
			av <- stats::na.omit(av)
			nas <- length(as.vector(attr(av, "na.action")))
			if (nas > 0) {
				if (nas >= 0.5 * avr) {
					stop("more than half of the background points have NA predictor values")
				} else {
					warning(nas, " (", round(100*nas/avr, 2), "%) of the background points have NA predictor values")
				}
			}
		} else { 
		# (bias-weighted) random background
			if (is.null(nbg)) {
				nbg <- 10000 
			} else {
				if (nbg < 100) {
					stop("number of background points too low")
				} else if (nbg < 1000) {
					warning("number of background points is very low")
				}
			}
			if (is.null(biasfile)) {
				av <- spatSample(x, nbg, "random", na.rm=TRUE, warn=FALSE)
			} else {
				av <- .biasBackground(x, nbg, biasfile)
			}
			if (nrow(av) < 100) {
				stop("only got: ", nrow(av), " random background point values; is there a layer with many NA values?")
			}
			if (nrow(av) < 1000) {
				warning("only got: ", nrow(av), " random background point values; Small exent? Or is there a layer with many NA values?")
			}
		}
		
		# Signature = data.frame, missing

		v <- rbind(pv, av)
		
#		if (!is.null(factors)) {
#		  for (f in factors) {
#		    x[,f] <- factor(x[,f])
#		  }
#		}
		
		p <- c(rep(1, nrow(pv)), rep(0, nrow(av)))
		MaxEnt(v, p, ...)	
	}
)


.getreps <- function(args) {
	if (is.null(args)) { return(1) } 
	args <- trimws(args)
	i <- which(substr(args,1,10) == "replicates")
	if (! isTRUE(i > 0)) {
		return(1)
	} else {
		i <- args[i]
		i <- strsplit(i, "=")[[1]][[2]]
		return(as.integer(i))
	}
}



setMethod("MaxEnt", signature(x="data.frame", p="numeric"), 
	function(x, p, args=NULL, path, silent=FALSE, ...) {
	
		stopifnot(MaxEnt(silent=TRUE))

		x <- cbind(p, x)
		x <- stats::na.omit(x)

		bi <- grep("^biasfile=", trimws(args), ignore.case=TRUE)
		if (length(bi) > 0) {
			stop("'biasfile=' cannot be used here. Either use the 'biasfile' argument\nin the MaxEnt<SpatRaster> method or supply a bias-weighted background via 'a'.", call.=FALSE)
		}

		nodata <- grep("nodata=", args, value=TRUE)
		if (length(nodata) == 1) {
			nodata = as.numeric(strsplit(nodata, "=")[[1]][2])
			x[is.na(x)] <- nodata			
		} else {
			x[is.na(x)] <- -9999 
		}

		p <- x[,1]
		x <- x[, -1 ,drop=FALSE]

		me <- new("MaxEnt_model")

		me@presence <- x[p==1, ,drop=FALSE]
		me@absence <- x[p==0, ,drop=FALSE]

		factors <- names(x)[sapply(x, is.factor)]
		for (f in factors) {
			me@levels[[f]] <- levels(x[[f]])
			x[f] <- as.integer(x[[f]])
		}
		
		if (!missing(path)) {
			path <- trimws(path)
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
		#file.remove(file.path(dirout, "maxent.html"))

		
		me@path <- dirout
		
		pv <- x[p==1, ,drop=FALSE]
		av <- x[p==0, ,drop=FALSE]
		#me@presence <- pv
		#me@absence <- av
		
		pv <- cbind(data.frame(species="species"), x=1:nrow(pv), y=1:nrow(pv), pv)
		av <- cbind(data.frame(species="background"), x=1:nrow(av), y=1:nrow(av), av)
		
		pfn <- paste(dirout, "/presence", sep="")
		afn <- paste(dirout, "/background", sep="")
		utils::write.table(pv, file=pfn, sep=",", row.names=FALSE)
		utils::write.table(av, file=afn, sep=",", row.names=FALSE)

		mxe <- rJava::.jnew("mebridge")
		
		names(args) = NULL
		replicates <- predicts:::.getreps(args) 
		args <- c("-z", args)

		# factors = NULL
		if (is.null(factors)) {
			str <- rJava::.jcall(mxe, "S", "fit", c("autorun", "-e", afn, "-o", dirout, "-s", pfn, args)) 
		} else {
			str <- rJava::.jcall(mxe, "S", "fit", c("autorun", "-e", afn, "-o", dirout, "-s", pfn, args), rJava::.jarray(factors))
		}
		if (!is.null(str)) {
			stop(paste0("args not understood:\n", str))
		}

	
		if (replicates > 1) {
		
			mer <- new("MaxEnt_model_replicates")
			d <- t(utils::read.csv(file.path(dirout, "maxentResults.csv")))
			d1 <- d[1,]
			d <- d[-1, ,drop=FALSE]
			dd <- matrix(as.numeric(d), ncol=ncol(d))
			rownames(dd) <- rownames(d)
			colnames(dd) <- d1
			mer@results <- dd
			f <- file.path(dirout, "species.html")
			html <- readLines(f)
			html[1] <- "<title>Maxent model</title>"
			html[2] <- "<CENTER><H1>Maxent model</H1></CENTER>"
			html[3] <- sub("model for species", "model result", html[3])
			newtext <- paste("using 'predicts' version ", utils::packageDescription("predicts")$Version, "& Maxent version")
			html[3] <- sub("using Maxent version", newtext, html[3])
			f <- file.path(dirout, "maxent.html")
			writeLines(html, f)	
			mer@html <- normalizePath(f)
						
			for (i in 0:(replicates-1)) {	
				mex <- me
				mex@lambdas <- unlist( readLines( file.path(dirout, paste0("species_", i, ".lambdas") ) ))
					
				f <- file.path(mex@path, paste0("species_", i, ".html"))
				html <- readLines(f)
				html[1] <- "<title>Maxent model</title>"
				html[2] <- "<CENTER><H1>Maxent model</H1></CENTER>"
				html[3] <- sub("model for species", "model result", html[3])
				newtext <- paste("using 'predicts' version ", utils::packageDescription("predicts")$Version, "& Maxent version")
				html[3] <- sub("using Maxent version", newtext, html[3])
				f <- paste(mex@path, "/maxent_", i, ".html", sep="")
				writeLines(html, f)
				mex@html <- normalizePath(f)
				mer@models[[i+1]] <- mex
				mer@models[[i+1]]@results <- dd[, 1+1, drop=FALSE]				
			}
			
			return(mer)
			
		} else {
			
			flambdas <- file.path(dirout, "species.lambdas")
			if (!file.exists(flambdas)) {
				stop(.maxentFailedMsg(dirout), call.=FALSE)
			}
			me@lambdas <- unlist( readLines( flambdas ) )
			d <- t(utils::read.csv(file.path(dirout, "maxentResults.csv") ))
			d <- d[-1, ,drop=FALSE]
			d[d=="na"] <- NA
			dd <- matrix(as.numeric(d))
			rownames(dd) <- rownames(d)
			me@results <- dd
			
			f <- file.path(me@path, "species.html")
			html <- readLines(f)
			html[1] <- "<title>Maxent model</title>"
			html[2] <- "<CENTER><H1>Maxent model</H1></CENTER>"
			html[3] <- sub("model for species", "model result", html[3])
			newtext <- paste("using 'predicts' version ", utils::packageDescription("predicts")$Version, "& Maxent version")
			html[3] <- sub("using Maxent version", newtext, html[3])
			f <- file.path(me@path, "maxent.html")
			writeLines(html, f)	
			me@html <- normalizePath(f)
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

setMethod("plot", signature(x="MaxEnt_model"), 
	function(x, ...) {
		r <- x@results
		rnames <- rownames(r)
		i <- grep(".contribution", rnames)
		r <- r[i, ]
		names(r) <- gsub(".contribution", "", names(r))
		r <- sort(r)
		dots <- list(...)
		if (is.null(dots$main)) dots$main="Variable contribution"
		if (is.null(dots$xlab)) dots$xlab="Percentage"
		dots$x <- r
		do.call(graphics::dotchart, dots)
		invisible(r)
	}
)



setMethod("predict", signature(object="MaxEnt_model_replicates"), 
	function(object, x, ext=NULL, filename="", args="", ...) {
		stopifnot(MaxEnt(silent=TRUE))

		n <- length(object@models)
		if (filename != "") {
			filename <- trimws(filename)
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
	x <- x[, variables, drop=FALSE]
	if (inherits(x, "data.frame")) {
		for (nm in names(object@levels)) {
			if (inherits(x[[nm]], "factor")) {
				if ((length(levels(x[[nm]])) == length(object@levels[[nm]])) && (all(levels(x[[nm]]) == object@levels[[nm]]))) {
					x[nm] <- as.integer(x[[nm]])
				} else {
					stop(paste0("levels of ", nm, " do not match the levels of the data used to fit model"))
				}
			} else {
				stop(paste("expecting a factor or character value for", nm))
			}	
		}
	}
	#else {
		#x[] <- as.numeric(x)
	#}
	
	out <- rep(NA, times=nrow(x))
	ok <- rowSums(is.na(x)) == 0
	if (sum(ok) > 0) {
		x <- as.matrix(x[ok, ,drop=FALSE])
		p <- rJava::.jcall(mxe, "[D", "predict", lambdas, rJava::.jarray(colnames(x)), rJava::.jarray(x, dispatch=TRUE), args)
		p[p < -9999] <- NA
		out[ok] <- p
	}
	out
}



setMethod("predict", signature(object="MaxEnt_model"), 
	function(object, x, ext=NULL, args="", filename="", ...) {

		stopifnot(MaxEnt(silent=TRUE))

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
		terra::readStart(x)
		on.exit(terra::readStop(x))
		b <- terra::writeStart(out, filename, ...)
		for (i in 1:b$n) {
			rowvals <- terra::readValues(x, b$row[i], b$nrows[i], 1, ncol(x), FALSE, TRUE)
			p <- .maxent_predict(object, mxe, args, rowvals)
			terra::writeValues(out, p, b$row[i], b$nrows[i])
		}
		terra::writeStop(out)
		return(out)
	}
)
		
		

setMethod ("show" , "MaxEnt_model", 
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



setMethod ("show" , "MaxEnt_model_replicates", 
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

		
