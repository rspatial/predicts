# Author: Robert J. Hijmans
# Date : December 2009
# Version 0.1
# Licence GPL v3


backgroundSample <- function(mask, n, p, ext=NULL, extf=1.1, excludep=TRUE, cellnumbers=FALSE, tryf=3, warn=2) {
	
	mask <- rast(mask)[[1]]	
	tryf <- max(round(tryf[1]), 1)
	
	if (missing(p)) { 
		excludep <- FALSE
	} else {
		if (inherits(p, 'SpatVector')) {
			p <- geom(p)[, c("x", "y")]
		}
	}
	
	if (inherits(ext, "character")) {
		if (! (ext %in% 'points') ) { 
			stop("if ext is a character variable it should be 'points'") 
		} else if (missing(p)) { 
			warning("if p is missing, 'ext=points' is meaningless") 
			ext <- ext(mask)  
		} else {
			ext <- ext(min(p[,1]), max(p[,1]), min(p[,2]), max(p[,2]))
		}
	} 

	if (! is.null(ext)) {
		ext <- ext(ext)
		ext <- ext * extf
		ext <- intersect(ext, ext(mask))
		mask2 <- crop(rast(mask), ext)
	}  else {
		mask2 <- rast(mask)
	}
	
	if (n > ncell(mask2)) {
		n <- ncell(mask2)
		if (warn>0) { warning('changed n to ncell of the mask (extent)') }
	}
	
	nn <- n * tryf
	nn <- max(nn, 10)

	nn <- min(ncell(mask2), nn)
	cells <- as.vector(spatSample(mask2, nn, cells=TRUE))
	if (hasValues(mask)) {
		vals <- cbind(cells, extract(mask, cells))
		cells <- stats::na.omit(vals)[,1]
	}
		
	if (excludep) {	
		pcells <- cellFromXY(mask, p)
		cells <- cells[ ! (cells %in% pcells) ] 	
	}

	if (length(cells) > n) { 
		
		cells <- sample(cells, n) 
		
	} else if (length(cells) < n) {
	
		frac <- length(cells) / n
		if (frac < 0.1) {
			stop("generated random points = ", frac," times requested number; Use a higher value for tryf" )
		}
		if (frac < 0.5  & warn==1) {
			warning("generated random points = ", frac," times requested number; Use a higher value for tryf" )
		} else if (warn > 1) {
			warning("generated random points = ", frac," times requested number")
		}
	}
	
	if (cellnumbers) {
		return(cells)
	} else {
		return(xyFromCell(mask, cells))
	}
}

