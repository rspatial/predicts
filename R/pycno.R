
# adapted from from R package pycno
# Chris Brunsdon <christopher.brunsdon at nuim.ie>
# GPL (≥ 2)


# x = SpatRaster
# v = SpatVector 
# pop = field in v, or vector
  
pycnophy <- function(x, v, pop, r=0.2, converge=3, verbose=TRUE) {


	s1d <- function(s) unclass(filter(s,c(0.5,0,0.5)))

	pycno <- function(x, r, zones, uzones, pops) {
		mval <- mean(x)
		pad <- rbind(mval, cbind(mval, x, mval), mval)
		pad <- (t(apply(pad, 1, s1d)) + apply(pad, 2, s1d))/2
		pad <- pad[2:(nrow(x)+1), 2:(ncol(x)+1)]
		x <- x * r + (1-r) * pad
		for (i in uzones) {
			z <- (zones == i)
			x[z] <- pmax(0, x[z] * pops[i]/sum(x[z]))
			x[z] <- x[z] + (pops[i] - sum(x[z]))/sum(z)
		}
		x
	}


# prep data 
	if (is.character(pop)) {
		stopifnot(pop %in% names(v))
		pop <- v[[pop, drop=T]]
	} else {
		stopifnot(length(pop) == nrow(v))
	}
# add zero for NAs
	pops <- c(pop, 0)

# zones
	x <- rasterize(v, x, 1:nrow(v))
	zones <- t(as.matrix(x, wide=TRUE))
	zones[is.na(zones)] <- max(zones, na.rm=TRUE)+1
	uzs <- sort(unique(as.vector(zones)))
   
# pop per cell
	xx <- zones*0
	for (i in uzs) {
		zone.set <- (zones == i)
		xx[zone.set] <- pops[i]/sum(zone.set)
	}

	stopper <- max(xx) * 10^(-converge)
	repeat {
		old.x <- xx
		xx <- pycno(xx, r, zones, uzs, pops)
		if (verbose) {
			cat(sprintf("Maximum Change: %12.5f - will stop at %12.5f\n", max(abs(old.x - xx)), stopper))
			flush.console()
		}
		if (max(abs(old.x - xx)) < stopper) break 
	}

	y <- setValues(x, xx)
	mask(y, x)
} 


