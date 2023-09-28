# author: Jean-Pierre Rossi <jean-pierre.rossi@supagro.inra.fr>
# modifications by Robert Hijmans and Paulo van Breugel
# rewritten for predicts by RH

.messi <- function(p, v) {
	
	v <- sort(v)
	f <- 100 * findInterval(p, v) / length(v)
	minv <- v[1]
	maxv <- v[length(v)]
	res <- 2*f 
	f[is.na(f)] <- -99
	i <- f>50 & f<100
	res[i] <- 200-res[i]

	i <- f==0 
	res[i] <- 100*(p[i]-minv)/(maxv-minv)
	i <- f==100
	res[i] <- 100*(maxv-p[i])/(maxv-minv)
	res
}


.messix <- function(p,v) {
# a little bit different, no negative values.
	a <- stats::ecdf(v)(p)
	a[a>0.5] <- 1-a[a>0.5]
	200 * a
}



setMethod("mess", signature(x="SpatRaster"), 
	function(x, v, full=FALSE, filename="", ...) {

		if (inherits(v, "SpatVector")) {
			if (geomtype(p) != "points") {
				stop("SpatVector v must have points geometry")
			}
			v <- extract(v, x)
		}
		v <- stats::na.omit(v)
		if (nrow(v) < 2) {
			stop("insufficient number of reference points")
		}
		stopifnot(NCOL(v) == nlyr(x))

		out <- rast(x)
		nl <- nlyr(x)
		nms <- paste0(names(x), "_mess")
		readStart(x)
		on.exit(readStop(x))
		if (nl == 1) {
			names(out) <- "mess"
			b <- writeStart(out, filename, ...)
			for (i in 1:b$n) {		
				vv <- terra::readValues(x, b$row[i], b$nrows[i])
				p <- .messi(vv, v)
				terra::writeValues(out, p, b$row[i], b$nrows[i])
			}
		} else {
			if (full) {
				nlyr(out) <- nl+1
				names(out) <- c(nms, "mess")
				b <- writeStart(out, filename, ...)
				for (i in 1:b$n) {
					vv <- terra::readValues(x, b$row[i], b$nrows[i], mat=TRUE)
					vv <- sapply(1:ncol(v), function(i) .messi(vv[,i], v[,i]))
					m <- apply(vv, 1, min, na.rm=TRUE)
					terra::writeValues(out, cbind(vv, m), b$row[i], b$nrows[i])
				}
			} else {			
				nlyr(out) <- 1
				names(out) <- "mess"
				b <- writeStart(out, filename, ...)
				for (i in 1:b$n) {
					vv <- terra::readValues(x, b$row[i], b$nrows[i], mat=TRUE)
					vv <- sapply(1:ncol(v), function(i) .messi(vv[,i], v[,i]))
					m <- apply(vv, 1, min, na.rm=TRUE)
					terra::writeValues(out, m, b$row[i], b$nrows[i])
				}
			}
		}
		writeStop(out)
		out
	}	
)

setMethod("mess", signature(x="data.frame"), 
	function(x, v, full=FALSE) {
		if (ncol(x) == 1) {
			data.frame(mess=.messi(x, v))
		} else {
			x <- sapply(1:ncol(x), function(i) .messi(x[,i], v[,i]))
			rmess <- apply(x, 1, min, na.rm=TRUE)
			if (full) {
				out <- data.frame(x, rmess)
				nms <- paste0(names(x), "_mess")
				names(out) <- c(nms, "mess")
				out
			} else {
				data.frame(mess=rmess)
			}
		}	
	}
)


