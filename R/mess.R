# author: Jean-Pierre Rossi <jean-pierre.rossi@supagro.inra.fr>
# modifications by Robert Hijmans and Paulo van Breugel


.messi3 <- function(p,v) {
	v <- stats::na.omit(v)
	f <- 100*findInterval(p, sort(v)) / length(v)
	minv <- min(v)
	maxv <- max(v)
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
	a <- ecdf(v)(p)
	a[a>0.5] <- 1-a[a>0.5]
	200 * a
}


if (!isGeneric("mess")) { setGeneric("mess", function(x, ...) standardGeneric("mess")) }	

setMethod("mess", signature(x="SpatRaster"), 
	function(x, v, full=FALSE, filename="", ...) {

		stopifnot(NCOL(v) == nlyr(x))
		out <- rast(x)
		nl <- nlyr(x)
		filename <- trim(filename)
		nms <- paste0(names(x), "_mess")
		readStart(x)
		on.exit(readStop(x))
		if (nl == 1) {
			names(out) <- "mess"
			b <- writeStart(out, filename, ...)
			for (i in 1:b$n) {		
				vv <- terra::readValues(x, b$row[i], b$nrows[i])
				p <- .messi3(vv, v)
				terra::writeValues(out, p, b$row[i], b$nrows[i])
			}
		} else {
			if (full) {
				nlyr(out) <- nl+1
				names(out) <- c(nms, "mess")
				b <- writeStart(out, filename, ...)
				for (i in 1:b$n) {
					vv <- terra::readValues(x, b$row[i], b$nrows[i])
					vv <- sapply(1:ncol(v), function(i) .messi3(vv[,i], v[,i]))
					m <- apply(vv, 1, min, na.rm=TRUE)
					terra::writeValues(out, cbind(vv, m), b$row[i], b$nrows[i])
				}
			} else {			
				names(out) <- "mess"
				b <- writeStart(out, filename, ...)
				for (i in 1:b$n) {
					vv <- terra::readValues(x, b$row[i], b$nrows[i])
					vv <- sapply(1:ncol(v), function(i) .messi3(vv[,i], v[,i]))
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
			data.frame(mess=.messi3(x, v))
		} else {
			x <- sapply(1:ncol(x), function(i) .messi3(x[,i], v[,i]))
			rmess <- apply(x, 1, min, na.rm=TRUE)
			if (full) {
				out <- data.frame(x, rmess)
				names(out) <- c(nms, "mess")
				out
			} else {
				data.frame(mess=rmess)
			}
		}	
	}
)


