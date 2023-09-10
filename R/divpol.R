
divider <- function(x, n, env=NULL, alpha=1, ...) {
	stopifnot(geomtype(x) == "polygons")
	n <- round(n)
	stopifnot(n > 0)
	if (n == 1) return(deepcopy(x))
	xcrs <- crs(x) 
	crs(x) <- "+proj=utm +zone=1"
	s <- terra::spatSample(x, max(n*4, 1000, log(n) * 100), "regular")
	xy <- terra::crds(s)
	if (!is.null(env)) {
		e <- extract(env, s, ID=FALSE)
		alpha <- rep_len(alpha, 2)
		xy[,1] <- xy[,1] * alpha[1]
		xy[,2] <- xy[,2] * alpha[2]
		xy <- na.omit(cbind(xy, e))
	} 
	k <- stats::kmeans(xy, centers = n, ...)
	ctrs <- k$centers[, 1:2]
	if (!is.null(env)) {
		ctrs[,1] <- ctrs[,1] / alpha[1]
		ctrs[,2] <- ctrs[,2] / alpha[2]
	}
	v <- terra::voronoi(vect(ctrs), bnd=x)
	v <- terra::crop(v, x)
	crs(v) <- xcrs
	v
}



stripper <- function(x, f=c(1/3, 2/3), vertical=TRUE){

	if ((any(f <= 0)) || (any(f >= 1))) {
		stop("f values must be > 0 and < 1")
	}
	u <- unique(f)
	if (length(u) < length(f)) {
		stop("f values must be unique")	
	}
	ord <- order(f)
	if (!all(ord == 1:length(ord))) {
		stop("f values must be in ascending order")		
	}
	
	totalArea <- expanse(x)
	e <- ext(x)
	ex <- data.frame(t(as.vector(e + 1)))
	e <- data.frame(t(as.vector(e)))
	
	if (vertical) {
		edges <- sapply(f, function(fraction){
			target <- totalArea * fraction
			target_function <- function(xm){
				expanse(crop(x, ext(ex$xmin, xm, ex$ymin, ex$ymax))) - target
			}
			stats::uniroot(target_function, lower=e$xmin, upper=e$xmax)$root
		})
		bnds <- matrix(c(ex$xmin, rep(edges,rep(2,length(edges))), ex$xmax), ncol=2, byrow=TRUE)
		a <- apply(bnds, 1, function(edges){ 
				r <- ext(edges[1], edges[2], ex$ymin, ex$ymax)
				crop(x,r)
			})
	} else {
		edges <- sapply(f, function(fraction){
			target <- totalArea * fraction
			target_function <- function(ym){
				expanse(crop(x, ext(ex$xmin, ex$xmax, ex$ymin, ym))) - target
			}
			stats::uniroot(target_function, lower=e$ymin, upper=e$ymax)$root
		})
		bnds <- matrix(c(ex$ymin, rep(edges,rep(2,length(edges))), ex$ymax), ncol=2, byrow=TRUE)
		a <- apply(bnds, 1, function(edges){ 
				r <- ext(ex$xmin, ex$xmax, edges[1], edges[2])
				crop(x,r)
			})
	}
				
	vect(a)
}


