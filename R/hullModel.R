
setClass("hull_model",
	contains = "SDM",
	representation (
		type = "character",
		polygons="SpatVector"
	),	
	prototype (	
		type="",
		polygons = vect()
	),
	validity = function(object)	{
		return(TRUE)
	}
)


setMethod("geometry", "hull_model",
	function(x) {
		x@polygons
	}
)

setMethod("plot", signature(x="hull_model", y="missing"), 
	function(x, ...) {
		terra::plot(x@polygons, ...)
	}
)

setMethod("lines", signature(x="hull_model"), 
	function(x, ...) {
		terra::lines(x@polygons, ...)
	}
)


.k_convexHulls <- function(xy, k) {

	if (k > (nrow(xy) / 2)) {
		stop('too many clusters (there should be at least twice as many points)')
	}

	cl <- stats::kmeans(xy, k, 100)$cluster
	clusters <- unique(cl)
	if (any(table(cl) < 3)) {
		stop("too many clusters (cluster with less than three points found)")
	}
	h <- vector(mode="list", length=k)
	xy <- data.frame(xy)
	for (i in clusters) {
		pts <- vect(xy[cl==i, ], geom=colnames(xy))
		h[[i]] <- convHull(pts)
	}
	vect(h)
}

.makeHullModel <- function(xy, crs="", type="circle", n=1) {
	
	type <- tolower(type)
	type <- match.arg(type, c("circle", "rectangle", "convex"))
	n <- round(n)
	
	m <- new("hull_model")
	
	stopifnot(nrow(xy) > 1)
	if (!inherits(xy, "SpatVector")) {
		m@presence <- data.frame(xy)
		lonlat <- is.lonlat(as.character(crs))
		xy <- vect(m@presence, geom=colnames(xy), crs=crs)
	} else {
		if (geomtype(p) != "points") {
			stop("SpatVector v must have points geometry")
		}
		m@presence <- data.frame(crds(xy))
		lonlat <- is.lonlat(xy)
	}
	xy <- na.omit(xy)

	if (is.na(lonlat)) {
		if (is.lonlat(xy, perhaps=TRUE)) {
			lonlat <- TRUE
			crs(xy) <- "+proj=longlat"
		} else {
			lonlat <- FALSE
			crs(xy) <- "planar"
		}
	} 

	if (type == "convex") {
		if (n > 1) {
			xy <- .k_convexHulls(xy, n)
		} else {
			xy <- convHull(xy)
		}
		if (lonlat) {
			xy <- densify(xy, 10000)
		}
		m@polygons <- xy
	} else if (type == "rectangle") {
		if (n > 1) {
			xy <- .k_convexHulls(xy, n)
			h <- vector(mode="list", length=n)
			for (i in 1:n) {
				h[[i]] <- minRect(xy[i,])
			}
			xy <- vect(h)
		} else {
			xy <- minRect(xy)
		}
		if (lonlat) {
			xy <- densify(xy, 10000)
		}
		m@polygons <- xy	
	} else if (type == "circle") {
		if (n > 1) {
			kxy <- .k_convexHulls(xy, n)
			h <- vector(mode="list", length=n)
			for (i in 1:n) {
				xy <- crds(kxy[i,])
				f <- function(p) { max(distance(rbind(p), xy, lonlat=lonlat)) }
				p <- stats::optim(colMeans(xy), f)
				v <- terra::vect(rbind(p$par), crs=crs)
				b <- buffer(v, width=p$value, quadsegs=45)
				values(b) <- data.frame(x=p$par[1], y=p$par[2], r=p$value)
				h[[i]] <- b
			}
			b <- vect(h)
		} else {
			xy <- convHull(xy)
			xy <- crds(xy)
			f <- function(p) { max(distance(rbind(p), xy, lonlat=lonlat)) }
			p <- stats::optim(colMeans(xy), f)
			v <- terra::vect(rbind(p$par), crs=crs)
			b <- buffer(v, width=p$value, quadsegs=45)
			values(b) <- data.frame(x=p$par[1], y=p$par[2], r=p$value)
		}
		m@polygons <- b
	}
	m
}



setMethod("hullModel", signature(p="data.frame"), 
	function(p, type="convex", crs=NA, n=1) {
		.makeHullModel(p, type=type, crs=crs, n=n)
	}
)

setMethod("hullModel", signature(p="matrix"), 
	function(p, type="convex", crs=NA, n=1) {
		.makeHullModel(data.frame(p), type=type, crs=crs, n=n)
	}
)

setMethod("hullModel", signature(p="SpatVector"), 
	function(p, type="convex", n=1) {
		.makeHullModel(p, type=type, crs=crs, n=n)
	}
)



setMethod("predict", signature(object="hull_model"), 
	function(object, x, ext=NULL, mask=FALSE, filename="",  ...) {
	
		nc <- nrow(object@polygons)
		if ( inherits(x, "SpatRaster"))  {
			if (! mask) {
				x <- rast(x)
			}
			if (! is.null(ext)) { 
				x <- crop(x, ext) 
			}
			
			xx <- rasterize(object@polygons, x, field=1, fun="sum")
			if (mask) {
				xx <- mask(xx, x, ...)
			}
			#if (nc > 1) {
			#	fun <- function(x){x / nc }
			#	xx <- app(xx, fun=fun, filename=filename, ...)
			#} else 
			if (filename != "") {
				xx <- writeRaster(xx, filename, ...)
			}
			return(xx)
		} else {
		
			if (!inherits(x, "SpatVector") )  {
				x <- data.frame(x)
				x <- vect(x, geom=colnames(x)[1:2])
			}
			as.integer(is.related(x, object@polygons, "coveredby"))
			# / nc
		}
	}
)

