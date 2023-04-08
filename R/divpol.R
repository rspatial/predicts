
divide_polygons <- function(x, n, ...) {
	stopifnot(geomtype(x) == "polygons")
	n <- round(n)
	stopifnot(n > 0)
	if (n == 1) return(deepcopy(x))
	xcrs <- crs(x) 
	crs(x) <- "+proj=utm +zone=1"
	s <- terra::spatSample(x, max(n*4, 1000, log(n) * 100), "regular")
	s <- terra::crds(s)
	k <- stats::kmeans(s, centers = n, ...)
	v <- terra::voronoi(vect(k$centers), bnd=x)
	v <- terra::crop(v, x)
	crs(v) <- xcrs
	v
}


