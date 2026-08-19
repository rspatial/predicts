# author: Jean-Pierre Rossi <jean-pierre.rossi@supagro.inra.fr>
# modifications by Robert Hijmans and Paulo van Breugel
# rewritten for predicts by RH

# Internal helper: compute MESS for one layer given pre-sorted reference vector.
# Keeping sort() out of this function avoids re-sorting on every block iteration.
.messi_sorted <- function(p, sv) {
  # sv must already be sorted with no NAs
  n <- length(sv)
  minv <- sv[1]
  maxv <- sv[n]
  f <- 100 * findInterval(p, sv) / n
  res <- 2 * f
  f[is.na(f)] <- -99
  i <- f > 50 & f < 100
  res[i] <- 200 - res[i]
  i <- f == 0
  res[i] <- 100 * (p[i] - minv) / (maxv - minv)
  i <- f == 100
  res[i] <- 100 * (maxv - p[i]) / (maxv - minv)
  res
}

# Public-facing wrapper that sorts before delegating; used by the data.frame
# method and any callers that pass unsorted reference vectors directly.
.messi <- function(p, v) {
  .messi_sorted(p, sort(v))
}


.messix <- function(p, v) {
  # a little bit different, no negative values.
  a <- stats::ecdf(v)(p)
  a[a > 0.5] <- 1 - a[a > 0.5]
  200 * a
}


setMethod(
  "mess",
  signature(x = "SpatRaster"),
  function(x, v, full = FALSE, filename = "", ...) {
    if (inherits(v, "SpatVector")) {
      if (geomtype(v) != "points") {
        stop("SpatVector v must have points geometry")
      }
      v <- extract(x, v, ID = FALSE)
    }
    v <- stats::na.omit(v)
    v <- as.matrix(v)
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
      # Pre-sort the single reference column once, outside the block loop
      sv <- sort(v[, 1])
      names(out) <- "mess"
      b <- writeStart(out, filename, ...)
      for (i in 1:b$n) {
        vv <- terra::readValues(x, b$row[i], b$nrows[i])
        terra::writeValues(out, .messi_sorted(vv, sv), b$row[i], b$nrows[i])
      }
    } else {
      # Pre-sort each reference column once, outside the block loop
      sv <- lapply(1:ncol(v), function(j) sort(v[, j]))

      if (full) {
        nlyr(out) <- nl + 1
        names(out) <- c(nms, "mess")
      } else {
        nlyr(out) <- 1
        names(out) <- "mess"
      }
      b <- writeStart(out, filename, ...)
      for (i in 1:b$n) {
        vv <- terra::readValues(x, b$row[i], b$nrows[i], mat = TRUE)
        mm <- vapply(
          1:ncol(v),
          function(j) .messi_sorted(vv[, j], sv[[j]]),
          numeric(nrow(vv))
        )
        suppressWarnings(m <- apply(mm, 1, min, na.rm = TRUE))
        m[!is.finite(m)] <- NA
        terra::writeValues(
          out,
          if (full) cbind(mm, m) else m,
          b$row[i],
          b$nrows[i]
        )
      }
    }
    writeStop(out)
    out
  }
)

setMethod("mess", signature(x = "data.frame"), function(x, v, full = FALSE) {
  if (ncol(x) == 1) {
    data.frame(mess = .messi(x[, 1], v[, 1]))
  } else {
    mm <- vapply(
      1:ncol(x),
      function(i) .messi(x[, i], v[, i]),
      numeric(nrow(x))
    )
    rmess <- apply(mm, 1, min, na.rm = TRUE)
    if (full) {
      out <- data.frame(mm, rmess)
      nms <- paste0(names(x), "_mess")
      names(out) <- c(nms, "mess")
      out
    } else {
      data.frame(mess = rmess)
    }
  }
})
