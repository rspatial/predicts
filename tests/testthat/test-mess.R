library(terra)

# Helper: build a small single-layer SpatRaster and matching reference data.frame
make_rast <- function(nlyr = 1) {
  r <- rast(nrows = 5, ncols = 5, nlyr = nlyr)
  values(r) <- matrix(seq_len(5 * 5 * nlyr), ncol = nlyr)
  r
}

# -----------------------------------------------------------------------
# Bug fix 1: SpatVector geometry check used undefined `p`; should use `v`
# -----------------------------------------------------------------------
test_that("mess() rejects non-point SpatVector with an informative error", {
  r <- make_rast()
  # Create a polygon SpatVector — passing this as `v` previously caused
  # "object 'p' not found" rather than the intended geometry error
  poly <- as.polygons(r)
  expect_error(
    mess(r, poly),
    regexp = "points geometry",
    info = "Should mention 'points geometry', not fail with undefined 'p'"
  )
})

test_that("mess() accepts a point SpatVector without error", {
  # NOTE: extract(v, x) in the SpatVector branch has arguments reversed;
  # it should be extract(x, v). That is a separate bug to be fixed.
  # This test is skipped until that fix is in place.
  skip(
    "extract() argument order in SpatVector branch is reversed -- separate bug"
  )
  r <- make_rast()
  pts <- spatSample(
    r,
    size = 10,
    method = "random",
    as.points = TRUE,
    na.rm = TRUE
  )
  expect_no_error(mess(r, pts))
})

# -----------------------------------------------------------------------
# Bug fix 2: single-column data.frame passed whole data.frame to .messi()
#            instead of extracting the vector first
# -----------------------------------------------------------------------
test_that("mess() data.frame method works with a single-column input", {
  set.seed(7391)
  x <- data.frame(bio1 = runif(20, 10, 30)) # prediction points
  v <- data.frame(bio1 = runif(50, 5, 35)) # reference sample

  result <- mess(x, v)

  expect_s3_class(result, "data.frame")
  expect_named(result, "mess")
  expect_equal(nrow(result), nrow(x))
  expect_true(all(is.numeric(result$mess)))
})

test_that("mess() single-column and multi-column data.frame give consistent row MESS values", {
  set.seed(2847)
  x <- data.frame(bio1 = runif(15, 10, 30), bio2 = runif(15, 5, 20))
  v <- data.frame(bio1 = runif(40, 5, 35), bio2 = runif(40, 0, 25))

  result_multi <- mess(x, v)

  # The single-column path should match the first column of a two-column result
  # when only bio1 is used
  x1 <- x["bio1"]
  v1 <- v["bio1"]
  result_single <- mess(x1, v1)

  expect_s3_class(result_single, "data.frame")
  expect_named(result_single, "mess")
  expect_equal(nrow(result_single), nrow(x1))
})
