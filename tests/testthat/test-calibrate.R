test_that("calibrated under the null with independent genes", {
  set.seed(42)
  e <- matrix(rnorm(3000 * 60), 3000, 60,
              dimnames = list(paste0("g", 1:3000), paste0("s", 1:60)))
  grp <- rep(c(TRUE, FALSE), each = 30)
  r <- calibrate_geneset(e, list(a = paste0("g", 1:100)), grp,
                         n_null = 100, n_split = 150)
  expect_gt(r$p_empirical, 0.05)
  expect_true(r$ratio > 0.5 && r$ratio < 2)
})

test_that("a real effect is detected", {
  set.seed(42)
  e <- matrix(rnorm(3000 * 60), 3000, 60,
              dimnames = list(paste0("g", 1:3000), paste0("s", 1:60)))
  grp <- rep(c(TRUE, FALSE), each = 30)
  e[1:100, grp] <- e[1:100, grp] + 1
  r <- calibrate_geneset(e, list(a = paste0("g", 1:100)), grp,
                         n_null = 100, n_split = 150)
  expect_lt(r$p_empirical, 0.05)
})

test_that("a coherent set inflates the variance ratio", {
  set.seed(42)
  lat <- rnorm(60)
  e <- matrix(rnorm(3000 * 60), 3000, 60,
              dimnames = list(paste0("g", 1:3000), paste0("s", 1:60)))
  e[1:100, ] <- 0.9 * matrix(rep(lat, each = 100), 100, 60) + 0.45 * e[1:100, ]
  grp <- rep(c(TRUE, FALSE), each = 30)
  r <- calibrate_geneset(e, list(a = paste0("g", 1:100)), grp,
                         n_null = 100, n_split = 150)
  expect_gt(r$ratio, 1.05)
  expect_true(r$sigma_across > r$sigma_within)
})

test_that("residual coherence ignores a pure group effect", {
  set.seed(42)
  e <- matrix(rnorm(2000 * 60), 2000, 60,
              dimnames = list(paste0("g", 1:2000), NULL))
  grp <- rep(c(TRUE, FALSE), each = 30)
  e[1:80, grp] <- e[1:80, grp] + 1
  raw <- set_coherence(e, paste0("g", 1:80))
  res <- set_coherence(e, paste0("g", 1:80), grp)
  expect_gt(raw, 0.1)
  expect_lt(abs(res), 0.05)
})

test_that("input validation works", {
  e <- matrix(rnorm(500 * 10), 500, 10, dimnames = list(paste0("g", 1:500), NULL))
  expect_error(calibrate_geneset(e, list(a = paste0("g", 1:20)), rep(TRUE, 5)))
  expect_error(calibrate_geneset(e, list(paste0("g", 1:20)), rep(c(TRUE, FALSE), 5)))
})
