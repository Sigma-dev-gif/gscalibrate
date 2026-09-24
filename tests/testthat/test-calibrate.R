test_that('independent genes give calibrated p-values under the null', {
  set.seed(42)
  e <- matrix(rnorm(3000 * 60), 3000, 60,
              dimnames = list(paste0('g', 1:3000), paste0('s', 1:60)))
  grp <- rep(c(TRUE, FALSE), each = 30)
  r <- calibrate_geneset(e, list(a = paste0('g', 1:100)), grp, n_null = 100)
  expect_gt(r$p_empirical, 0.05)
  expect_true(r$reliable)
})

test_that('a real effect is detected', {
  set.seed(42)
  e <- matrix(rnorm(3000 * 60), 3000, 60,
              dimnames = list(paste0('g', 1:3000), paste0('s', 1:60)))
  grp <- rep(c(TRUE, FALSE), each = 30)
  e[1:100, grp] <- e[1:100, grp] + 1
  r <- calibrate_geneset(e, list(a = paste0('g', 1:100)), grp, n_null = 100)
  expect_lt(r$p_empirical, 0.05)
})

test_that('a coherent set is flagged unreliable', {
  set.seed(42)
  lat <- rnorm(60)
  e <- matrix(rnorm(3000 * 60), 3000, 60,
              dimnames = list(paste0('g', 1:3000), paste0('s', 1:60)))
  e[1:100, ] <- 0.7 * matrix(rep(lat, each = 100), 100, 60) + 0.7 * e[1:100, ]
  grp <- rep(c(TRUE, FALSE), each = 30)
  r <- suppressWarnings(
    calibrate_geneset(e, list(a = paste0('g', 1:100)), grp, n_null = 100))
  expect_false(r$reliable)
  expect_gt(r$rho_set, r$rho_null_max)
})

test_that('input validation works', {
  e <- matrix(rnorm(500 * 10), 500, 10, dimnames = list(paste0('g', 1:500), NULL))
  expect_error(calibrate_geneset(e, list(a = paste0('g', 1:20)), rep(TRUE, 5)))
  expect_error(calibrate_geneset(e, list(paste0('g', 1:20)), rep(c(TRUE, FALSE), 5)))
})
