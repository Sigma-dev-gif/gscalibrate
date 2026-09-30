
test_that("cox model runs and returns a coefficient", {
  set.seed(7)
  e <- matrix(rnorm(2000 * 80), 2000, 80,
              dimnames = list(paste0("g", 1:2000), paste0("s", 1:80)))
  tm <- rexp(80, 0.05); ev <- rbinom(80, 1, 0.6)
  r <- calibrate_geneset(e, list(a = paste0("g", 1:60)), model = "cox",
                         time = tm, event = ev, n_null = 50)
  expect_equal(r$model, "cox")
  expect_true(is.finite(r$beta))
  expect_true(r$p_empirical > 0 && r$p_empirical <= 1)
})
