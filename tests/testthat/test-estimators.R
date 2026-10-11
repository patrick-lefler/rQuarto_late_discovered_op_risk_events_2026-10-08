test_that("project_chainladder reproduces the ChainLadder package on a dense triangle", {
  tri <- ChainLadder::RAA
  mine <- project_chainladder(tri)
  pkg <- stats::predict(ChainLadder::chainladder(tri, weights = 1))
  expect_equal(unname(mine[, 10]), unname(pkg[, 10]), tolerance = 1e-8)
})

test_that("project_chainladder uses factor 1 where a development column has zero volume", {
  m <- matrix(c(5, 8, 8, 8,
                4, 6, 6, NA,
                7, 9, NA, NA,
                3, NA, NA, NA), 4, byrow = TRUE)
  full <- project_chainladder(ChainLadder::as.triangle(m))
  expect_equal(full[4, 4], 3 * 23 / 16)

  z <- matrix(c(0, 2, 2,
                0, 3, NA,
                0, NA, NA), 3, byrow = TRUE)
  expect_true(all(is.finite(project_chainladder(ChainLadder::as.triangle(z)))))
})

test_that("run_safely converts errors into a failed result and never throws", {
  r <- run_safely("boom", \() stop("bad triangle"))
  expect_equal(r$status, "failed")
  expect_match(r$message, "bad triangle")
})

test_that("estimators run on a dense triangle and are close to the truth on average", {
  # one database is dominated by Poisson noise (~40 undiscovered events), so average several
  p <- apply_scenario(params, "baseline")
  runs <- purrr::map(1:6, \(rep) {
    ev <- simulate_database(p, rep = rep)
    ic <- make_inputs(ev, "trade_error", "count", p, period = 4)
    ia <- make_inputs(ev, "trade_error", "amount", p, period = 4)
    list(ic = ic, ia = ia, res = list(est_chainladder(ic), est_mack(ic), est_boot(ic, R = 200)), fs = est_freqsev(ic, ia))
  })
  truth <- mean(purrr::map_dbl(runs, \(r) sum(r$ic$truth$ibnr_window)))
  for (k in 1:3) {
    expect_true(all(purrr::map_chr(runs, \(r) r$res[[k]]$status) == "ok"))
    expect_lt(abs(mean(purrr::map_dbl(runs, \(r) r$res[[k]]$total$ibnr)) / truth - 1), 0.3)
  }
  expect_true(all(purrr::map_chr(runs, \(r) r$fs$status) == "ok"))
  truth_amt <- mean(purrr::map_dbl(runs, \(r) sum(r$ia$truth$ibnr_window)))
  expect_lt(abs(mean(purrr::map_dbl(runs, \(r) r$fs$total$ibnr)) / truth_amt - 1), 0.5)
})

test_that("sparse triangles are reported as failed or suspect, never as silent zeros", {
  p <- apply_scenario(params, "baseline")
  ev <- simulate_database(p)
  inp <- make_inputs(ev, "unauthorized_trading", "count", p, period = 1)
  expect_equal(est_mack(inp)$status, "failed")
  expect_true(est_boot(inp, R = 100)$status %in% c("suspect", "failed"))
})

test_that("lognormal interval brackets the mean and is NA for degenerate input", {
  ci <- interval_lognormal(100, 30)
  expect_lt(ci[1], 100)
  expect_gt(ci[2], 100)
  expect_true(all(is.na(interval_lognormal(0, 1))))
})

test_that("make_inputs never leaks unobserved events into the triangle", {
  p <- apply_scenario(params, "baseline")
  ev <- simulate_database(p)
  inp <- make_inputs(ev, "nav_error", "amount", p, period = 4)
  m <- unclass(inp$tri_incr)
  expect_true(all(is.na(m[!observed_mask(10, 10)])))
  expect_equal(sum(m, na.rm = TRUE), sum(ev$severity[ev$type == "nav_error" & ev$origin + ev$lag <= 40]), tolerance = 1e-9)
})
