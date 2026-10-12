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

test_that("mack_manual reproduces ChainLadder::MackChainLadder on a dense triangle", {
  tri <- ChainLadder::RAA
  pkg <- ChainLadder::MackChainLadder(tri, est.sigma = "Mack")
  mine <- mack_manual(tri)
  expect_equal(unname(mine$ibnr[-1]), unname(summary(pkg)$ByOrigin$IBNR[-1]), tolerance = 1e-6)
  expect_equal(mine$total_ibnr, unname(summary(pkg)$Totals["IBNR:", 1]), tolerance = 1e-6)
  expect_equal(mine$total_se, unname(summary(pkg)$Totals["Mack S.E.:", 1]), tolerance = 1e-4)
})

test_that("mack_manual runs where the package fails (empty development columns) and reports factor 1 there", {
  m <- matrix(c(0, 0, 0, 5, 5,
                0, 0, 0, 0, NA,
                0, 0, 4, NA, NA,
                0, 0, NA, NA, NA,
                0, NA, NA, NA, NA), 5, byrow = TRUE)
  tri <- ChainLadder::as.triangle(m)
  expect_error(suppressWarnings(ChainLadder::MackChainLadder(tri, est.sigma = "Mack")))
  r <- mack_manual(tri)
  expect_true(all(is.finite(r$ibnr)) && is.finite(r$total_se))
  expect_equal(r$factors[1], 1)
  expect_equal(est_mack_manual(list(tri_cum = tri))$status, "ok")
})

test_that("restricting to the latest diagonals changes factors only when development drifts", {
  tri <- ChainLadder::RAA
  expect_equal(mack_manual(tri, n_diag = Inf)$factors, mack_manual(tri)$factors)
  expect_false(isTRUE(all.equal(mack_manual(tri, n_diag = 3)$factors, mack_manual(tri)$factors)))
})

test_that("benchmark unreported share is a proper share, rises for recent origins and matches simulation", {
  p <- apply_scenario(params, "baseline")
  b <- benchmark_unreported(p, "nav_error", "amount", period = 4)
  expect_true(all(b >= 0 & b <= 1))
  expect_true(all(diff(b) >= -1e-12))
  expect_lt(b[1], 0.01)
  sim <- purrr::map(1:40, \(r) {
    ev <- simulate_database(p, rep = r)
    ia <- make_inputs(ev, "nav_error", "amount", p, period = 4)
    c(ibnr = sum(ia$truth$ibnr_window), tot = sum(ia$truth$ibnr_window + ia$truth$reported))
  })
  frac <- sum(purrr::map_dbl(sim, "ibnr")) / sum(purrr::map_dbl(sim, "tot"))
  e <- make_inputs(simulate_database(p, rep = 1), "nav_error", "amount", p, period = 4)$exposure
  expect_equal(sum(b * e) / sum(e), frac, tolerance = 0.25)
})

test_that("BF equals chain-ladder when the prior is the chain-ladder ultimate, and benktander lies between BF and CL", {
  p <- apply_scenario(params, "baseline")
  ia <- make_inputs(simulate_database(p, rep = 2), "nav_error", "amount", p, period = 4)
  cl <- est_chainladder(ia)$total$ibnr
  bf <- est_bf(ia, "own", "bf")
  gb <- est_bf(ia, "own", "benktander")
  expect_equal(bf$status, "ok")
  expect_true(all(is.finite(c(bf$total$ibnr, gb$total$ibnr))))
  lo <- min(bf$total$ibnr, cl)
  hi <- max(bf$total$ibnr, cl)
  expect_gte(gb$total$ibnr, lo - 1e-6)
  expect_lte(gb$total$ibnr, hi + 1e-6)
  expect_equal(est_bf(ia, "benchmark", "bf")$status, "ok")
})

test_that("a zero bootstrap estimate is flagged suspect, never presented as a result", {
  tri <- ChainLadder::as.triangle(matrix(c(0, 0, 0, 0,
                                           0, 0, 0, NA,
                                           0, 0, NA, NA,
                                           0, NA, NA, NA), 4, byrow = TRUE))
  inp <- list(tri_cum = tri)
  r <- est_boot(inp, R = 50)
  expect_true(r$status %in% c("suspect", "failed"))
})

test_that("uplift-pooled severity gives every lag the same donor pool and a non-negative uplift", {
  p <- apply_scenario(params, "baseline")
  ev <- simulate_database(p, rep = 1)
  ia <- make_inputs(ev, "nav_error", "amount", p, period = 4)
  sp <- severity_pools(ia$events_obs, ia$n_lag, pooled = "uplift")
  expect_true(sp$e >= 0)
  expect_true(all(sp$n_donors == nrow(ia$events_obs)))
  nb <- severity_pools(ia$events_obs, ia$n_lag, pooled = "nearby")
  expect_length(nb$pools, ia$n_lag)
})

test_that("calendar_z is near zero on a stable triangle and negative when the latest diagonals under-develop", {
  set.seed(1)
  n <- 10
  f <- c(2, 1.4, 1.2, 1.1, 1.05, 1.03, 1.02, 1.01, 1.005)
  gen <- \(shock) {
    out <- matrix(NA_real_, n, n)
    for (i in 1:n) {
      out[i, 1] <- 1000 * exp(rnorm(1, 0, 0.1))
      for (k in 1:(n - 1)) {
        fk <- if (i + k >= n) 1 + (f[k] - 1) * shock else f[k]
        out[i, k + 1] <- out[i, k] * fk * exp(rnorm(1, 0, 0.02))
      }
    }
    out[outer(1:n, 1:n, `+`) > n + 1] <- NA
    ChainLadder::as.triangle(out)
  }
  z0 <- mean(replicate(40, calendar_z(gen(1))), na.rm = TRUE)
  z1 <- mean(replicate(40, calendar_z(gen(0.3))), na.rm = TRUE)
  expect_lt(abs(z0), 1)
  expect_lt(z1, z0 - 1)
})
