fake <- function(est, lo, hi, status = "ok", truth = 100) {
  tibble::tibble(
    scenario = "s", granularity = "annual", type = "t", kind = "amount", method = "m",
    status = status, est = est, lo = lo, hi = hi, est_recent = est,
    truth_window = truth, truth_total = truth * 2, truth_window_recent = truth, truth_total_recent = truth * 2,
    reported_recent = 50
  )
}

test_that("scoring counts coverage, bias and unusable results correctly", {
  df <- dplyr::bind_rows(
    fake(110, 90, 130),
    fake(80, 60, 90),
    fake(120, 105, 140),
    fake(NA_real_, NA_real_, NA_real_, status = "failed"),
    fake(0, 0, 0, status = "suspect")
  )
  s <- score_estimates(df)
  expect_equal(s$n, 5)
  expect_equal(s$n_failed, 1)
  expect_equal(s$n_suspect, 1)
  expect_equal(s$usable_rate, 3 / 5)
  expect_equal(s$n_usable, 3)
  expect_equal(s$rel_bias, (10 - 20 + 20) / 300)
  expect_equal(s$coverage_usable, 1 / 3)
  expect_equal(s$coverage_all, 1 / 5)
})

test_that("tail variants are scored against total, others against the window", {
  a <- fake(100, 80, 120) |> dplyr::mutate(method = "clark_ldf_loglogistic_tail")
  b <- fake(100, 80, 120)
  expect_equal(score_estimates(a)$rel_bias, -0.5)
  expect_equal(score_estimates(b)$rel_bias, 0)
})

test_that("a method with no intervals has NA coverage, not zero", {
  expect_true(is.na(score_estimates(fake(100, NA_real_, NA_real_))$coverage_usable))
})

test_that("apply_control shifts the lag distribution and clamps w_fast", {
  p <- apply_scenario(params, "baseline")
  c1 <- apply_control(p, "independent_pnl_verification")
  expect_equal(c1$event_types$unauthorized_trading$lag$slow_meanlog, p$event_types$unauthorized_trading$lag$slow_meanlog - 0.8)
  expect_lte(c1$event_types$unauthorized_trading$lag$w_fast, 1)
  p$controls$recon_automation$lag_shift$nav_error$w_fast <- 5
  expect_equal(apply_control(p, "recon_automation")$event_types$nav_error$lag$w_fast, 1)
  expect_error(apply_control(p, "nope"))
})

test_that("a faster-detection control lowers undiscovered loss on paired seeds", {
  p <- apply_scenario(params, "baseline")
  base <- control_metrics(simulate_database(p, rep = 3), 40)$summary
  ctl <- control_metrics(simulate_database(apply_control(p, "independent_pnl_verification"), rep = 3), 40)$summary
  expect_lt(sum(ctl$undiscovered), sum(base$undiscovered))
})

test_that("replications are identical regardless of the number of cores", {
  p <- params
  p$replication$n_boot <- 50
  p$replication$runs <- list(baseline = list("annual"))
  one <- run_replications(p, n_reps = 2, n_cores = 1, verbose = FALSE)
  two <- run_replications(p, n_reps = 2, n_cores = 2, verbose = FALSE)
  expect_equal(one, two)
  expect_false("task_error" %in% one$status)
})

test_that("incurred_gap reports the share of recent incurred loss that is unreported", {
  df <- fake(100, 80, 120) |> dplyr::mutate(method = "chainladder")
  g <- incurred_gap(df)
  expect_equal(g$mean_understatement, 200 / (50 + 200))
})

test_that("apply_control scale multiplies the shifts; scale 0 changes nothing", {
  p <- apply_scenario(params, "baseline")
  half <- apply_control(p, "recon_automation", scale = 0.5)
  expect_equal(half$event_types$unauthorized_trading$lag$slow_meanlog, p$event_types$unauthorized_trading$lag$slow_meanlog - 0.2)
  expect_equal(apply_control(p, "recon_automation", scale = 0)$event_types, p$event_types)
})

test_that("break_even_scale interpolates the zero crossing of net benefit", {
  curve <- tibble::tibble(control = "c", scale = c(0, 1, 2), net = c(-100, -20, 60))
  expect_equal(break_even_scale(curve)$break_even_scale, 1 + 20 / 80)
  expect_true(is.na(break_even_scale(dplyr::mutate(curve, net = -1))$break_even_scale))
})
