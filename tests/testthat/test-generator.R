test_that("same seed gives identical databases, different seeds differ", {
  p <- apply_scenario(params, "baseline")
  expect_identical(simulate_database(p, rep = 1), simulate_database(p, rep = 1))
  expect_false(identical(simulate_database(p, rep = 1), simulate_database(p, rep = 2)))
})

test_that("events are valid", {
  ev <- simulate_database(apply_scenario(params, "baseline"))
  expect_true(all(ev$lag >= 0 & ev$origin >= 1 & ev$origin <= params$horizon$n_origin))
  expect_true(all(ev$severity > 0))
  expect_setequal(unique(ev$type), names(params$event_types))
})

test_that("fast-lag share is close to its parameter", {
  p <- params$event_types$trade_error$lag
  lag <- draw_lag(2e5, p)
  expect_equal(mean(lag == 0), p$w_fast * p$p_fast + (1 - p$w_fast) * mean(floor(rlnorm(2e5, p$slow_meanlog, p$slow_sdlog)) == 0), tolerance = 0.01)
})

test_that("tail severities exceed the splice point", {
  s <- params$event_types$unauthorized_trading$severity
  x <- draw_severity(1e5, modifyList(s, list(p_tail = 1, lag_elasticity = 0)), lag = 0)
  expect_gte(min(x), qlnorm(s$tail_q, s$meanlog, s$sdlog))
})

test_that("control shift never produces negative lags and shortens discovery", {
  base <- simulate_database(apply_scenario(params, "baseline"))
  shifted <- simulate_database(apply_scenario(params, "control_improvement"))
  expect_true(all(shifted$lag >= 0))
  expect_lt(mean(shifted$lag), mean(base$lag))
})

test_that("injected events appear with exact severity", {
  ev <- simulate_database(apply_scenario(params, "heavy_tail_shock"))
  expect_true(any(ev$severity == 6e7 & ev$origin == 12 & ev$lag == 15))
})

test_that("invalid parameters are rejected", {
  bad <- params
  bad$event_types$nav_error$lag$w_fast <- 1.5
  expect_error(validate_params(bad))
  expect_error(apply_scenario(params, "nope"))
})

test_that("control shift halves the time remaining measured from max(origin, shift quarter)", {
  cq <- params$scenarios$control_improvement$control_shift$calendar_quarter
  ev <- tibble::tibble(origin = rep(c(10, 30), each = 20000), lag = rep(c(30, 6), each = 20000))
  out <- apply_control_shift(ev, list(calendar_quarter = cq, reduction = 0.5))
  # pending at cq: discovery 40, remaining 16 -> 8, so discovery 32 and lag 22 (exact for even remaining time)
  expect_true(all(out$lag[ev$origin == 10] == 22))
  # originates after cq: the whole lag falls by half, not collapses to zero
  expect_equal(mean(out$lag[ev$origin == 30]), 3, tolerance = 0.02)
  # odd remaining time rounds stochastically, preserving the mean reduction
  odd <- apply_control_shift(tibble::tibble(origin = rep(30, 20000), lag = 3), list(calendar_quarter = cq, reduction = 0.5))
  expect_equal(mean(odd$lag), 1.5, tolerance = 0.02)
  # events already discovered before the shift quarter are untouched
  early <- apply_control_shift(tibble::tibble(origin = 2, lag = 5), list(calendar_quarter = cq, reduction = 0.5))
  expect_equal(early$lag, 5)
})

test_that("control improvement leaves a visible share of events undiscovered at the evaluation date", {
  n_o <- params$horizon$n_origin
  und <- \(sc) {
    ev <- purrr::map(1:10, \(r) simulate_database(apply_scenario(params, sc), rep = r, scenario_index = 1) |> dplyr::filter(type == "trade_error"))
    mean(purrr::map_dbl(ev, \(e) mean(e$origin + e$lag > n_o)))
  }
  b <- und("baseline")
  ci <- und("control_improvement")
  expect_gt(b, 0)
  expect_gt(ci, 0.3 * b)
  expect_lt(ci, b)
})

test_that("negative-binomial counts are overdispersed relative to Poisson with the same mean", {
  pn <- apply_scenario(params, "overdispersion")
  expect_identical(pn$count_dist, "negbin")
  counts <- \(p) purrr::map_dbl(1:150, \(r) sum(simulate_database(p, rep = r)$type == "trade_error"))
  cb <- counts(apply_scenario(params, "baseline"))
  cn <- counts(pn)
  expect_equal(mean(cn), mean(cb), tolerance = 0.05)
  expect_gt(var(cn), 2 * var(cb))
})

test_that("reporting threshold removes small losses and nothing else", {
  p <- apply_scenario(params, "reporting_threshold_10k")
  ev <- simulate_database(p, rep = 1)
  full <- simulate_database(apply_scenario(params, "baseline"), rep = 1)
  expect_true(all(ev$severity >= 1e4))
  expect_equal(nrow(ev), sum(full$severity >= 1e4))
})
