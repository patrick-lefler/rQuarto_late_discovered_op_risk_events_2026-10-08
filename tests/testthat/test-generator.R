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
