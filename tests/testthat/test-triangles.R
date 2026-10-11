test_that("square plus tail reconciles exactly to event-level data", {
  p <- apply_scenario(params, "heavy_tail_shock")
  ev <- simulate_database(p)
  n_o <- p$horizon$n_origin; n_l <- p$horizon$n_lag
  for (type in names(p$event_types)) {
    for (value in c("count", "amount")) {
      b <- build_square(ev, type, value, n_o, n_l)
      expected <- if (value == "count") sum(ev$type == type) else sum(ev$severity[ev$type == type])
      expect_equal(sum(b$square) + sum(b$tail), expected, tolerance = 1e-9)
    }
  }
})

test_that("true IBNR plus reported equals total per origin", {
  p <- apply_scenario(params, "baseline")
  ev <- simulate_database(p)
  b <- build_square(ev, "nav_error", "amount", p$horizon$n_origin, p$horizon$n_lag)
  t <- true_ibnr(b$square, b$tail)
  expect_equal(t$reported + t$ibnr_total, unname(rowSums(b$square)) + b$tail)
  expect_equal(t$ibnr_total[1], b$tail[1])
  expect_equal(t$ibnr_window[1], 0)
})

test_that("observed triangle has NA exactly in the unobserved region", {
  p <- apply_scenario(params, "baseline")
  b <- build_square(simulate_database(p), "trade_error", "count", p$horizon$n_origin, p$horizon$n_lag)
  tri <- observed_triangle(b$square)
  expect_true(all(is.na(as.matrix(unclass(tri))) == !observed_mask(40, 40)))
  expect_true(inherits(tri, "triangle"))
})

test_that("annual coarsening reconciles and keeps the observed region consistent", {
  p <- apply_scenario(params, "baseline")
  ev <- simulate_database(p)
  b <- build_square(ev, "nav_error", "amount", 40, 40, period = 4)
  expect_equal(dim(b$square), c(10, 10))
  expect_equal(sum(b$square) + sum(b$tail), sum(ev$severity[ev$type == "nav_error"]), tolerance = 1e-9)
  # nothing may be recorded in the unobserved region of what the database has seen
  seen <- ev[ev$origin + ev$lag <= 40, ]
  bs <- build_square(seen, "nav_error", "amount", 40, 40, period = 4)
  expect_equal(sum(bs$square[!observed_mask(10, 10)]), 0)
  expect_equal(sum(bs$tail), 0)
})
