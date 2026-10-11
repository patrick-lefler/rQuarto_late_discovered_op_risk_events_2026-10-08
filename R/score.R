# Tail variants extrapolate beyond the window, so they are scored against total undiscovered loss;
# everything else is scored against the undiscovered loss inside the development window.
add_target <- function(df) {
  df |>
    dplyr::mutate(
      target = dplyr::if_else(grepl("_tail$", method), "total", "window"),
      truth = dplyr::if_else(target == "total", truth_total, truth_window),
      truth_recent = dplyr::if_else(target == "total", truth_total_recent, truth_window_recent)
    )
}

# status "ok" results are usable; "suspect" and "failed" count as unusable (a zero presented as an estimate is worse than an error).
score_estimates <- function(df) {
  df |>
    add_target() |>
    dplyr::mutate(usable = status == "ok" & is.finite(est), covered = usable & is.finite(lo) & is.finite(hi) & lo <= truth & truth <= hi) |>
    dplyr::summarise(
      n = dplyr::n(),
      n_failed = sum(status == "failed"),
      n_suspect = sum(status == "suspect"),
      n_usable = sum(usable),
      usable_rate = mean(usable),
      rel_bias = if (sum(truth[usable]) > 0) sum(est[usable] - truth[usable]) / sum(truth[usable]) else NA_real_,
      median_rel_err = stats::median((est[usable] - truth[usable]) / pmax(truth[usable], 1e-9), na.rm = TRUE),
      rmse = if (any(usable)) sqrt(mean((est[usable] - truth[usable])^2)) else NA_real_,
      coverage_usable = if (any(usable & is.finite(lo))) mean(covered[usable & is.finite(lo)]) else NA_real_,
      coverage_all = mean(covered),
      has_interval = any(usable & is.finite(lo)),
      rel_bias_recent = if (sum(truth_recent[usable]) > 0) sum(est_recent[usable] - truth_recent[usable]) / sum(truth_recent[usable]) else NA_real_,
      .by = c(scenario, granularity, type, kind, method, target)
    )
}

# How much of the latest-origin incurred loss is missing from the reported figure?
incurred_gap <- function(df) {
  df |>
    dplyr::filter(method == "chainladder", kind == "amount") |>
    dplyr::mutate(understatement = truth_total_recent / (reported_recent + truth_total_recent)) |>
    dplyr::summarise(
      mean_understatement = mean(understatement, na.rm = TRUE),
      median_understatement = stats::median(understatement, na.rm = TRUE),
      p90_understatement = stats::quantile(understatement, 0.9, na.rm = TRUE, names = FALSE),
      .by = c(scenario, granularity, type)
    )
}
