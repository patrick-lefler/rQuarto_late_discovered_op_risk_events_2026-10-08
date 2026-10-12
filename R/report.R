# Helpers shared by the Quarto pages (labels, theme, loaders, plots).
type_lab <- c(trade_error = "Trade errors", nav_error = "NAV errors", unauthorized_trading = "Unauthorized trading")
type_col <- c("Trade errors" = "#0072B2", "NAV errors" = "#E69F00", "Unauthorized trading" = "#D55E00")
method_lab <- c(
  chainladder = "Chain-ladder (vol-weighted)", mack = "Mack (package)", mack_ec = "Mack (empty columns = 1)",
  mack_last3 = "Mack (latest 3 diagonals)", bf = "Bornhuetter-Ferguson", benktander = "Benktander",
  bf_benchmark = "BF (benchmark pattern)", benktander_benchmark = "Benktander (benchmark pattern)", boot = "Bootstrap",
  clark_ldf_loglogistic = "Clark LDF loglogistic", clark_ldf_weibull = "Clark LDF Weibull",
  clark_ldf_loglogistic_tail = "Clark LDF loglogistic + tail", clark_cc_weibull = "Clark Cape Cod Weibull",
  freq_x_sev = "Frequency x severity", freq_x_sev_uplift = "Frequency x severity (pooled, uplift)"
)
method_order <- unname(method_lab)
ctl_lab <- c(baseline = "Baseline", recon_automation = "Reconciliation automation", independent_pnl_verification = "Independent P&L verification")
ctl_col <- c(recon_automation = "#0072B2", independent_pnl_verification = "#E69F00")

fmt_usd <- function(x, digits = 1) {
  dplyr::case_when(
    is.na(x) ~ "NA",
    abs(x) >= 1e9 ~ sprintf("$%.*fB", digits, x / 1e9),
    abs(x) >= 1e6 ~ sprintf("$%.*fM", digits, x / 1e6),
    abs(x) >= 1e3 ~ sprintf("$%.*fK", digits, x / 1e3),
    TRUE ~ sprintf("$%.0f", x)
  )
}
fmt_pct <- function(x, digits = 0) ifelse(is.na(x), "NA", sprintf("%.*f%%", digits, 100 * x))

theme_op <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title.position = "plot",
      plot.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )
}

load_results <- function(dir = "data") {
  list(
    reps = readRDS(file.path(dir, "replications.rds")),
    ctl = readRDS(file.path(dir, "control_study.rds")),
    be = readRDS(file.path(dir, "control_breakeven.rds")),
    ctl_sens = readRDS(file.path(dir, "control_sensitivity.rds")),
    ctl_tail = readRDS(file.path(dir, "control_tail_ci.rds")),
    ctl_tail_diff = readRDS(file.path(dir, "control_tail_diff.rds")),
    meta = readRDS(file.path(dir, "run_meta.rds"))
  )
}

label_scores <- function(s) {
  s |>
    dplyr::mutate(type = factor(type_lab[type], levels = type_lab), method = factor(method_lab[method], levels = method_order))
}

# Observed vs true triangle: observed cells coloured by size, unobserved (not yet discovered) cells hatched grey with true values faintly shown.
plot_triangle <- function(square, title, fill_lab = "Value") {
  n_o <- nrow(square)
  n_l <- ncol(square)
  obs <- observed_mask(n_o, n_l)
  d <- tibble::tibble(
    origin = rep(seq_len(n_o), n_l), lag = rep(seq_len(n_l) - 1, each = n_o),
    value = as.vector(square), observed = as.vector(obs)
  )
  ggplot2::ggplot(d, ggplot2::aes(lag, origin)) +
    ggplot2::geom_tile(data = ~ dplyr::filter(.x, observed), ggplot2::aes(fill = value), colour = "white") +
    ggplot2::geom_tile(data = ~ dplyr::filter(.x, !observed), ggplot2::aes(alpha = value), fill = "#D55E00", colour = "white") +
    ggplot2::scale_fill_gradient(low = "#EAF2F8", high = "#0072B2", name = fill_lab, labels = scales::label_number(scale_cut = scales::cut_short_scale()),
                                 guide = ggplot2::guide_colourbar(barwidth = 12)) +
    ggplot2::scale_alpha(range = c(0.1, 0.9), guide = "none") +
    ggplot2::scale_y_reverse(breaks = seq_len(n_o)) +
    ggplot2::scale_x_continuous(breaks = seq_len(n_l) - 1) +
    ggplot2::labs(title = title, subtitle = "Blue: recorded to date. Orange: true but not yet discovered (invisible to the database).", x = "Discovery lag (periods after origin)", y = "Origin period") +
    theme_op() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}

# Methods with fewer usable fits than this are shown as "unusable" in scorecards.
min_usable <- 20

plot_scorecard <- function(s, kind_, scen = "baseline", gran = "annual") {
  d <- s |>
    dplyr::filter(kind == kind_, scenario == scen, granularity == gran) |>
    label_scores() |>
    dplyr::mutate(
      shown = n_usable >= min_usable,
      fill_val = dplyr::if_else(shown, pmin(abs(rel_bias), 1), NA_real_),
      label = dplyr::if_else(shown, sprintf("%+.0f%%\nn=%d", 100 * rel_bias, n_usable), sprintf("unusable\nn=%d", n_usable))
    )
  ggplot2::ggplot(d, ggplot2::aes(type, forcats::fct_rev(method), fill = fill_val)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 2.8, lineheight = 0.9) +
    ggplot2::scale_fill_gradient(low = "#EAF2F8", high = "#D55E00", na.value = "grey85", limits = c(0, 1), name = "|relative bias| (capped at 100%)", labels = scales::label_percent()) +
    ggplot2::labs(x = NULL, y = NULL, title = "Relative bias of the total undiscovered estimate", subtitle = "n = replications (of 200) with a usable result; grey = unusable") +
    theme_op()
}

plot_coverage <- function(s, kind_, scen = "baseline", gran = "annual") {
  d <- s |>
    dplyr::filter(kind == kind_, scenario == scen, granularity == gran, has_interval, n_usable >= min_usable) |>
    label_scores()
  ggplot2::ggplot(d, ggplot2::aes(coverage_usable, forcats::fct_rev(method), colour = type)) +
    ggplot2::geom_vline(xintercept = 0.95, linetype = 2, colour = "grey40") +
    ggplot2::geom_point(size = 3, position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::scale_colour_manual(values = type_col, name = NULL) +
    ggplot2::scale_x_continuous(limits = c(0, 1), labels = scales::label_percent()) +
    ggplot2::labs(x = "Share of usable fits whose 95% interval contained the truth", y = NULL, title = "Interval coverage (dashed line = nominal 95%)") +
    theme_op()
}

plot_errors <- function(reps, kind_, scen = "baseline", gran = "annual", methods = NULL) {
  d <- reps |>
    add_target() |>
    dplyr::filter(kind == kind_, scenario == scen, granularity == gran, status == "ok", is.finite(est)) |>
    dplyr::mutate(err = est - truth)
  if (!is.null(methods)) d <- dplyr::filter(d, method %in% methods)
  d <- d |>
    dplyr::mutate(cap = 3 * mean(truth), .by = c(type, kind)) |>
    dplyr::mutate(clipped = abs(err) > cap, err = pmax(pmin(err, cap), -cap)) |>
    label_scores()
  ggplot2::ggplot(d, ggplot2::aes(err, forcats::fct_rev(method), fill = type)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey40") +
    ggplot2::geom_boxplot(outlier.size = 0.6, linewidth = 0.3) +
    ggplot2::facet_wrap(~type, scales = "free_x", nrow = 1) +
    ggplot2::scale_fill_manual(values = type_col, guide = "none") +
    ggplot2::scale_x_continuous(labels = if (kind_ == "amount") fmt_usd else scales::label_number()) +
    ggplot2::labs(x = "Estimate minus truth (clipped at +/- 3x the mean true value)", y = NULL, title = "Estimation error by method (usable fits only)") +
    theme_op()
}

score_gt <- function(s, kind_, scen = "baseline", gran = "annual") {
  s |>
    dplyr::filter(kind == kind_, scenario == scen, granularity == gran) |>
    label_scores() |>
    dplyr::arrange(type, method) |>
    dplyr::mutate(
      # fewer than min_usable fits is not interpretable: show a dash rather than a bias, coverage or error
      dplyr::across(c(rel_bias, median_rel_err, coverage_usable, rel_bias_recent), \(x) dplyr::if_else(n_usable >= min_usable, x, NA_real_)),
      coverage_all = dplyr::if_else(has_interval, coverage_all, NA_real_),
      median_rel_err = dplyr::if_else(is.finite(median_rel_err) & abs(median_rel_err) < 0.005, 0, median_rel_err)
    ) |>
    dplyr::transmute(
      Type = type, Method = method, `Usable (of 200)` = n_usable, Failed = n_failed, Suspect = n_suspect,
      `Rel. bias` = rel_bias, `Median signed rel. error` = median_rel_err,
      `Coverage (usable fits)` = coverage_usable, `Coverage (all fits)` = coverage_all,
      `Recent-origin bias` = rel_bias_recent
    ) |>
    gt::gt(groupname_col = "Type") |>
    gt::fmt_percent(c(`Rel. bias`, `Median signed rel. error`, `Coverage (usable fits)`, `Coverage (all fits)`, `Recent-origin bias`), decimals = 0) |>
    gt::sub_missing(missing_text = "-") |>
    gt::tab_source_note(paste(
      "Relative bias = sum(estimate - truth) / sum(truth) over usable fits.",
      "Coverage (usable fits) = share of usable fits whose 95% interval contains the truth; Coverage (all fits) counts failed, suspect and interval-free fits as misses.",
      "Recent-origin bias = relative bias restricted to the most recent origin period (latest 4 quarters).",
      "Failed or suspect fits are excluded from bias. A dash means fewer than 20 usable fits (not interpretable) or no interval."
    ))
}

# Share of the true undiscovered loss that lies beyond the development window (discovery lag >= n_lag),
# pooled over replications. Uses the chain-ladder rows only because they carry one truth row per database.
tail_share <- function(reps, type_, scen = "baseline", gran = "annual") {
  d <- reps |>
    dplyr::filter(method == "chainladder", kind == "amount", scenario == scen, granularity == gran, type == type_)
  1 - sum(d$truth_window) / sum(d$truth_total)
}

# Mean estimate vs. mean truth of undiscovered dollar loss (usable fits only), by scenario, type and method.
stress_means <- function(reps, methods = c("chainladder", "mack", "mack_last3", "bf_benchmark", "boot", "freq_x_sev", "freq_x_sev_uplift")) {
  reps |>
    add_target() |>
    dplyr::filter(kind == "amount", granularity == "annual", status == "ok", is.finite(est), method %in% methods) |>
    dplyr::summarise(Estimate = mean(est), Truth = mean(truth), n = dplyr::n(), .by = c(scenario, type, method)) |>
    dplyr::mutate(
      scenario = factor(scenario, levels = c("baseline", "heavy_tail_shock", "control_improvement", "exposure_growth", "overdispersion", "reporting_threshold_10k")),
      type = factor(type_lab[type], levels = type_lab),
      method = factor(method_lab[method], levels = method_order)
    )
}

# Share of events whose discovery lag falls at or beyond the window, from a large fixed-seed draw of the lag
# mixture (so the Data and Overview pages quote the same number).
beyond_window_share <- function(params, n = 2e5, seed = 1) {
  set.seed(seed)
  purrr::map_dbl(params$event_types, \(spec) mean(draw_lag(n, spec$lag) >= params$horizon$n_lag))
}
