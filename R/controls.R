# True (known) loss under a control: incurred = every event originating in the window,
# reported = those discovered by the evaluation date, undiscovered = the difference.
control_metrics <- function(events, n_origin) {
  by_year <- events |>
    dplyr::mutate(year = ceiling(origin / 4)) |>
    dplyr::summarise(incurred = sum(severity), .by = year)
  list(
    summary = events |>
      dplyr::summarise(
        incurred = sum(severity),
        reported = sum(severity[origin + lag <= n_origin]),
        events = dplyr::n(),
        .by = type
      ) |>
      dplyr::mutate(undiscovered = incurred - reported),
    annual = by_year
  )
}

run_control_study <- function(params, n_reps = params$replication$n_reps, n_cores = 1L) {
  arms <- c("baseline", names(params$controls))
  grid <- tidyr::expand_grid(arm = arms, rep = seq_len(n_reps))
  one <- \(i) {
    g <- grid[i, ]
    p <- apply_scenario(params, "baseline")
    if (g$arm != "baseline") p <- apply_control(p, g$arm)
    m <- control_metrics(simulate_database(p, rep = g$rep, scenario_index = 1L), p$horizon$n_origin)
    list(
      summary = dplyr::mutate(m$summary, arm = g$arm, rep = g$rep),
      annual = dplyr::mutate(m$annual, arm = g$arm, rep = g$rep)
    )
  }
  out <- parallel::mclapply(seq_len(nrow(grid)), one, mc.cores = n_cores)
  list(
    summary = purrr::map(out, "summary") |> purrr::list_rbind(),
    annual = purrr::map(out, "annual") |> purrr::list_rbind()
  )
}

# Paired comparison against the baseline arm (same seeds, so counts and severity draws are shared).
summarise_controls <- function(study, params) {
  n_years <- params$horizon$n_origin / 4
  per_rep <- study$summary |>
    dplyr::summarise(incurred = sum(incurred), undiscovered = sum(undiscovered), .by = c(arm, rep)) |>
    dplyr::mutate(incurred_per_year = incurred / n_years)
  base <- per_rep |> dplyr::filter(arm == "baseline") |> dplyr::select(rep, base_inc = incurred_per_year, base_und = undiscovered)
  cost <- purrr::map_dbl(params$controls, "annual_cost_usd")
  tail_by_arm <- study$annual |>
    dplyr::summarise(
      annual_var95 = stats::quantile(incurred, 0.95, names = FALSE),
      annual_es95 = mean(incurred[incurred >= stats::quantile(incurred, 0.95, names = FALSE)]),
      .by = arm
    )
  per_rep |>
    dplyr::left_join(base, by = "rep") |>
    dplyr::mutate(saving_per_year = base_inc - incurred_per_year, undiscovered_reduction = base_und - undiscovered) |>
    dplyr::summarise(
      mean_incurred_per_year = mean(incurred_per_year),
      mean_saving_per_year = mean(saving_per_year),
      saving_se = stats::sd(saving_per_year) / sqrt(dplyr::n()),
      p_saving_positive = mean(saving_per_year > 0),
      mean_undiscovered = mean(undiscovered),
      mean_undiscovered_reduction = mean(undiscovered_reduction),
      .by = arm
    ) |>
    dplyr::left_join(tail_by_arm, by = "arm") |>
    dplyr::mutate(annual_cost = unname(cost[arm]), net_benefit_per_year = mean_saving_per_year - annual_cost)
}

# Break-even effectiveness: scale each control's lag shifts by a factor and find where the
# paired mean annual saving equals the control's annual cost (linear interpolation on the grid).
run_breakeven <- function(params, scales = c(0, 0.25, 0.5, 0.75, 1, 1.5, 2), n_reps = params$replication$n_reps, n_cores = 1L) {
  n_years <- params$horizon$n_origin / 4
  grid <- tidyr::expand_grid(control = names(params$controls), scale = scales, rep = seq_len(n_reps))
  base_inc <- unlist(parallel::mclapply(seq_len(n_reps), \(r) {
    sum(simulate_database(apply_scenario(params, "baseline"), rep = r)$severity) / n_years
  }, mc.cores = n_cores))
  inc <- unlist(parallel::mclapply(seq_len(nrow(grid)), \(i) {
    g <- grid[i, ]
    p <- apply_control(apply_scenario(params, "baseline"), g$control, g$scale)
    sum(simulate_database(p, rep = g$rep)$severity) / n_years
  }, mc.cores = n_cores))
  grid$saving <- base_inc[grid$rep] - inc
  curve <- grid |>
    dplyr::summarise(mean_saving = mean(saving), se = stats::sd(saving) / sqrt(dplyr::n()), .by = c(control, scale)) |>
    dplyr::mutate(annual_cost = purrr::map_dbl(params$controls, "annual_cost_usd")[control], net = mean_saving - annual_cost)
  list(curve = curve, break_even = break_even_scale(curve))
}

break_even_scale <- function(curve) {
  curve |>
    dplyr::arrange(control, scale) |>
    dplyr::summarise(
      break_even_scale = {
        above <- which(net >= 0)
        if (length(above) == 0) NA_real_ else if (above[1] == 1) scale[1] else {
          i <- above[1]
          scale[i - 1] + (0 - net[i - 1]) * (scale[i] - scale[i - 1]) / (net[i] - net[i - 1])
        }
      },
      .by = control
    )
}
