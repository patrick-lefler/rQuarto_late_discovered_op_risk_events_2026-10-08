period_of <- c(annual = 4L, quarterly = 1L)

flatten_estimates <- function(res, kind) {
  purrr::imap(res, \(r, method) {
    n <- length(r$by_origin$ibnr %||% NA)
    tibble::tibble(
      kind = kind, method = method, status = r$status, message = substr(r$message %||% "", 1, 80),
      est = r$total$ibnr %||% NA_real_, se = r$total$se %||% NA_real_,
      lo = r$total$lo %||% NA_real_, hi = r$total$hi %||% NA_real_,
      est_recent = if (is.null(r$by_origin)) NA_real_ else sum(utils::tail(r$by_origin$ibnr, attr(r, "recent_n") %||% 1))
    )
  }) |> purrr::list_rbind()
}

run_one <- function(params, scenario, scenario_index, rep, granularity) {
  period <- period_of[[granularity]]
  p <- apply_scenario(params, scenario)
  events <- simulate_database(p, rep = rep, scenario_index = scenario_index)
  recent_n <- max(1L, p$replication$recent_quarters %/% period)
  seed <- derive_seed(p$master_seed, rep, scenario_index, 99L)
  purrr::map(names(p$event_types), \(type) {
    ic <- make_inputs(events, type, "count", p, period)
    ia <- make_inputs(events, type, "amount", p, period)
    est <- all_estimators(ic, ia, boot_R = p$replication$n_boot, seed = seed)
    est$amount <- purrr::map(est$amount, \(r) structure(r, recent_n = recent_n))
    est$count <- purrr::map(est$count, \(r) structure(r, recent_n = recent_n))
    truth <- \(inp) tibble::tibble(
      truth_window = sum(inp$truth$ibnr_window), truth_total = sum(inp$truth$ibnr_total),
      truth_window_recent = sum(utils::tail(inp$truth$ibnr_window, recent_n)),
      truth_total_recent = sum(utils::tail(inp$truth$ibnr_total, recent_n)),
      reported_recent = sum(utils::tail(inp$truth$reported, recent_n))
    )
    dplyr::bind_rows(
      dplyr::bind_cols(flatten_estimates(est$amount, "amount"), truth(ia)),
      dplyr::bind_cols(flatten_estimates(est$count, "count"), truth(ic))
    ) |> dplyr::mutate(type = type)
  }) |>
    purrr::list_rbind() |>
    dplyr::mutate(scenario = scenario, rep = rep, granularity = granularity, .before = 1)
}

run_replications <- function(params, n_reps = params$replication$n_reps, n_cores = 1L, verbose = TRUE) {
  scen_names <- names(params$scenarios)
  grid <- purrr::imap(params$replication$runs, \(grans, scenario) {
    tidyr::expand_grid(scenario = scenario, granularity = unlist(grans), rep = seq_len(n_reps))
  }) |>
    purrr::list_rbind() |>
    dplyr::mutate(scenario_index = match(scenario, scen_names))
  if (verbose) message(sprintf("%d tasks on %d cores", nrow(grid), n_cores))
  one <- \(i) {
    g <- grid[i, ]
    tryCatch(
      run_one(params, g$scenario, g$scenario_index, g$rep, g$granularity),
      error = \(e) tibble::tibble(scenario = g$scenario, rep = g$rep, granularity = g$granularity, status = "task_error", message = conditionMessage(e))
    )
  }
  out <- parallel::mclapply(seq_len(nrow(grid)), one, mc.cores = n_cores, mc.preschedule = FALSE)
  dplyr::bind_rows(out)
}
