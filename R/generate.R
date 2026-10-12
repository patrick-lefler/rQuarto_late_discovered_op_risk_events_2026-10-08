draw_lag <- function(n, lag) {
  fast <- stats::runif(n) < lag$w_fast
  fast_lag <- stats::rgeom(n, lag$p_fast)
  slow_lag <- floor(stats::rlnorm(n, lag$slow_meanlog, lag$slow_sdlog))
  pmin(ifelse(fast, fast_lag, slow_lag), 1e4)
}

draw_severity <- function(n, sev, lag) {
  u <- stats::qlnorm(sev$tail_q, sev$meanlog, sev$sdlog)
  is_tail <- stats::runif(n) < sev$p_tail
  body <- stats::qlnorm(stats::runif(n) * sev$tail_q, sev$meanlog, sev$sdlog)
  tail <- u + sev$beta * ((1 - stats::runif(n))^(-sev$xi) - 1) / sev$xi
  ifelse(is_tail, tail, body) * (1 + sev$lag_elasticity * lag)
}

# "Remaining time to discovery" is measured from the later of the shift quarter and the event's origin, so an event
# that originates after the shift has its whole lag reduced, never a lag that was partly spent before it existed.
# Fractional quarters are rounded stochastically so the mean reduction is exactly `reduction`.
apply_control_shift <- function(events, shift) {
  if (is.null(shift) || nrow(events) == 0) return(events)
  start <- pmax(events$origin, shift$calendar_quarter)
  discovery <- events$origin + events$lag
  remaining <- pmax(discovery - start, 0) * (1 - shift$reduction)
  rounded <- floor(remaining) + (stats::runif(nrow(events)) < remaining - floor(remaining))
  events$lag <- ifelse(discovery > start, start + rounded - events$origin, events$lag)
  events
}

draw_counts <- function(n, mu, spec, params) {
  if (identical(params$count_dist, "negbin")) {
    stats::rnbinom(n, size = params$count_size[[spec$name]], mu = mu)
  } else {
    stats::rpois(n, mu)
  }
}

simulate_type <- function(type, params, seed) {
  spec <- params$event_types[[type]]
  set.seed(seed)
  n_origin <- params$horizon$n_origin
  exposure <- params$exposure$base * (1 + params$exposure$growth_per_quarter)^(seq_len(n_origin) - 1)
  spec$name <- type
  counts <- draw_counts(n_origin, spec$lambda * exposure, spec, params)
  origin <- rep(seq_len(n_origin), counts)
  lag <- draw_lag(length(origin), spec$lag)
  severity <- draw_severity(length(origin), spec$severity, lag)
  tibble::tibble(type = type, origin = origin, lag = lag, severity = severity) |>
    apply_control_shift(params$control_shift)
}

inject_events <- function(events, injected) {
  if (length(injected) == 0) return(events)
  extra <- purrr::map(injected, tibble::as_tibble) |> purrr::list_rbind()
  dplyr::bind_rows(events, dplyr::mutate(extra, dplyr::across(c(origin, lag), as.numeric)))
}

simulate_database <- function(params, rep = 1L, scenario_index = 1L) {
  types <- names(params$event_types)
  events <- purrr::imap(types, \(type, i) {
    simulate_type(type, params, derive_seed(params$master_seed, rep, scenario_index, i))
  }) |>
    purrr::list_rbind() |>
    inject_events(params$injected_events)
  events <- events[events$severity >= (params$reporting_threshold_usd %||% 0), ]
  events |>
    dplyr::mutate(discovery = origin + lag, event_id = dplyr::row_number()) |>
    dplyr::arrange(type, origin, lag)
}
