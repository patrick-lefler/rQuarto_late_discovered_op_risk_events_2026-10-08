load_params <- function(path = "params.yaml") {
  params <- yaml::read_yaml(path)
  validate_params(params)
  params
}

validate_params <- function(params) {
  stopifnot(
    params$horizon$n_origin > 0,
    params$horizon$n_lag > 0,
    length(params$event_types) > 0
  )
  for (nm in names(params$event_types)) {
    spec <- params$event_types[[nm]]
    stopifnot(
      spec$lambda > 0,
      spec$lag$w_fast >= 0, spec$lag$w_fast <= 1,
      spec$lag$p_fast > 0, spec$lag$p_fast <= 1,
      spec$severity$p_tail >= 0, spec$severity$p_tail <= 1,
      spec$severity$tail_q > 0, spec$severity$tail_q < 1,
      spec$severity$xi > 0, spec$severity$beta > 0
    )
  }
  invisible(params)
}

apply_scenario <- function(params, scenario = "baseline") {
  overrides <- params$scenarios[[scenario]]
  if (is.null(overrides)) stop("Unknown scenario: ", scenario)
  params$scenario_name <- scenario
  params$injected_events <- overrides$injected_events
  params$control_shift <- overrides$control_shift
  params$count_dist <- overrides$count_dist %||% "poisson"
  params$count_size <- overrides$count_size
  if (!is.null(overrides$reporting_threshold_usd)) params$reporting_threshold_usd <- overrides$reporting_threshold_usd
  if (!is.null(overrides$exposure)) {
    params$exposure <- utils::modifyList(params$exposure, overrides$exposure)
  }
  params
}

derive_seed <- function(master_seed, rep, scenario_index, type_index = 0L) {
  as.integer((master_seed + 1009 * rep + 100003 * scenario_index + 7919 * type_index) %% .Machine$integer.max)
}

# Additive shifts to the lag distribution of selected event types; w_fast is clamped to [0, 1].
apply_control <- function(params, control, scale = 1) {
  shifts <- params$controls[[control]]$lag_shift
  if (is.null(shifts)) stop("Unknown control: ", control)
  for (type in names(shifts)) {
    for (field in names(shifts[[type]])) {
      new <- params$event_types[[type]]$lag[[field]] + scale * shifts[[type]][[field]]
      if (field == "w_fast") new <- min(max(new, 0), 1)
      params$event_types[[type]]$lag[[field]] <- new
    }
  }
  validate_params(params)
  params$control_name <- control
  params
}
