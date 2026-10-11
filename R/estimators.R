# Every estimator takes the same `inp` list (see make_inputs) and returns
#   list(method, status, message, by_origin, total)
# status: "ok" | "failed" (error or non-finite output) | "suspect" (ran but failed a sanity check).
# by_origin: origin, ibnr, se.  total: ibnr, se, lo, hi (95% interval), sims (optional).
# ibnr is the undiscovered amount *inside the development window*, except Clark with
# maxage = Inf, which also extrapolates a tail (reported as method "clark_ldf_tail").

make_inputs <- function(events_all, type, value, params, period = 1) {
  n_o <- params$horizon$n_origin
  n_l <- params$horizon$n_lag
  seen <- events_all[events_all$origin + events_all$lag <= n_o, ]
  full <- build_square(events_all, type, value, n_o, n_l, period)
  obs <- build_square(seen, type, value, n_o, n_l, period)
  exposure <- params$exposure$base * (1 + params$exposure$growth_per_quarter)^(seq_len(n_o) - 1)
  exposure_p <- as.numeric(tapply(exposure, ceiling(seq_len(n_o) / period), sum))
  list(
    type = type, value = value, period = period,
    n_origin = n_o / period, n_lag = n_l / period,
    tri_incr = observed_triangle(obs$square),
    tri_cum = ChainLadder::incr2cum(observed_triangle(obs$square)),
    exposure = exposure_p,
    events_obs = add_periods(seen[seen$type == type, ], period),
    truth = true_ibnr(full$square, full$tail)
  )
}

interval_lognormal <- function(mean, se, level = 0.95) {
  if (!is.finite(mean) || !is.finite(se) || mean <= 0) return(c(NA_real_, NA_real_))
  s2 <- log(1 + (se / mean)^2)
  z <- stats::qnorm(0.5 + level / 2)
  c(exp(log(mean) - s2 / 2 - z * sqrt(s2)), exp(log(mean) - s2 / 2 + z * sqrt(s2)))
}

run_safely <- function(method, expr_fn) {
  warns <- character()
  res <- withCallingHandlers(
    tryCatch(expr_fn(), error = \(e) e),
    warning = \(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = \(m) invokeRestart("muffleMessage")
  )
  if (inherits(res, "error")) {
    return(list(method = method, status = "failed", message = conditionMessage(res), by_origin = NULL, total = NULL))
  }
  res$method <- method
  res$n_warnings <- length(warns)
  res
}

check_output <- function(res) {
  ok <- is.finite(res$total$ibnr) && all(is.finite(res$by_origin$ibnr))
  if (!ok) {
    res$status <- "failed"
    res$message <- "non-finite output"
  }
  res
}

new_result <- function(by_origin, total_ibnr, total_se, sims = NULL, interval = NULL) {
  if (is.null(interval)) {
    interval <- if (is.null(sims)) interval_lognormal(total_ibnr, total_se) else stats::quantile(sims, c(0.025, 0.975), names = FALSE)
  }
  list(
    status = "ok", message = "",
    by_origin = by_origin,
    total = list(ibnr = total_ibnr, se = total_se, lo = interval[1], hi = interval[2], sims = sims)
  )
}

# Volume-weighted development factors; a factor with an all-zero denominator is set to 1.
# Unlike ChainLadder::chainladder() this cannot fail on sparse triangles; it returns point estimates only.
project_chainladder <- function(tri_cum) {
  m <- unclass(tri_cum)
  n_o <- nrow(m)
  n_l <- ncol(m)
  f <- vapply(seq_len(n_l - 1), \(k) {
    rows <- !is.na(m[, k + 1])
    den <- sum(m[rows, k])
    if (den > 0) sum(m[rows, k + 1]) / den else 1
  }, numeric(1))
  full <- m
  for (i in seq_len(n_o)) for (k in seq_len(n_l - 1)) {
    if (is.na(full[i, k + 1])) full[i, k + 1] <- full[i, k] * f[k]
  }
  full
}

est_chainladder <- function(inp) {
  run_safely("chainladder_vw", \(){
    full <- project_chainladder(inp$tri_cum)
    latest <- vapply(seq_len(nrow(full)), \(i) full[i, ncol(full) - i + 1], numeric(1))
    ibnr <- unname(full[, ncol(full)] - latest)
    check_output(new_result(tibble::tibble(origin = seq_along(ibnr), ibnr = ibnr, se = NA_real_), sum(ibnr), NA_real_))
  })
}

est_mack <- function(inp) {
  run_safely("mack", \(){
    m <- ChainLadder::MackChainLadder(inp$tri_cum, est.sigma = "Mack")
    s <- summary(m)
    by <- tibble::tibble(origin = seq_len(nrow(inp$tri_cum)), ibnr = unname(s$ByOrigin$IBNR), se = unname(s$ByOrigin$Mack.S.E))
    tot <- s$Totals
    check_output(new_result(by, tot["IBNR:", 1], tot["Mack S.E.:", 1]))
  })
}

est_boot <- function(inp, R = 1000, seed = 1L) {
  run_safely("boot", \(){
    b <- ChainLadder::BootChainLadder(inp$tri_cum, R = R, process.distr = "gamma", seed = seed)
    s <- summary(b)$ByOrigin
    idx <- match(rownames(inp$tri_cum), rownames(s))
    by <- tibble::tibble(origin = seq_len(nrow(inp$tri_cum)), ibnr = s[["Mean IBNR"]][idx], se = s[["SD IBNR"]][idx])
    by[is.na(idx), c("ibnr", "se")] <- 0
    sims <- as.numeric(b$IBNR.Totals)
    res <- check_output(new_result(by, mean(sims), stats::sd(sims), sims = sims))
    ref <- est_chainladder(inp)
    if (res$status == "ok" && ref$status == "ok" && ref$total$ibnr > 0 && res$total$ibnr < 0.1 * ref$total$ibnr) {
      res$status <- "suspect"
      res$message <- sprintf("bootstrap mean IBNR %.3g is <10%% of deterministic chain-ladder %.3g", res$total$ibnr, ref$total$ibnr)
    }
    res
  })
}

clark_core <- function(fit) {
  by <- tibble::tibble(origin = seq_along(fit$FutureValue), ibnr = unname(fit$FutureValue), se = unname(fit$StdError))
  res <- check_output(new_result(by, sum(by$ibnr), as.numeric(fit$Total["StdError"])))
  par <- fit$PAR[c("omega", "theta")]
  at_bound <- any(is.finite(par) & (par <= 0.011 | abs(par - 0.5) < 1e-6 & names(par) == "theta"))
  if (res$status == "ok" && isTRUE(at_bound)) {
    res$status <- "suspect"
    res$message <- sprintf("growth-curve parameter at optimiser bound (omega = %.3g, theta = %.3g)", par[["omega"]], par[["theta"]])
  }
  res
}

est_clark <- function(inp, G = "loglogistic", tail = FALSE, cape_cod = FALSE) {
  name <- paste0(if (cape_cod) "clark_cc_" else "clark_ldf_", G, if (tail) "_tail" else "")
  run_safely(name, \(){
    maxage <- if (tail) Inf else inp$n_lag
    tri <- inp$tri_cum
    # Clark ages are period ends (1..n), not our zero-based lags
    dimnames(tri) <- list(origin = seq_len(nrow(tri)), dev = seq_len(ncol(tri)))
    # Cape Cod needs exposure on the loss scale so the expected loss ratio starts near 1
    latest <- vapply(seq_len(nrow(tri)), \(i) tri[i, ncol(tri) - i + 1], numeric(1))
    premium <- inp$exposure * sum(latest) / sum(inp$exposure)
    fit <- NULL
    utils::capture.output(fit <- if (cape_cod) {
      ChainLadder::ClarkCapeCod(tri, Premium = premium, G = G, maxage = maxage)
    } else {
      ChainLadder::ClarkLDF(tri, G = G, maxage = maxage)
    })
    clark_core(fit)
  })
}

observed_mean_severity <- function(events_obs, n_lag, bandwidth = 1, min_n = 5) {
  pooled <- mean(events_obs$severity)
  purrr::map(0:(n_lag - 1), \(k) {
    near <- events_obs$severity[abs(events_obs$lag_p - k) <= bandwidth]
    if (length(near) >= min_n) near else events_obs$severity
  })
}

# Frequency x severity: project undiscovered counts per future cell with chain-ladder on the
# *count* triangle, then price each projected event by resampling observed severities near that lag.
# Parameter uncertainty is ignored by design (so expect under-coverage).
est_freqsev <- function(inp_count, inp_amount, n_sim = 500, seed = 1L) {
  run_safely("freq_x_sev", \(){
    full_cum <- project_chainladder(inp_count$tri_cum)
    full_incr <- cbind(full_cum[, 1], full_cum[, -1] - full_cum[, -ncol(full_cum)])
    obs <- observed_mask(inp_count$n_origin, inp_count$n_lag)
    lambda <- full_incr * !obs
    lambda[!is.finite(lambda) | lambda < 0] <- 0
    pools <- observed_mean_severity(inp_amount$events_obs, inp_count$n_lag)
    set.seed(seed)
    sims <- matrix(0, n_sim, inp_count$n_origin)
    for (i in seq_len(inp_count$n_origin)) for (k in seq_len(inp_count$n_lag)) {
      if (lambda[i, k] > 0) {
        n <- stats::rpois(n_sim, lambda[i, k])
        sims[, i] <- sims[, i] + vapply(n, \(m) if (m > 0) sum(sample(pools[[k]], m, replace = TRUE)) else 0, numeric(1))
      }
    }
    mean_sev <- vapply(pools, mean, numeric(1))
    by <- tibble::tibble(origin = seq_len(inp_count$n_origin), ibnr = as.numeric(lambda %*% mean_sev), se = apply(sims, 2, stats::sd))
    tot_sims <- rowSums(sims)
    check_output(new_result(by, sum(by$ibnr), stats::sd(tot_sims), sims = tot_sims))
  })
}

all_estimators <- function(inp_count, inp_amount, boot_R = 1000, seed = 1L) {
  list(
    amount = list(
      chainladder = est_chainladder(inp_amount),
      mack = est_mack(inp_amount),
      boot = est_boot(inp_amount, R = boot_R, seed = seed),
      clark_ldf_loglogistic = est_clark(inp_amount, "loglogistic"),
      clark_ldf_weibull = est_clark(inp_amount, "weibull"),
      clark_ldf_loglogistic_tail = est_clark(inp_amount, "loglogistic", tail = TRUE),
      clark_cc_weibull = est_clark(inp_amount, "weibull", cape_cod = TRUE),
      freq_x_sev = est_freqsev(inp_count, inp_amount, seed = seed)
    ),
    count = list(
      chainladder = est_chainladder(inp_count),
      mack = est_mack(inp_count),
      boot = est_boot(inp_count, R = boot_R, seed = seed),
      clark_ldf_loglogistic = est_clark(inp_count, "loglogistic"),
      clark_ldf_weibull = est_clark(inp_count, "weibull")
    )
  )
}

summarise_estimators <- function(res_list) {
  purrr::imap(res_list, \(r, nm) tibble::tibble(
    method = nm, status = r$status, message = r$message %||% "",
    ibnr = r$total$ibnr %||% NA_real_, se = r$total$se %||% NA_real_,
    lo = r$total$lo %||% NA_real_, hi = r$total$hi %||% NA_real_
  )) |> purrr::list_rbind()
}
