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
    truth = true_ibnr(full$square, full$tail),
    bench_unreported = benchmark_unreported(params, type, value, period)
  )
}

# Share of each origin period's in-window loss (or count) that is still undiscovered at the evaluation date,
# implied by the lag distribution in params. Used as an externally supplied (benchmark) development pattern for
# Bornhuetter-Ferguson; in this study it is the generator's own pattern, so it is an optimistic benchmark.
benchmark_unreported <- function(params, type, value, period = 1) {
  n_o <- params$horizon$n_origin
  n_l <- params$horizon$n_lag
  spec <- params$event_types[[type]]
  lag <- 0:(n_o + n_l)
  slow <- stats::plnorm(lag + 1, spec$lag$slow_meanlog, spec$lag$slow_sdlog) - stats::plnorm(lag, spec$lag$slow_meanlog, spec$lag$slow_sdlog)
  pmf <- spec$lag$w_fast * stats::dgeom(lag, spec$lag$p_fast) + (1 - spec$lag$w_fast) * slow
  w_lag <- pmf * if (value == "amount") 1 + spec$severity$lag_elasticity * lag else 1
  n_p <- n_o / period
  n_lp <- n_l / period
  exposure <- params$exposure$base * (1 + params$exposure$growth_per_quarter)^(seq_len(n_o) - 1)
  vapply(seq_len(n_p), \(i) {
    q <- ((i - 1) * period + 1):(i * period)
    d <- outer(q, lag, `+`)
    in_win <- ceiling(d / period) <= i + n_lp - 1
    w <- outer(exposure[q], w_lag)
    sum(w[in_win & d > n_o]) / sum(w[in_win])
  }, numeric(1))
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
    if (res$status == "ok" && res$total$ibnr <= 0) {
      res$status <- "suspect"
      res$message <- "bootstrap mean IBNR is zero (no development information)"
    } else if (res$status == "ok" && ref$status == "ok" && ref$total$ibnr > 0 && res$total$ibnr < 0.1 * ref$total$ibnr) {
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
  res$detail <- list(omega = unname(par[["omega"]]), theta = unname(par[["theta"]]), elr = unname(fit$PAR["ELR"]))
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

# Mack (1993) chain-ladder with an explicit factor of 1 for development columns with no volume, so it runs on sparse
# triangles where ChainLadder::MackChainLadder() stops on a non-finite factor. `n_diag` restricts the development
# factors (and variance parameters) to the latest n calendar diagonals. Verified against the package on RAA.
mack_manual <- function(tri_cum, n_diag = Inf) {
  m <- unclass(tri_cum)
  n <- nrow(m)
  stopifnot(ncol(m) == n)
  diag_idx <- outer(seq_len(n), seq_len(n), `+`)
  keep <- diag_idx > n + 1 - n_diag
  f <- numeric(n - 1)
  sig2 <- numeric(n - 1)
  for (k in seq_len(n - 1)) {
    rows <- which(!is.na(m[, k + 1]) & keep[, k + 1])
    den <- sum(m[rows, k])
    f[k] <- if (den > 0) sum(m[rows, k + 1]) / den else 1
    pos <- rows[m[rows, k] > 0]
    if (length(rows) > 1 && length(pos) > 0) {
      sig2[k] <- sum(m[pos, k] * (m[pos, k + 1] / m[pos, k] - f[k])^2) / (length(rows) - 1)
    }
  }
  if (n >= 4) {
    a <- sig2[n - 2]
    b <- sig2[n - 3]
    sig2[n - 1] <- if (b > 0) min(a^2 / b, a, b) else 0
  }
  full <- m
  for (i in seq_len(n)) for (k in seq_len(n - 1)) if (is.na(full[i, k + 1])) full[i, k + 1] <- full[i, k] * f[k]
  latest <- vapply(seq_len(n), \(i) m[i, n + 1 - i], numeric(1))
  ult <- full[, n]
  ibnr <- ult - latest
  S <- vapply(seq_len(n - 1), \(k) sum(m[seq_len(n - k), k]), numeric(1))
  inv0 <- \(x) ifelse(x > 0, 1 / x, 0)
  term <- sig2 / f^2
  mse <- numeric(n)
  for (i in 2:n) {
    ks <- (n + 1 - i):(n - 1)
    mse[i] <- ult[i]^2 * sum(term[ks] * (inv0(full[i, ks]) + inv0(S[ks])))
  }
  cross <- 0
  for (i in 2:(n - 1)) {
    ks <- (n + 1 - i):(n - 1)
    cross <- cross + ult[i] * sum(ult[(i + 1):n]) * sum(term[ks] * inv0(S[ks]))
  }
  # standardised age-to-age residuals on every observed transition, indexed by calendar diagonal i + k
  resid <- matrix(NA_real_, n, n - 1)
  for (k in seq_len(n - 1)) for (i in seq_len(n)) {
    if (!is.na(m[i, k + 1]) && m[i, k] > 0 && sig2[k] > 0) resid[i, k] <- (m[i, k + 1] - f[k] * m[i, k]) / sqrt(sig2[k] * m[i, k])
  }
  list(ibnr = unname(ibnr), se = unname(sqrt(mse)), total_ibnr = unname(sum(ibnr)), total_se = unname(sqrt(sum(mse) + 2 * cross)),
       factors = f, resid = resid)
}

# Calendar-diagonal shift diagnostic: mean standardised residual on the latest `n_last` diagonals, as a z-score.
# Strongly negative = recent development is below what the earlier factors predict (e.g. faster detection).
calendar_z <- function(tri_cum, n_last = 2) {
  r <- mack_manual(tri_cum)$resid
  n <- nrow(r)
  diag_idx <- outer(seq_len(n), seq_len(n - 1), `+`)
  z <- r[diag_idx > n - n_last & diag_idx <= n]
  z <- z[is.finite(z)]
  if (length(z) < 3) return(NA_real_)
  mean(z) * sqrt(length(z))
}

est_mack_manual <- function(inp, name = "mack_ec", n_diag = Inf) {
  run_safely(name, \(){
    mk <- mack_manual(inp$tri_cum, n_diag)
    by <- tibble::tibble(origin = seq_along(mk$ibnr), ibnr = mk$ibnr, se = mk$se)
    out <- check_output(new_result(by, mk$total_ibnr, mk$total_se))
    if (name == "mack_ec") out$detail <- list(cal_z = calendar_z(inp$tri_cum))
    out
  })
}

# Bornhuetter-Ferguson and Benktander (k = 1) on the in-window amount triangle.
# A-priori ultimate = expected loss per unit exposure (from the 5 oldest origin periods, whose development is nearly
# complete) x exposure. Percent unreported comes either from the triangle's own chain-ladder pattern ("own") or from
# an externally supplied pattern ("benchmark", inp$bench_unreported). Point estimates only: no interval is produced.
est_bf <- function(inp, pattern = c("own", "benchmark"), method = c("bf", "benktander")) {
  pattern <- match.arg(pattern)
  method <- match.arg(method)
  name <- paste0(method, if (pattern == "benchmark") "_benchmark" else "")
  run_safely(name, \(){
    full <- project_chainladder(inp$tri_cum)
    n <- nrow(full)
    latest <- vapply(seq_len(n), \(i) full[i, n + 1 - i], numeric(1))
    f <- vapply(seq_len(n - 1), \(k) {
      m <- unclass(inp$tri_cum)
      rows <- !is.na(m[, k + 1])
      den <- sum(m[rows, k])
      if (den > 0) sum(m[rows, k + 1]) / den else 1
    }, numeric(1))
    cdf <- vapply(seq_len(n), \(i) if (i == 1) 1 else prod(f[(n + 1 - i):(n - 1)]), numeric(1))
    unrep <- if (pattern == "own") 1 - 1 / cdf else inp$bench_unreported
    oldest <- seq_len(min(5, n))
    elr <- sum(full[oldest, n]) / sum(inp$exposure[oldest])
    prior <- elr * inp$exposure
    bf_ibnr <- prior * unrep
    ibnr <- if (method == "bf") bf_ibnr else unrep * (latest + bf_ibnr)
    check_output(new_result(tibble::tibble(origin = seq_len(n), ibnr = ibnr, se = NA_real_), sum(ibnr), NA_real_))
  })
}

# Severity donors for each lag: observed events within `bandwidth` periods of the lag; if fewer than `min_n`, all
# observed events are used (which drops the lag uplift). `pooled = "uplift"` instead de-uplifts every observed severity
# by an estimated linear lag effect e (OLS slope/intercept of severity on lag, floored at 0) and re-applies it at the
# projected lag, so all events inform every lag.
severity_pools <- function(events_obs, n_lag, bandwidth = 1, min_n = 5, pooled = c("nearby", "uplift")) {
  pooled <- match.arg(pooled)
  if (pooled == "uplift") {
    fit <- stats::lm(severity ~ lag_p, data = events_obs)
    e <- max(unname(stats::coef(fit)[2] / stats::coef(fit)[1]), 0)
    if (!is.finite(e)) e <- 0
    base <- events_obs$severity / (1 + e * events_obs$lag_p)
    return(list(pools = purrr::map(0:(n_lag - 1), \(k) base * (1 + e * (k + 0.5))), n_donors = rep(length(base), n_lag), fallback = rep(FALSE, n_lag), e = e))
  }
  idx <- purrr::map(0:(n_lag - 1), \(k) which(abs(events_obs$lag_p - k) <= bandwidth))
  fb <- lengths(idx) < min_n
  list(
    pools = purrr::map2(idx, fb, \(i, f) if (f) events_obs$severity else events_obs$severity[i]),
    n_donors = lengths(idx), fallback = fb, e = NA_real_
  )
}

# Frequency x severity: project undiscovered counts per future cell with chain-ladder on the
# *count* triangle, then price each projected event by resampling observed severities (see severity_pools).
# Parameter uncertainty is ignored by design (so expect under-coverage).
est_freqsev <- function(inp_count, inp_amount, n_sim = 500, seed = 1L, pooled = c("nearby", "uplift")) {
  pooled <- match.arg(pooled)
  run_safely(if (pooled == "nearby") "freq_x_sev" else "freq_x_sev_uplift", \(){
    full_cum <- project_chainladder(inp_count$tri_cum)
    full_incr <- cbind(full_cum[, 1], full_cum[, -1] - full_cum[, -ncol(full_cum)])
    obs <- observed_mask(inp_count$n_origin, inp_count$n_lag)
    lambda <- full_incr * !obs
    lambda[!is.finite(lambda) | lambda < 0] <- 0
    sp <- severity_pools(inp_amount$events_obs, inp_count$n_lag, pooled = pooled)
    pools <- sp$pools
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
    out <- check_output(new_result(by, sum(by$ibnr), stats::sd(tot_sims), sims = tot_sims))
    # share of projected events priced from the all-event fallback pool (no nearby donors), and the fewest donors used
    w <- colSums(lambda)
    out$detail <- list(fallback_share = if (sum(w) > 0) sum(w[sp$fallback]) / sum(w) else NA_real_,
                       min_donors = if (any(w > 0)) min(sp$n_donors[w > 0]) else NA_integer_, uplift = sp$e)
    out
  })
}

all_estimators <- function(inp_count, inp_amount, boot_R = 1000, seed = 1L) {
  list(
    amount = list(
      chainladder = est_chainladder(inp_amount),
      mack = est_mack(inp_amount),
      mack_ec = est_mack_manual(inp_amount, "mack_ec"),
      mack_last3 = est_mack_manual(inp_amount, "mack_last3", n_diag = 3),
      bf = est_bf(inp_amount, "own", "bf"),
      benktander = est_bf(inp_amount, "own", "benktander"),
      bf_benchmark = est_bf(inp_amount, "benchmark", "bf"),
      benktander_benchmark = est_bf(inp_amount, "benchmark", "benktander"),
      boot = est_boot(inp_amount, R = boot_R, seed = seed),
      clark_ldf_loglogistic = est_clark(inp_amount, "loglogistic"),
      clark_ldf_weibull = est_clark(inp_amount, "weibull"),
      clark_ldf_loglogistic_tail = est_clark(inp_amount, "loglogistic", tail = TRUE),
      clark_cc_weibull = est_clark(inp_amount, "weibull", cape_cod = TRUE),
      freq_x_sev = est_freqsev(inp_count, inp_amount, seed = seed),
      freq_x_sev_uplift = est_freqsev(inp_count, inp_amount, seed = seed, pooled = "uplift")
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
