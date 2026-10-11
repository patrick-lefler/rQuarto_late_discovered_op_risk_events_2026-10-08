# Cell (i, k): origin quarter i, discovery lag k (column index k + 1).
# Observed iff i + k <= n_origin. Lags >= n_lag go to the tail vector.
# Cell (i, k): origin period i, discovery lag k (column index k + 1).
# Observed iff i + k <= n_origin. Lags >= n_lag go to the tail vector.
# period > 1 coarsens to calendar blocks (e.g. 4 = annual): lag is the number of
# block boundaries crossed between origin and discovery.
add_periods <- function(events, period = 1) {
  events |>
    dplyr::mutate(
      origin_p = ceiling(origin / period),
      lag_p = ceiling((origin + lag) / period) - origin_p
    )
}

build_square <- function(events, type, value = c("count", "amount"), n_origin, n_lag, period = 1) {
  value <- match.arg(value)
  n_origin <- n_origin / period
  n_lag <- n_lag / period
  ev <- add_periods(events[events$type == type, ], period)
  ev$value <- if (value == "count") 1 else ev$severity
  sq <- matrix(0, n_origin, n_lag, dimnames = list(origin = seq_len(n_origin), lag = seq_len(n_lag) - 1))
  cells <- ev[ev$lag_p < n_lag, ] |> dplyr::summarise(value = sum(value), .by = c(origin_p, lag_p))
  sq[cbind(cells$origin_p, cells$lag_p + 1)] <- cells$value
  tail <- numeric(n_origin)
  tail_cells <- ev[ev$lag_p >= n_lag, ] |> dplyr::summarise(value = sum(value), .by = origin_p)
  tail[tail_cells$origin_p] <- tail_cells$value
  list(square = sq, tail = tail)
}

observed_mask <- function(n_origin, n_lag) {
  outer(seq_len(n_origin), seq_len(n_lag) - 1, \(i, k) i + k <= n_origin)
}

observed_triangle <- function(square) {
  tri <- square
  tri[!observed_mask(nrow(square), ncol(square))] <- NA
  ChainLadder::as.triangle(tri)
}

true_ibnr <- function(square, tail) {
  obs <- observed_mask(nrow(square), ncol(square))
  tibble::tibble(
    origin = seq_len(nrow(square)),
    reported = unname(rowSums(square * obs)),
    ibnr_window = unname(rowSums(square * !obs)),
    ibnr_tail = tail,
    ibnr_total = ibnr_window + ibnr_tail
  )
}
