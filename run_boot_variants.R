# Bootstrap chain-ladder settings check: process distribution "gamma" (used in the main study) vs "od.pois".
# Baseline, annual dollar triangles, the first n_reps databases. Usage: Rscript run_boot_variants.R [n_reps]
suppressPackageStartupMessages({library(dplyr); library(ChainLadder)})
for (f in list.files("R", full.names = TRUE)) source(f)
args <- commandArgs(trailingOnly = TRUE)
params <- load_params()
n_reps <- if (length(args) >= 1) as.integer(args[1]) else 100L
base <- apply_scenario(params, "baseline")
n_cores <- max(1L, parallel::detectCores() - 2L)

one <- \(g) {
  events <- simulate_database(base, rep = g$rep)
  inp <- make_inputs(events, g$type, "amount", base, period = 4)
  out <- run_safely("boot", \(){
    b <- ChainLadder::BootChainLadder(inp$tri_cum, R = params$replication$n_boot, process.distr = g$process, seed = g$rep)
    sims <- as.numeric(b$IBNR.Totals)
    list(est = mean(sims), lo = stats::quantile(sims, 0.025, names = FALSE), hi = stats::quantile(sims, 0.975, names = FALSE), status = "ok", message = "")
  })
  truth <- sum(inp$truth$ibnr_window)
  tibble::tibble(rep = g$rep, type = g$type, process = g$process, status = out$status,
                 est = out$est %||% NA_real_, lo = out$lo %||% NA_real_, hi = out$hi %||% NA_real_, truth = truth)
}
grid <- tidyr::expand_grid(rep = seq_len(n_reps), type = c("trade_error", "nav_error"), process = c("gamma", "od.pois"))
out <- dplyr::bind_rows(parallel::mclapply(seq_len(nrow(grid)), \(i) one(grid[i, ]), mc.cores = n_cores))
saveRDS(out, "data/boot_variants.rds")
cat(sprintf("%d fits, %d failed\n", nrow(out), sum(out$status != "ok")))
