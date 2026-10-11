# Experimental comparator. Uses the pinned mlReserve fork from an isolated library (lib_mlreserve/),
# so the CRAN ChainLadder used everywhere else is untouched.
# Usage: Rscript run_mlreserve.R   (after: Rscript install_mlreserve.R)
lib <- file.path(getwd(), "lib_mlreserve")
if (!dir.exists(file.path(lib, "ChainLadder"))) stop("Run install_mlreserve.R first")
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages({library(dplyr); library(ChainLadder)})
for (f in list.files("R", full.names = TRUE)) source(f)
params <- load_params()
cfg <- params$mlreserve
learner <- list(ranger = ranger::ranger, svm = e1071::svm)
base <- apply_scenario(params, "baseline")
n_cores <- max(1L, parallel::detectCores() - 2L)

one <- \(g) {
  events <- simulate_database(base, rep = g$rep)
  inp <- make_inputs(events, g$type, "amount", base, period = 4)
  c <- cfg$configs[[g$config]]
  t <- system.time(res <- tryCatch({
    m <- mlReserve(inp$tri_cum, fit_func = learner[[c$learner]], features = c$features, transform = c$transform,
                   mse.method = "none", seed = g$rep)
    list(est = sum(m$summary$IBNR[rownames(m$summary) != "total"], na.rm = TRUE), status = "ok", message = "")
  }, error = \(e) list(est = NA_real_, status = "failed", message = conditionMessage(e))))
  tibble::tibble(rep = g$rep, type = g$type, config = g$config, est = res$est, status = res$status,
                 message = substr(res$message, 1, 80), truth = sum(inp$truth$ibnr_window), seconds = t[["elapsed"]])
}
grid <- tidyr::expand_grid(rep = seq_len(cfg$n_reps), type = names(params$event_types), config = names(cfg$configs))
out <- dplyr::bind_rows(parallel::mclapply(seq_len(nrow(grid)), \(i) one(grid[i, ]), mc.cores = n_cores))
saveRDS(list(results = out, commit = cfg$commit, fork_version = as.character(packageVersion("ChainLadder")),
             ranger = as.character(packageVersion("ranger")), e1071 = as.character(packageVersion("e1071"))),
        "data/mlreserve.rds")
cat(sprintf("%d fits, %d failed\n", nrow(out), sum(out$status != "ok")))
