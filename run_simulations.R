# Runs the replication study and the control study, writing results to data/.
# Usage: Rscript run_simulations.R [n_reps] [n_cores]
suppressPackageStartupMessages({library(dplyr); library(ChainLadder)})
for (f in list.files("R", full.names = TRUE)) source(f)
args <- commandArgs(trailingOnly = TRUE)
params <- load_params()
n_reps <- if (length(args) >= 1) as.integer(args[1]) else params$replication$n_reps
n_cores <- if (length(args) >= 2) as.integer(args[2]) else max(1L, parallel::detectCores() - 2L)
dir.create("data", showWarnings = FALSE)

t_rep <- system.time(reps <- run_replications(params, n_reps, n_cores))
t_ctl <- system.time({
  ctl <- run_control_study(params, n_reps, n_cores)
  be <- run_breakeven(params, n_reps = n_reps, n_cores = n_cores)
})

saveRDS(reps, "data/replications.rds")
saveRDS(ctl, "data/control_study.rds")
saveRDS(be, "data/control_breakeven.rds")
saveRDS(list(
  n_reps = n_reps, n_cores = n_cores, replication_seconds = t_rep[["elapsed"]], control_seconds = t_ctl[["elapsed"]],
  params_md5 = unname(tools::md5sum("params.yaml")), finished = format(Sys.time(), tz = "UTC", usetz = TRUE)
), "data/run_meta.rds")
cat(sprintf("replications: %d rows, %.0fs | controls: %.0fs\n", nrow(reps), t_rep[["elapsed"]], t_ctl[["elapsed"]]))
cat("task errors:", sum(reps$status == "task_error", na.rm = TRUE), "\n")
