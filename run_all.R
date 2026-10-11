# One command from a clean clone:  Rscript run_all.R
# Prerequisite: restore the R environment once with  Rscript -e 'renv::restore()'
stopifnot(file.exists("params.yaml"))
run <- function(...) {
  status <- system2("Rscript", c(...))
  if (status != 0) stop("Step failed: Rscript ", paste(c(...), collapse = " "))
}
# 1. main simulations (CRAN ChainLadder)
run("run_simulations.R")
# 2. experimental mlReserve fork, installed in an isolated library outside renv
if (!file.exists("lib_mlreserve/ChainLadder/DESCRIPTION")) {
  run("--no-init-file", "install_mlreserve.R", shQuote(getwd()))
}
run("run_mlreserve.R")
# 3. render the website to docs/
if (system2("quarto", "render") != 0) stop("quarto render failed")
