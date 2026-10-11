# Installs the pinned mlReserve fork (and ranger, e1071) into lib_mlreserve/, separate from the renv library.
# Run from OUTSIDE the project's renv activation, e.g.:  cd /tmp && Rscript /path/to/install_mlreserve.R /path/to/project
project <- commandArgs(trailingOnly = TRUE)[1]
lib <- file.path(project, "lib_mlreserve")
dir.create(lib, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
commit <- yaml::read_yaml(file.path(project, "params.yaml"))$mlreserve$commit
options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages(c("ranger", "e1071", "proxy"), lib = lib)
remotes::install_github(paste0("thierrymoudiki/ChainLadder@", commit), lib = lib, upgrade = "never", force = TRUE)
