# Late-Discovered Operational Risk Events

Estimating incurred-but-not-yet-recorded (IBNR-style) losses in a buy-side loss event database, using run-off triangle methods from the R `ChainLadder` package. All data are synthetic; see `project_spec.md` for the design and the site's *Limitations* page for what the results do and do not show.

## Run

Requirements: R >= 4.5, Quarto >= 1.8 (bundles Typst for the PDF brief), a C/Fortran toolchain for packages.

```sh
Rscript -e 'renv::restore()'   # once, restores package versions from renv.lock
Rscript run_all.R               # simulations -> mlReserve experiment -> render site to docs/
```

The one-page board brief (`brief.qmd`, Typst PDF) is deliberately not part of the website render; build it on demand with `quarto render brief.qmd`.

Individual steps: `Rscript run_simulations.R [n_reps] [n_cores]`, `Rscript run_mlreserve.R`, `quarto render`. Tests: `Rscript -e 'testthat::test_dir("tests/testthat")'`.

Runtime on a 12-core laptop (10 workers): about 8-10 minutes for the simulations, a few minutes for mlReserve and the render. Results are bit-reproducible for the same platform and package versions because every random stream is seeded explicitly (see the *Reproducibility* page).

## Layout

| Path | Purpose |
|---|---|
| `params.yaml` | Every simulator parameter, labelled as an assumption, with ranges |
| `R/` | Generator, triangle construction, estimators, scoring, controls, report helpers |
| `tests/testthat/` | Unit tests (triangle reconciliation, seeds, scoring, estimator behaviour) |
| `run_simulations.R`, `run_mlreserve.R`, `install_mlreserve.R`, `run_all.R` | Pipeline |
| `*.qmd`, `_quarto.yml` | Website; output goes to `docs/` |
| `data/` | Generated results (not committed; regenerate with the pipeline) |
| `lib_mlreserve/` | Isolated library holding the pinned `mlReserve` fork (not committed) |
