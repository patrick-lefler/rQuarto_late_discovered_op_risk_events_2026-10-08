# Project Spec: Late-Discovered Operational Risk Events — Estimating Incurred-But-Not-Yet-Recorded (IBNR) Losses for a Buy-Side Loss Event Database

**Language / stack:** R + Quarto website (ChainLadder, tidyverse, data.table, fitdistrplus, ggplot2, gt, yaml, testthat; pak for installs, renv for the environment). Native pipe and `.by` grouping throughout. Static ggplot2 charts by default; interactive plots (plotly) only with your explicit approval.
**Output:** Quarto website rendered to `docs/` (GitHub Pages ready)
**Audience:** Buy-side COOs, CROs, operational risk and ODD teams, prime-broker relationship and risk teams
**Positioning:** Demonstration (a test under stated assumptions, not a proof) of method for Seeking Eureka's buy-side operational risk advisory work
**Status:** Draft v1 for approval. No code is written until this spec is approved.

---

## 1. Problem Statement

A loss event database records events by the date they are *discovered* (or booked), not the date they *occurred*. Trade errors, NAV pricing errors, and unauthorized trading are often found weeks, quarters, or years after the origin date. Consequently:

- Recent origin periods always look cleaner than they are. A falling loss trend in the latest quarters may be an artifact of detection lag, not improved controls.
- Loss figures reported to boards, allocators, and in operational risk capital or loss-data submissions are biased low for recent periods.
- Managers cannot answer: *"How much loss has already been incurred but not yet recorded, and how uncertain is that number?"*

This is structurally the actuarial IBNR problem. An event originates at period t₀ and emerges across lags t₁ … tₖ. Run-off triangle mathematics applies directly, with discovery lag in the role of the claims reporting delay.

## 2. Proposed Solution

A reproducible R/Quarto study that:

1. Builds a **synthetic buy-side loss event database** with known ground truth, for three Basel-aligned event types.
2. Converts it into **count triangles** and **amount triangles** (origin quarter × discovery lag), one pair per event type.
3. Estimates undiscovered counts and losses with several methods from the ChainLadder package family.
4. **Scores each method against the known truth** over many simulated databases, which is only possible because the data are synthetic.
5. Stress-tests the methods where the chain-ladder assumptions are known to fail.
6. Quantifies the value of faster detection controls.

Core chain:

> Event origin → discovery lag distribution → loss recorded at discovery → triangle → IBNR estimate and predictive distribution → validation against truth → control value

## 3. Scope and Assumptions

**In scope**
- One stylized mid-size buy-side firm (about $5B AUM), quarterly origin periods over 10 years (40 origin quarters), quarterly discovery lags 0–39 (a 40×40 triangle; collapsible to annual for presentation).
- Three event types, each with its own triangles (Basel Level 1 mapping in brackets):
  | Event type | Basel category | Frequency | Typical lag | Severity |
  |---|---|---|---|---|
  | Trade errors / breaks | Execution, Delivery & Process Management | High | Short (days to weeks) | Small, mild tail |
  | NAV / pricing errors | Execution, Delivery & Process Management; Clients, Products & Business Practices | Medium | Medium (quarters) | Medium, moderate tail |
  | Unauthorized trading / fraud | Internal Fraud | Rare | Long (quarters to years) | Very large, heavy tail |
- Loss is **recognized in full at the discovery date** (see Section 4.2 for the extension that separates settlement development).
- Prime-broker reconciliation is modeled as a detection channel that shortens the lag distribution (see Section 4.3).

**Out of scope (documented extensions)**
- Settlement/payment development after discovery (a second triangle layer)
- Cyber and third-party outage events (see the separate ransomware spec)
- Fund-complex or multi-firm views
- Regulatory capital calculation (the Basel AMA is superseded by the SMA; this project does not compute capital)
- Real internal loss data

**Data policy:** Fully synthetic. Every simulator parameter is stored in a single `params.yaml` with a note labeling it as an **assumption** or as **anchored to a public source**, plus a plausible range. Public sources (for example ORX publications and public rogue-trading cases) may be used to set orders of magnitude only. No citation appears in the report unless it was verified; unverified values are labeled "assumption, to be replaced with client data."

## 4. Model Architecture

### 4.1 Synthetic data generator
For event type *e* and origin quarter *t*:

- **Frequency:** N(e,t) ~ Poisson(λₑ × exposure_t) (negative binomial as a robustness option). Exposure follows a growth path (trade volume or AUM proxy) so origin-period effects exist.
- **Discovery lag:** discrete mixture per event type (a fast-detection component plus a slow tail). The slow component is a discretized lognormal or Weibull. Lags are right-truncated at the evaluation date, which produces the triangle.
- **Severity:** lognormal body with a generalized Pareto tail, per event type. Severity can increase with lag (undetected errors compound), controlled by a single parameter, which makes the amount triangle differ from the count triangle.
- **Loss recorded at discovery:** event loss enters the cell for its discovery quarter.

The generator produces the **full square** (including the unobserved lower-right), so true IBNR is exactly known.

### 4.2 Triangles
- Incremental and cumulative **count** triangles and **amount** triangles by event type.
- Exposure-normalized variants for diagnosing origin-period effects.
- Tail beyond lag 39 is handled explicitly. Unauthorized trading can exceed the window, so a tail factor assumption is stated and tested.

### 4.3 Estimation methods
| Method | Package / function | Role |
|---|---|---|
| Mack chain-ladder | `ChainLadder::MackChainLadder` | Baseline with analytic standard error |
| Bootstrap chain-ladder | `ChainLadder::BootChainLadder` | Predictive distribution of IBNR |
| Clark LDF / Cape Cod | `ChainLadder::ClarkLDF`, `ClarkCapeCod` | Parametric development curves, better with sparse tails, uses exposure |
| Frequency × severity | counts from count triangle × fitted severity | Separates "how many are missing" from "how big" |
| Machine-learning reserving | `mlReserve()` from the `thierrymoudiki/ChainLadder` GitHub fork | Optional comparator (own page; see 4.5) |

Methods are applied to count triangles, amount triangles, or both, as appropriate. Methods that cannot run on sparse or zero-heavy triangles (for example the unauthorized-trading amount triangle) are reported as failing, with the reason, not silently dropped.

### 4.4 Scenarios
1. **Baseline:** stable detection and proportional development.
2. **Heavy-tail shock:** one or two very large unauthorized-trading losses inside the window.
3. **Control improvement (calendar effect):** detection lag shortens partway through the sample (for example after a new reconciliation control). Chain-ladder assumes stable development, so this tests the key assumption.
4. **Exposure growth:** volume growth that distorts origin-period comparisons.

### 4.5 Machine-learning comparator
`mlReserve` is an experimental fork that fits random forest or SVM learners to triangles. Published evidence for it is limited to simulated data. It goes on its own page, labeled experimental, and is excluded from headline results. If it fails to install or run, the page reports that and the rest of the site is unaffected.

### 4.6 Validation against truth
- Replicate each scenario over many simulated databases (target 200; adjustable in `params.yaml`).
- Metrics per method, event type, and origin quarter: bias, RMSE, relative error, and **coverage of the 95% interval**.
- Question answered: *which method can be trusted for which event type, and where does it break?*

### 4.7 Control value
- Shorten the lag distribution (for example, a faster reconciliation or pre-trade control) and measure the change in reported loss, in true loss, and in IBNR.
- Express the result as the reduction in expected undiscovered loss and in the 95th-percentile IBNR. Control costs are labeled assumptions.

## 5. Website Structure (Quarto website, Large tier: six pages, approved 2026-10-11)

| Page | File | Content |
|---|---|---|
| Overview | `index.qmd` | Executive summary, headline findings, the discovery-lag problem in pictures |
| Synthetic Data | `data.qmd` | Generator design, parameters, ground truth, scale sanity check, scenarios |
| Results | `results.qmd` | Summary by event type, count triangles, dollar triangles (incl. Bornhuetter-Ferguson), method validation, stress tests |
| Control Value | `controls.qmd` | Value of faster detection, break-even effectiveness |
| Reproducibility and Limitations | `reproducibility.qmd` | Section 6 below, plus limitations and model risk |
| ML appendix | `ml-appendix.qmd` | Experimental `mlReserve` comparator, moved off Results |

The one-page board brief (`brief.qmd`, PDF) is kept in the repository but excluded from the website render.

## 6. Reproducibility Page (dedicated)

The reproducibility section of `reproducibility.qmd` documents, and where possible *verifies live*:

- **One-command run** from a clean clone (for example `Rscript run_all.R`, or `make`).
- **Environment:** `renv.lock`, R and Quarto versions, `sessionInfo()` output rendered into the page, and the GitHub install pin (commit hash) for the `mlReserve` fork.
- **Randomness:** a single master seed in `params.yaml`, with a documented scheme for deriving per-replication seeds, so results are bit-reproducible on the same platform.
- **Parameters:** the full `params.yaml` rendered as a table with assumption/source labels and ranges.
- **Data lineage:** generated data are not committed as source of truth. The generator is. A hash (for example SHA-256) of each generated dataset is rendered, so a reader can confirm they reproduced the same data.
- **Tests:** `testthat` results rendered on the page, covering the triangle construction (a triangle built from the full square must reconcile exactly to the event-level data), the generator (marginal distributions, truncation), and the scoring functions.
- **Caching:** which computations are cached via Quarto `freeze`, how to invalidate them, and the runtime of a full from-scratch run.
- **Known limits of reproducibility:** cross-platform floating-point differences, and package version drift.

## 7. Deliverables

1. **Quarto website** in `docs/`, structured per Section 5.
2. **Reproducible repo:** `params.yaml`, `R/` functions, `tests/`, `renv.lock`, `README.md` with run instructions, `run_all.R`.
3. **One-page brief** (PDF, via Quarto's bundled Typst): three charts and three takeaways for a COO or CRO.

## 8. Key Charts (target set)
- Heatmap of observed vs. true triangle, with the unobserved region shaded.
- "Reported vs. true incurred loss by origin quarter," showing the apparent improvement in recent quarters as an artifact.
- Discovery-lag distributions by event type.
- IBNR predictive distributions by method (against the true value).
- Coverage and bias scorecard by method × event type × scenario.
- Control-shift stress: method estimates before and after the detection change.
- Value of faster detection: reduction in expected and 95th-percentile IBNR.

## 9. Work Plan (about 2–3 weeks part-time)

| Stage | Tasks |
|---|---|
| 1 | `params.yaml`; generator; ground-truth square; triangle construction; unit tests |
| 2 | Estimation wrappers (Mack, Boot, Clark, freq × sev); failure handling |
| 3 | Replication and scoring framework; scenarios; control-value analysis |
| 4 | `mlReserve` page; website pages; reproducibility page; one-page brief |
| 5 | Review: reconciliation checks, full clean-clone run, render to `docs/` |

## 10. Success Criteria
- The reproducibility page's checks pass on a clean clone, and the site renders from one command.
- Every parameter appears in `params.yaml` with an assumption/source label and a range.
- Triangles reconcile exactly to the event-level data (tested).
- The report states clearly, per event type, which methods are reliable (acceptable bias and near-nominal coverage) and which fail, with evidence.
- At least one quotable, defensible finding, stated with its conditions, for example: "In the baseline scenario, reported losses for the latest four origin quarters understate true incurred loss by about X%." The figure comes from the results; none is assumed in advance.

## 11. Commercial Hooks
- **Offer:** a "Loss Data Completeness Review," replacing synthetic parameters with the client's loss event database and measuring its discovery lags.
- **Targets:** buy-side COOs and CROs, ODD teams, and prime-broker risk teams that care about how fast breaks and errors are caught.
- **Linkage:** the same method extends to cyber and third-party outage events (see the separate ransomware spec).

## 12. Open Decisions
1. Whether to add a fourth event type, prime-broker settlement/reconciliation breaks, in v1 or defer it (recommend defer).
2. Replication count (200 proposed) against total runtime, and whether results are cached between renders.
3. Whether the brief's three takeaways are fixed in advance or drawn from the results (recommend drawn from the results).
4. GitHub Pages hosting details (repo name, whether `docs/` is committed), which you control.

## 13. Known Risks
- `mlReserve` is an unreleased GitHub fork. Install may fail, or the API may change, so it is pinned to a commit and isolated on one page.
- Chain-ladder methods can fail on triangles with zero columns or tiny counts, which is expected for unauthorized trading, and failures are reported as findings.
- Synthetic results show how methods behave *under the assumed generator*. They do not prove performance on any real firm's data. The limitations page says so.
