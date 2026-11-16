# Cluster Compliance: Supplementary Code

This repository accompanies the manuscript on Bayesian analysis of pragmatic
cluster randomized trials with individual noncompliance and latent
implementation strata. It contains the six simulation case studies and two
synthetic examples of the METRIcAL analysis.

## Requirements

Install the R packages used by the scripts:

```r
install.packages(c(
  "boot", "mixtools", "coda", "dplyr", "MCMCpack", "DirichletReg",
  "mvtnorm", "MASS", "scales", "invgamma", "truncnorm", "sn"
))
```

`rstudioapi` is optional for RStudio use. The scripts and their inputs do not
require network access after package installation.

The scripts were checked with R 4.4.1 and these direct dependency versions:
`boot` 1.3-30, `mixtools` 2.0.0.1, `coda` 0.19-4.1, `dplyr` 1.1.4,
`MCMCpack` 1.7-1, `DirichletReg` 0.7-2, `mvtnorm` 1.3-3,
`MASS` 7.3-60.2, `scales` 1.4.0, `invgamma` 1.2,
`truncnorm` 1.0-9, and `sn` 2.1.1. Exact seeded draws may depend on the
R and package versions and numerical libraries.

## Run the supplement

From the repository root in a terminal:

```sh
Rscript --vanilla Simulated_Case_Studies/Simulated_Case_Studies_main.R case1
Rscript --vanilla Simulated_Case_Studies/Simulated_Case_Studies_main.R case6
Rscript --vanilla Continuous_Outcome_Analysis/Continuous_Outcome_Analysis_main.R
Rscript --vanilla Binary_Outcome_Analysis/Binary_Outcome_Analysis_main.R
```

Replace `case1` with `case2` through `case6` to select another simulation.
Each script locates its inputs and helpers relative to its own file. Results
are written under `outputs/` within the corresponding analysis or simulation
folder. Full simulation runs are computationally intensive.

## Contents

- [Simulated Case Studies](Simulated_Case_Studies/README.md): six simulation
  settings, default run sizes, operating characteristics, and a standalone
  Monte Carlo comparison with the appendix formulas.
- [Continuous Outcome Analysis](Continuous_Outcome_Analysis/README.md):
  multichain analysis of synthetic change in CMAI.
- [Binary Outcome Analysis](Binary_Outcome_Analysis/README.md):
  multichain probit analysis of synthetic antipsychotic use.

The analysis inputs are synthetic and are not intended to reproduce the
manuscript's METRIcAL results. The simulation uses the age and ADL covariate
input included in the original supplement. No private METRIcAL analysis data,
saved application fits, or reviewer materials are included.
