# Continuous Outcome Analysis

This folder contains the synthetic continuous outcome analysis for the
METRIcAL application. The analysis combines a latent implementation-strata
mixture model, a multilevel compliance model, and a continuous outcome model
fit by Gibbs sampling. The principal estimands are super-population
intention-to-treat effects and complier-average causal effects, marginally and
within the latent
implementation strata.

## File Structure

- `Continuous_Outcome_Analysis_main.R`
  Main script for running the multichain analysis.
- `Continuous_Outcome_Analysis_helpers.R`
  Workflow-level functions for single-chain and multichain sampling, post hoc
  relabeling of latent strata, posterior summaries, convergence diagnostics,
  trace plots, and posterior predictive checks.
- `Continuous_Outcome_Analysis_model_helpers.R`
  Conditional-update functions used within each Gibbs iteration.
- `CZ_Matrix_Synth_OG.csv`, `X_Synth_OG.csv`, and
  `Outcome_Synth_Continuous_OG.csv`
  Synthetic data used by the analysis.

The facility CSV is keyed by `Facility`. Its `C` (on the logit scale)
and `Z` (centered census scale) values were generated for illustration,
not matched to the older unkeyed file or to the private METRIcAL data.
A fixed seed (20261113) was used to assign two illustrative mixture groups
within each arm and generate correlated Gaussian facility variables with
group-specific means. The resident covariates, outcomes, and assignments
were not regenerated. These examples demonstrate the analysis workflow;
they are not intended to reproduce the manuscript results or to be samples
from the fitted joint model.

Treatment assignment comes from the resident outcome file. The scripts
check that assignment is constant within each facility, that `Z` is
complete, and that `C` is available for treated facilities only. The two
resident CSVs are paired by row order; do not sort them independently.

## Required Packages

The analysis uses:

- `boot`
- `mixtools`
- `coda`
- `MCMCpack`
- `DirichletReg`
- `mvtnorm`
- `MASS`
- `scales`
- `invgamma`
- `dplyr`
- `truncnorm`

## Running the Analysis

From R or RStudio, run:

```r
source("Continuous_Outcome_Analysis_main.R")
```

The main script runs four chains by default and saves the fitted object to:

```text
outputs/Continuous_Outcome_multichain_fit.rds
```

Trace plots for the primary estimands are saved to:

```text
outputs/Continuous_Outcome_multichain_traces.pdf
```

Tabular posterior summaries, diagnostics, and predictive checks are saved
as CSV files in `outputs/reviewer_tables/`.

## Main Analysis Function

The central function is:

```r
run_continuous_outcome_multichain(config, n_chains, chain_seeds, make_trace_pdf)
```

Inputs:

- `config`: controls the seed, burn-in, thinning, retained posterior draws,
  latent-strata dimension, and output directory.
- `n_chains`: number of independent MCMC chains.
- `chain_seeds`: chain-specific seeds for reproducible multichain sampling.
- `make_trace_pdf`: whether to save trace plots for the primary relabelled
  estimands.

Outputs:

- `chain_results`: retained posterior output from each chain.
- `pooled_estimand_draws`: combined relabelled posterior estimand draws across
  chains.
- `pooled_posterior_summary`: posterior means, standard deviations, and 95%
  credible intervals.
- `diagnostics`: Gelman--Rubin statistics and effective sample sizes.
- `latent_class_diagnostics`: posterior facility-membership probabilities,
  ambiguous-membership summaries, and posterior class-size summaries after
  relabeling.
- `pp_check_results`: posterior predictive summaries for the principal signal
  and noise contrasts used to assess model fit.

## Reproducibility Notes

The analysis is controlled by explicit seeds. The main script uses the seed in
the configuration object and then offsets it across chains.
