# Binary Outcome Analysis

## Overview

This folder contains the synthetic binary outcome analysis for the METRIcAL
application. The outcome is antipsychotic use at follow-up, modeled with a
probit regression. Individual compliance is also modeled with a multilevel
probit regression, and clusters are classified into latent implementation
strata through a finite Normal mixture model for the facility-level
implementation measure and baseline facility characteristic.

The analysis targets super-population ITT and CACE estimands, marginally and
within latent implementation strata. At each retained posterior draw it
integrates the Normal random effects in the probit models and averages the
resulting outcome risks over the empirical covariates within each stratum.
Marginal risks use the mixture weights, adjusted by stratum compliance
probabilities among compliers. Posterior predictive checks still use
simulated replicate datasets.

## File Structure

- `Binary_Outcome_Analysis_main.R`: main script for running the multichain analysis.
- `Binary_Outcome_Analysis_helpers.R`: single-chain and multichain workflow, relabeling, posterior summaries, diagnostics, trace plots, and posterior predictive checks.
- `Binary_Outcome_Analysis_model_helpers.R`: Gibbs update functions, latent probit augmentation, imputation steps, analytic estimand calculations, and PPC simulation.
- `CZ_Matrix_Synth_OG.csv`, `X_Synth_OG.csv`, and `Outcome_Synth_Binary_OG.csv`: synthetic data used by the analysis.
- `README.md`: this file.

## Required Data

- `CZ_Matrix_Synth_OG.csv`
- `X_Synth_OG.csv`
- `Outcome_Synth_Binary_OG.csv`

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

## Required R Packages

The workflow uses:

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
source("Binary_Outcome_Analysis_main.R")
```

The main script runs four chains by default and saves the full fitted object to:

```text
outputs/Binary_Outcome_multichain_fit.rds
```

Trace plots are saved to:

```text
outputs/Binary_Outcome_multichain_traces.pdf
```

Tabular posterior summaries, diagnostics, and predictive checks are saved
as CSV files in `outputs/reviewer_tables/`.

## Main Analysis Function

The central function is:

```r
run_binary_outcome_multichain(config, n_chains, chain_seeds, make_trace_pdf)
```

Inputs:

- `config`: controls the seed, burn-in, thinning, retained posterior draws, and output directory.
- `n_chains`: number of independent MCMC chains.
- `chain_seeds`: chain-specific seeds for reproducible multichain sampling.
- `make_trace_pdf`: whether to save trace plots for the primary relabelled estimands.

Outputs:

- `chain_results`: retained posterior output from each chain.
- `pooled_estimand_draws`: combined relabelled posterior estimand draws across chains.
- `pooled_posterior_summary`: posterior means, standard deviations, and 95% credible intervals.
- `diagnostics`: Gelman--Rubin statistics and effective sample sizes.
- `latent_class_diagnostics`: posterior facility-membership probabilities,
  ambiguous-membership summaries, and posterior class-size summaries after
  relabeling.
- `pp_check_results`: posterior predictive p-values for signal, noise, and signal-to-noise discrepancy measures.

## Label Switching

Latent mixture labels are not identifiable without an ordering convention. The sampler is left unconstrained, and labels are corrected after sampling. For each retained posterior draw, strata are ordered by the first component of the mixture mean vector, corresponding to the facility-level implementation measure. This makes stratum-specific estimands and PPC summaries comparable across retained draws and chains. Marginal estimands such as ITT and CACE are unchanged by relabeling.

## Posterior Predictive Checks

Posterior predictive checks are computed on the combined retained draws across chains. The checks compare replicated and observed signal, noise, and signal-to-noise statistics among compliers within each latent implementation stratum.

## Reproducibility Notes

The analysis is controlled by explicit seeds. The main script uses the seed in
the configuration object and then offsets it across chains.
