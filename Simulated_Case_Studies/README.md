# Simulated Case Studies

The simulations evaluate posterior inference for super-population
intent-to-treat (ITT) and complier average causal effects (CACE), marginally
and within latent implementation strata.

| Argument | Data-generating mechanism |
| --- | --- |
| `case1` | Correctly specified model |
| `case2` | Skew-$t$ mixture for facility-level variables |
| `case3` | Burr link for individual compliance |
| `case4` | Covariate interactions in the outcome model |
| `case5` | Weaker separation between latent implementation strata |
| `case6` | Case 1 mechanism with cluster sizes of 10 and 30, balanced across arms |

## Files

- `Simulated_Case_Studies_main.R`: evaluates the case-specific truth, fits
  the joint Bayesian model by Gibbs sampling, and summarizes operating
  characteristics across simulated datasets.
- `Simulated_Case_Studies_utils.R`: defines the data-generating mechanisms,
  simulation settings, empirical standardization, and performance summaries.
- `Simulated_Case_Studies_gibbs_helpers.R`: contains conditional updates and
  posterior estimand calculations.
- `CC_Sim_Covariates.csv`: paired age and MDS ADL covariate input used in
  all six case studies.

## Run

See the [repository README](../README.md) for packages. From the repository root:

```sh
Rscript --vanilla Simulated_Case_Studies/Simulated_Case_Studies_main.R case1
```

If no case is supplied, the main script uses `manual_case_study` (initially
`case2`). In RStudio, set that value and source the main script. The default
settings use 10,000 complete datasets to approximate truth and 2,500
datasets for operating characteristics. Within each fitted dataset, the
sampler discards 200 iterations and retains every fifth subsequent draw
for 1,000 posterior draws. For a short example, reduce `truth_runs`,
`iterations`, and the MCMC settings in `make_simulation_config()`.

## Check the appendix formulas without fitting the model

```sh
Rscript --vanilla Simulated_Case_Studies/Simulated_Case_Studies_main.R case1 --truth-only
```

This compares the analytic expressions with Monte Carlo approximations
obtained by pooling outcome contrasts and eligible-resident counts across
complete datasets. The check is available for Cases 1, 2, 5, and 6; the
script does not use these closed-form formulas for the alternative link
or outcome interactions in Cases 3 and 4.

## Reproducibility and outputs

The simulation seed is 1. Within each case, one covariate pool is drawn
with replacement and held fixed. Its rows are permuted across resident
positions in each dataset. Operating-characteristic replication `m` uses
seed `1 + m`. Simulation labels are ordered by the first mixture-mean
coordinate in ascending order. Super-population standardization uses the
pooled empirical covariate distribution.

The script writes results to `outputs/<case>/`: truth, replication-level
posterior summaries and performance metrics, uncertainty ratios, and an
RDS with the run configuration, priors, seeds, and results. The median SD
ratio divides the median posterior SD by the empirical SD of posterior
means across replications. When available, the formula comparison is
written to `<case>_truth_formula_check.csv`. A new run of the same case
replaces files for that case.
