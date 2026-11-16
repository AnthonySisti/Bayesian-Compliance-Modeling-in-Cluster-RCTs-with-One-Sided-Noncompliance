#
# This script is the main script for the continuous-outcome analysis. It
# runs the multichain Gibbs sampler, writes the fitted object to disk, and
# prints the principal posterior summaries, MCMC diagnostics, and posterior
# predictive checks used to assess model fit.
#
# Files associated with this analysis:
# - `Continuous_Outcome_Analysis_main.R`
#   Main execution script. It defines the analysis configuration, launches the
#   multichain fit, and saves the fitted object and printed summaries.
# - `Continuous_Outcome_Analysis_helpers.R`
#   Workflow-level helper file. It defines the configuration object, the
#   single-chain and multichain analysis functions, post hoc relabeling of the
#   latent strata, posterior summaries, diagnostics, and posterior predictive
#   checks.
# - `Continuous_Outcome_Analysis_model_helpers.R`
#   Gibbs-step helper file. It contains the conditional-update functions used
#   within each Gibbs iteration.
#
# Data inputs expected by the helper functions:
# - `CZ_Matrix_Synth_OG.csv`
#   Facility-level compliance and covariate summaries.
# - `X_Synth_OG.csv`
#   Individual-level baseline covariates.
# - `Outcome_Synth_Continuous_OG.csv`
#   Individual-level treatment, compliance, and continuous outcome data.
#
# Main outputs produced by this script:
# - `outputs/Continuous_Outcome_multichain_fit.rds`
#   Full fitted object containing chain-specific retained draws, pooled
#   posterior summaries, diagnostics, and posterior predictive checks.
# - `outputs/Continuous_Outcome_multichain_traces.pdf`
#   Trace plots for the primary ITT and CACE estimands across chains.
#
# The central analysis function is `run_continuous_outcome_multichain()`.
# Inputs:
# - `config`: Gibbs-sampler settings, including the random seed, burn-in,
#   thinning, number of retained posterior draws, and output location.
# - `n_chains`: number of independently initialized Markov chains.
# - `chain_seeds`: chain-specific seeds for reproducible multichain execution.
# - `make_trace_pdf`: indicator for whether to save trace plots for the primary
#   estimands.
# Outputs:
# - `chain_results`: retained posterior output for each chain.
# - `pooled_estimand_draws` and `pooled_posterior_summary`: combined-chain
#   inference for the primary causal estimands.
# - `diagnostics`: Gelman-Rubin statistics and effective sample sizes for the
#   primary estimands.
# - `pp_check_results`: posterior predictive checks summarizing whether the
#   fitted model reproduces the observed signal and noise structure within the
#   latent compliance strata.
#
# Substantively, these outputs summarize the posterior distribution of the
# super-population treatment and complier effects, assess between-chain
# convergence, and evaluate whether the fitted model can reproduce the main
# empirical features targeted by the application.

get_script_dir <- function() {
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    path <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(path)) {
      return(dirname(normalizePath(path)))
    }
  }

  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_match <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_match) > 0) {
    script_path <- gsub("~\\+~", " ", sub("^--file=", "", file_match[1]))
    if (file.exists(script_path)) {
      return(dirname(normalizePath(script_path)))
    }
  }

  frame_files <- vapply(sys.frames(), function(frame) {
    if (!is.null(frame$ofile)) normalizePath(frame$ofile) else NA_character_
  }, character(1))
  frame_files <- frame_files[!is.na(frame_files)]
  if (length(frame_files) > 0) {
    return(dirname(frame_files[length(frame_files)]))
  }

  stop("Could not determine the script directory.")
}

script_dir <- get_script_dir()
source(file.path(script_dir, "Continuous_Outcome_Analysis_helpers.R"))

config <- make_continuous_outcome_config(
  seed = 20260422L,
  burn = 500L,
  thin = 3L,
  n_imps = 1000L,
  K = 2L,
  plot_initial_clusters = FALSE,
  output_dir = file.path(script_dir, "outputs")
)

fit <- run_continuous_outcome_multichain(
  config = config,
  n_chains = 4L,
  make_trace_pdf = TRUE
)

saveRDS(
  fit,
  file.path(config$output_dir, "Continuous_Outcome_multichain_fit.rds")
)

export_continuous_outcome_tables(fit)

print(fit$pooled_posterior_summary)
print(fit$diagnostics)
print(fit$pp_check_results)
print(fit$expanded_ppc_results)
