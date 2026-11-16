# Binary outcome analysis main script for the synthetic METRIcAL example.
#
# This script runs the multichain analysis for the binary
# antipsychotic-use outcome. The model combines a latent implementation-strata
# mixture model, a multilevel probit compliance model, and a multilevel probit
# outcome model. Super-population ITT and CACE estimands are evaluated by
# averaging integrated probit risks at each retained posterior draw.
#
# Associated files:
# - `Binary_Outcome_Analysis_model_helpers.R` contains the Gibbs updating
#   functions, latent probit augmentation steps, imputation steps, and
#   analytic estimand calculations and posterior predictive simulation.
# - `Binary_Outcome_Analysis_helpers.R` contains the single-chain and
#   multichain analysis workflow, posterior relabeling, pooled summaries,
#   convergence diagnostics, and trace-plot utilities.
# Main analysis function:
# `run_binary_outcome_multichain(config, n_chains, chain_seeds,
# make_trace_pdf)` runs independent MCMC chains and returns chain-specific
# posterior draws, pooled estimand draws, pooled posterior summaries,
# Gelman--Rubin and effective-size diagnostics, and posterior predictive checks.
#
# Main outputs:
# - `outputs/Binary_Outcome_multichain_fit.rds` stores the full fitted
#   analysis object.
# - `outputs/Binary_Outcome_multichain_traces.pdf` stores trace plots
#   for the primary relabelled estimands when trace output is requested.

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
    script_path <- sub("^--file=", "", file_match[1])
    script_path <- gsub("~\\+~", " ", script_path)
    script_path <- path.expand(script_path)
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
source(file.path(script_dir, "Binary_Outcome_Analysis_helpers.R"))

output_dir <- file.path(script_dir, "outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

config <- make_binary_outcome_config(
  output_dir = output_dir
)

fit <- run_binary_outcome_multichain(
  config = config,
  n_chains = 4L,
  make_trace_pdf = TRUE
)

saveRDS(
  fit,
  file.path(output_dir, "Binary_Outcome_multichain_fit.rds")
)

export_binary_outcome_tables(fit)

print(fit$pooled_posterior_summary)
print(fit$diagnostics)
print(fit$pp_check_results)
print(fit$expanded_ppc_results)
