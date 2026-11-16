############################################################################
## Simulated case studies described in the
## supplementary materials.
##
## For the selected case study, this script:
## 1. approximates the super-population truth by repeated simulation under the
##    data-generating mechanism,
## 2. fits the Bayesian model to each replicated dataset by Gibbs sampling,
## 3. computes posterior draws of the target estimands, and
## 4. summarizes operating characteristics across simulation replications.
############################################################################

# Locate the companion scripts and the resident covariate input.
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
    path <- gsub("~\\+~", " ", sub("^--file=", "", file_match[1]))
    return(dirname(normalizePath(path, mustWork = TRUE)))
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
setwd(script_dir)
source("Simulated_Case_Studies_utils.R")
source("Simulated_Case_Studies_gibbs_helpers.R")
load_case_study_packages()

plot_initial_clusters <- interactive()

# Select the case-study setting to run when launching the script interactively.
# If the script is run from the command line with a case-study argument, that
# argument overrides this value.
manual_case_study <- "case2"
run_args <- commandArgs(trailingOnly = TRUE)
selected_case_study <- if (length(run_args) >= 1 && nzchar(run_args[1])) {
  run_args[1]
} else {
  manual_case_study
}
config <- make_simulation_config(case_study = selected_case_study)
truth_only <- "--truth-only" %in% run_args[-1L]
if (truth_only && !supports_analytic_truth_check(config)) {
  stop("The analytic truth comparison is available for Cases 1, 2, 5, and 6.")
}
simulation_seed <- config$simulation_seed

############################################################################
## Simulation settings
############################################################################

# The outer loop indexes simulated datasets used to evaluate operating
# characteristics. Within each dataset, the Gibbs sampler runs for `Sims`
# iterations, of which `n.imps` draws are retained after burn-in and thinning.
iterations <- config$iterations
truth_runs <- config$truth_runs
n.imps <- config$n_imps
thin <- config$thin
burn <- config$burn
Sims <- burn + thin * n.imps
samps_use <- burn + 1:n.imps * thin
samp_lookup <- integer(Sims)
samp_lookup[samps_use] <- seq_along(samps_use)

############################################################################
## Empirical covariate distribution and study layout
############################################################################

# The resident-level baseline covariates are drawn from the empirical
# distribution in the METRIcAL application and then reused across all
# simulation replications.
Covariates <- load_simulation_covariates()

J <- config$J
K <- config$K
pi_k <- config$pi_k
P <- config$P
L <- config$L
YCN <- ncol(Covariates)
DCN <- ncol(Covariates)
N <- sum(simulation_cluster_sizes(config))

# Replication-level summaries for the principal estimands reported in the
# manuscript tables.
ITTPostMean <- rep(NA_real_, iterations)
ITTPostSD <- rep(NA_real_, iterations)
ITTCover <- rep(NA_real_, iterations)
ITTBias <- rep(NA_real_, iterations)
ITTIW <- rep(NA_real_, iterations)

CACEPostMean <- rep(NA_real_, iterations)
CACEPostSD <- rep(NA_real_, iterations)
CACECover <- rep(NA_real_, iterations)
CACEBias <- rep(NA_real_, iterations)
CACEIW <- rep(NA_real_, iterations)

ITT1PostMean <- rep(NA_real_, iterations)
ITT1PostSD <- rep(NA_real_, iterations)
ITT1Cover <- rep(NA_real_, iterations)
ITT1Bias <- rep(NA_real_, iterations)
ITT1IW <- rep(NA_real_, iterations)

ITT2PostMean <- rep(NA_real_, iterations)
ITT2PostSD <- rep(NA_real_, iterations)
ITT2Cover <- rep(NA_real_, iterations)
ITT2Bias <- rep(NA_real_, iterations)
ITT2IW <- rep(NA_real_, iterations)

CACE1PostMean <- rep(NA_real_, iterations)
CACE1PostSD <- rep(NA_real_, iterations)
CACE1Cover <- rep(NA_real_, iterations)
CACE1Bias <- rep(NA_real_, iterations)
CACE1IW <- rep(NA_real_, iterations)

CACE2PostMean <- rep(NA_real_, iterations)
CACE2PostSD <- rep(NA_real_, iterations)
CACE2Cover <- rep(NA_real_, iterations)
CACE2Bias <- rep(NA_real_, iterations)
CACE2IW <- rep(NA_real_, iterations)

# Case-study indicators determine which data-generating mechanism is active.
skew <- config$skew
mix.df <- config$mix_df
useBurr <- config$use_burr
interaction_coeff_treat <- config$interaction_coeff_treat
interaction_coeff_control <- config$interaction_coeff_control

# Draw a finite covariate population from the empirical distribution once, then
# resample from this finite population within each simulated dataset.
set.seed(simulation_seed)
X_full <- sample_covariate_population(Covariates, config)

static_structure <- build_static_simulation_structure(config)
Facility <- static_structure$facility
treat_Fac_bin <- static_structure$treat_fac_bin
OutcomeData <- static_structure$outcome_data
treat_Fac_bool <- treat_Fac_bin == 1
treat_Fac_ind <- which(treat_Fac_bin == 1)
treat_Person_bin <- OutcomeData$Treat
treat_Person_bool <- treat_Person_bin == 1
treat_person_ind <- which(treat_Person_bin == 1)

############################################################################
## Monte Carlo approximation to the super-population truth
############################################################################

# Pool outcome contrasts and eligible residents across replications before
# taking ratios to approximate the resident-level population estimands.
estimand_names <- c("ITT", "ITT1", "ITT2", "CACE", "CACE1", "CACE2")
truth_sums <- truth_counts <- matrix(NA_real_, nrow = truth_runs, ncol = 6L,
                                    dimnames = list(NULL, estimand_names))

# Repeatedly simulate complete datasets under the data-generating mechanism and
# accumulate contributions from the realized potential outcomes.
for (t in seq_len(truth_runs)) {
  dgm_sample <- generate_case_study_dataset(X_full, static_structure, config)
  OutcomeData <- dgm_sample$OutcomeData
  totals <- summarize_truth_totals(OutcomeData)
  truth_sums[t, ] <- totals["contrast_sum", ]
  truth_counts[t, ] <- totals["count", ]

  if (t %% 1000 == 0) {
    print(t)
  }
}

pooled_counts <- colSums(truth_counts)
if (any(pooled_counts == 0)) {
  stop("No eligible residents in the truth simulations for: ",
       paste(names(pooled_counts)[pooled_counts == 0], collapse = ", "))
}
pooled_truth <- colSums(truth_sums) / pooled_counts
TrueITT <- pooled_truth[["ITT"]]
TrueITT1 <- pooled_truth[["ITT1"]]
TrueITT2 <- pooled_truth[["ITT2"]]
TrueCACE <- pooled_truth[["CACE"]]
TrueCACE1 <- pooled_truth[["CACE1"]]
TrueCACE2 <- pooled_truth[["CACE2"]]

TrueEstimandData <- build_true_estimand_table(
  case_study = config$case_study,
  case_label = config$case_label,
  true_itt = TrueITT,
  true_itt1 = TrueITT1,
  true_itt2 = TrueITT2,
  true_cace = TrueCACE,
  true_cace1 = TrueCACE1,
  true_cace2 = TrueCACE2
)

# For case studies with closed-form truth expressions, compare the Monte Carlo
# truth approximation to the corresponding analytic formulas.
TruthFormulaCheck <- NULL
if (supports_analytic_truth_check(config)) {
  monte_carlo_truth <- c(
    ITT = TrueITT,
    ITT1 = TrueITT1,
    ITT2 = TrueITT2,
    CACE = TrueCACE,
    CACE1 = TrueCACE1,
    CACE2 = TrueCACE2
  )

  TruthFormulaCheck <- compare_monte_carlo_and_analytic_truth(
    x_population = X_full,
    config = config,
    dgm_sample = dgm_sample,
    monte_carlo_truth = monte_carlo_truth
  )
}

# The appendix comparison can be run independently of posterior sampling.
if (truth_only) {
  print(TruthFormulaCheck, row.names = FALSE)
  output_dir <- file.path(script_dir, "outputs", config$case_study)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(
    TruthFormulaCheck,
    file.path(output_dir, paste0(config$case_study, "_truth_formula_check.csv")),
    row.names = FALSE
  )
  quit(save = "no")
}

############################################################################
## Operating-characteristic simulation
############################################################################

# The static study layout is held fixed across replications; only the latent
# variables, compliance indicators, and outcomes are regenerated.
OutcomeData <- static_structure$outcome_data

interval_columns <- as.vector(rbind(paste0(estimand_names, "Lower"),
                                    paste0(estimand_names, "Upper")))
credible_intervals <- matrix(NA_real_, nrow = iterations,
                             ncol = length(interval_columns),
                             dimnames = list(NULL, interval_columns))

for (m in seq_len(iterations)) {
  # Each simulated dataset is generated under a deterministic seed so the full
  # operating-characteristic study is reproducible.
  set.seed(simulation_seed + m)
  dgm_sample <- generate_case_study_dataset(X_full, static_structure, config)
  X <- dgm_sample$X
  OutcomeData <- dgm_sample$OutcomeData
  Facility_Strata <- dgm_sample$Facility_Strata
  True_CZ <- dgm_sample$True_CZ
  True_C <- dgm_sample$True_C
  True_Z <- dgm_sample$True_Z
  Coeff_Alpha <- dgm_sample$Coeff_Alpha
  PhiD_Fac <- dgm_sample$PhiD_Fac
  Tau2_PhiD <- dgm_sample$Tau2_PhiD
  Coeff_BetaD0 <- dgm_sample$Coeff_BetaD0
  Coeff_BetaD1 <- dgm_sample$Coeff_BetaD1
  Coeff_D0 <- dgm_sample$Coeff_D0
  Coeff_D1 <- dgm_sample$Coeff_D1
  Coeff_muY <- dgm_sample$Coeff_muY
  PhiY_Fac <- dgm_sample$PhiY_Fac
  Tau2_PhiY <- dgm_sample$Tau2_PhiY
  Sigma2Y <- dgm_sample$Sigma2Y
  Sigma <- dgm_sample$Sigma
  Facility_df <- data.frame(Facility = Facility)

  CZ_Matrix <- True_CZ
  CZ_Matrix[treat_Fac_bin == 0, 1:L] <- NA
  C_Matrix <- as.matrix(CZ_Matrix[, 1:L, drop = FALSE])
  Z_Matrix <- as.matrix(CZ_Matrix[, (L + 1):(L + P), drop = FALSE])

  ######################################################################
  ## Posterior computation for a single simulated dataset
  ######################################################################

  C <- C_Matrix
  Z <- Z_Matrix
  K <- 2
  D <- OutcomeData$D
  X.use <- X
  Y1 <- OutcomeData$Y1
  Y0 <- OutcomeData$Y0
  Facility <- as.numeric(OutcomeData$Facility)
  L <- ncol(C_Matrix)
  P <- ncol(Z_Matrix)
  if (is.null(L)) {
    L <- 1
  }
  if (is.null(P)) {
    P <- 1
  }
  J <- nrow(CZ_Matrix)

  # The latent implementation strata are initialized by fitting a Gaussian
  # mixture model to the treated facilities, where the observed cluster-level
  # implementation summaries and baseline cluster covariates enter jointly.
  CZ <- CZ_Matrix[treat_Fac_bin == 1, , drop = FALSE]
  tot_S <- rep(0, K)
  while (min(tot_S) < 3) {
    check_mix <- try(mvnormalmixEM(CZ, arbvar = FALSE, k = K))
    init_clust_treat <- apply(check_mix$posterior, 1, which.max)
    init_clust <- rep(NA_integer_, J)
    init_clust[treat_Fac_bin == 1] <- init_clust_treat

    for (k in seq_len(K)) {
      tot_S[k] <- sum(init_clust_treat == k)
    }
  }

  # Control facilities do not have observed implementation summaries, so they
  # are assigned an initial latent stratum by proximity in the observed
  # baseline cluster covariate.
  cont_fac <- which(treat_Fac_bin == 0)
  for (i in seq_along(cont_fac)) {
    mean1 <- abs(CZ_Matrix[cont_fac[i], 2] - check_mix$mu[[1]][2])
    mean2 <- abs(CZ_Matrix[cont_fac[i], 2] - check_mix$mu[[2]][2])
    init_clust[cont_fac[i]] <- ifelse(mean1 < mean2, 1L, 2L)
  }

  if (plot_initial_clusters) {
    pairs(CZ, col = apply(check_mix$posterior, 1, which.max))
  }

  # Storage objects for posterior states and retained draws within this dataset.
  mu_Sim <- lapply(seq_len(K), function(k) matrix(NA_real_, nrow = L + P, ncol = Sims))
  Sigma_Sim <- vector("list", Sims)
  C_Sim <- lapply(seq_len(L), function(l) matrix(NA_real_, nrow = J, ncol = Sims))
  for (l in seq_len(L)) {
    C_Sim[[l]][treat_Fac_bin == 1, ] <- C_Matrix[treat_Fac_bin == 1, l]
  }
  S_Sim <- matrix(NA_integer_, nrow = J, ncol = Sims)
  pi_Sim <- matrix(NA_real_, nrow = K, ncol = Sims)

  ######################################################################
  ## Initialization of latent variables and model parameters
  ######################################################################
  S_Sim[, 1] <- init_clust
  for (k in seq_len(K)) {
    mu_Sim[[k]][, 1] <- colMeans(CZ_Matrix[treat_Fac_bin == 1 & S_Sim[, 1] == k, , drop = FALSE])
  }

  Sigma_Sim[[1]] <- diag(0, nrow = L + P)
  for (k in seq_len(K)) {
    Sigma_Sim[[1]] <- Sigma_Sim[[1]] + cov(CZ_Matrix[treat_Fac_bin == 1 & S_Sim[, 1] == k, , drop = FALSE])
  }

  for (l in seq_len(L)) {
    for (k in seq_len(K)) {
      C_Sim[[l]][treat_Fac_bin == 0 & init_clust == k, 1] <- mean(C_Matrix[treat_Fac_bin == 1 & init_clust == k, l])
    }
  }

  pi_Sim[, 1] <- 1 / K

  # The individual-level compliance model is a probit regression on baseline
  # covariates, fit separately within each latent implementation stratum.
  X_D <- build_compliance_design_matrix(X)

  # These objects store latent compliance indicators, missing potential
  # outcomes, and model parameters across Gibbs iterations.
  D_Sim <- matrix(NA_integer_, nrow = N, ncol = Sims)
  D_Sim[treat_Person_bool, 1] <- D[treat_Person_bool]
  for (s in seq_len(Sims)) {
    D_Sim[treat_Person_bool, s] <- D[treat_Person_bool]
  }

  facility_cache <- build_facility_index_cache(Facility_df$Facility, J)
  provider_Fac <- facility_cache$provider_fac
  person_Fac <- facility_cache$person_fac
  person_Fac_factor <- facility_cache$person_fac_factor
  treat_Fac_bool <- seq_len(J) %in% unique(OutcomeData$Facility[OutcomeData$Treat == 1])
  treat_Fac_ind <- as.numeric(unique(OutcomeData$Facility[OutcomeData$Treat == 1]))

  Y0_Sim <- matrix(NA_real_, nrow = N, ncol = Sims)
  Y1_Sim <- matrix(NA_real_, nrow = N, ncol = Sims)

  muY_Sim <- matrix(NA_real_, nrow = K, ncol = Sims)
  PhiY_Sim <- matrix(NA_real_, nrow = J, ncol = Sims)
  Sigma2_Sim <- matrix(NA_real_, nrow = K, ncol = Sims)
  Tau2_PhiY_Sim <- rep(NA_real_, Sims)
  PhiY_Person <- rep(0, nrow(X))

  Y1_Sim[treat_Person_bool, 1] <- Y1[treat_Person_bool]
  Y0_Sim[!treat_Person_bool, 1] <- Y0[!treat_Person_bool]
  for (s in seq_len(Sims)) {
    Y1_Sim[treat_Person_bool, s] <- Y1[treat_Person_bool]
    Y0_Sim[!treat_Person_bool, s] <- Y0[!treat_Person_bool]
  }

  PhiY_Sim[, 1] <- 0
  Tau2_PhiY_Sim[1] <- 1
  Sigma2_Sim[, 1] <- 1

  PhiD_Sim <- matrix(0, nrow = J, ncol = Sims)
  Tau2_PhiD_Sim <- rep(0.1, Sims)
  PhiD_Person <- rep(0, nrow(X_D))

  Strata_Fac <- S_Sim[, 1]
  Strata_Person <- Strata_Fac[Facility]

  Alpha_Sim <- lapply(seq_len(K), function(k) matrix(NA_real_, nrow = ncol(X_D), ncol = Sims))

  # Initialize the compliance coefficients from stratum-specific probit models
  # fit to the treated-arm data, and use these to impute latent compliance for
  # control-arm individuals.
  CheckData_comp <- as.data.frame(cbind(D = D_Sim[, 1], X_D))
  CheckData_comp$D[!treat_Person_bool] <- NA
  prob_D <- rep(NA_real_, length(D))
  for (k in seq_len(K)) {
    reg_compk <- glm(D ~ . - 1, data = CheckData_comp[Strata_Person == k, , drop = FALSE], family = binomial("probit"))
    Alpha_Sim[[k]][, 1] <- summary(reg_compk)$coefficients[, 1]
    prob_D[Strata_Person == k] <- predict(reg_compk, newdata = CheckData_comp[Strata_Person == k, , drop = FALSE], type = "response")
  }

  D_Sim[!treat_Person_bool, 1] <- rbinom(sum(!treat_Person_bool), 1, prob_D[!treat_Person_bool])

  # Initialize the Albert-Chib latent Gaussian variables underlying the probit
  # compliance model.
  U_Sim <- matrix(NA_real_, nrow = nrow(X_D), ncol = Sims)
  U_Sim[D_Sim[, 1] == 1, 1] <- truncnorm::rtruncnorm(sum(D_Sim[, 1]), a = 0, b = Inf)
  U_Sim[D_Sim[, 1] == 0, 1] <- truncnorm::rtruncnorm(sum(1 - D_Sim[, 1]), a = -Inf, b = 0)
  X_D <- as.matrix(X_D)

  # Build the design matrices for the outcome model under observed treatment,
  # under treatment for all individuals, and under control for all individuals.
  outcome_design_cache <- build_outcome_design_cache(X, OutcomeData$Treat)
  outcome_designs <- update_outcome_designs_from_cache(outcome_design_cache, D_Sim[, 1])
  X_Main_obs <- outcome_designs$observed
  X_Main_treat <- outcome_designs$treat_all
  X_Main_control <- outcome_designs$control_current
  X_Main_control_comp <- outcome_designs$control_comp
  X_Main_control_nocomp <- outcome_designs$control_nocomp

  YCN <- ncol(X)
  BetaD0_Sim <- lapply(seq_len(K), function(k) matrix(NA_real_, nrow = YCN, ncol = Sims))
  BetaD1_Sim <- lapply(seq_len(K), function(k) matrix(NA_real_, nrow = YCN, ncol = Sims))
  Delta0_Sim <- matrix(NA_real_, nrow = K, ncol = Sims)
  Delta1_Sim <- matrix(NA_real_, nrow = K, ncol = Sims)
  Beta_full_Sim <- lapply(seq_len(K), function(k) matrix(NA_real_, nrow = ncol(X_Main_obs), ncol = Sims))

  Y_obs <- build_observed_outcomes(Y0, Y1, OutcomeData$Treat)
  CheckData <- as.data.frame(cbind(Y_obs, X_Main_obs))
  colnames(CheckData)[1] <- "Y"

  # Initialize the stratum-specific outcome regression coefficients from linear
  # models fit to the observed outcomes.
  for (k in seq_len(K)) {
    reg1k <- lm(Y ~ . - 1, data = CheckData[Strata_Person == k, , drop = FALSE])
    muY_Sim[k, 1] <- summary(reg1k)$coefficients[1, 1]
    BetaD0_Sim[[k]][, 1] <- summary(reg1k)$coefficients[1 + seq_len(nrow(BetaD0_Sim[[1]])), 1]
    Delta0_Sim[k, 1] <- summary(reg1k)$coefficients[2 + nrow(BetaD0_Sim[[1]]), 1]
    BetaD1_Sim[[k]][, 1] <- summary(reg1k)$coefficients[(2 + YCN) + seq_len(YCN), 1]
    Delta1_Sim[k, 1] <- summary(reg1k)$coefficients[ncol(X_Main_obs), 1]
    Y1_Sim[Strata_Person == k & treat_Person_bool == 0, 1] <- predict(reg1k, newdata = as.data.frame(X_Main_treat[Strata_Person == k & treat_Person_bool == 0, , drop = FALSE]))
    Y0_Sim[Strata_Person == k & treat_Person_bool == 1, 1] <- predict(reg1k, newdata = as.data.frame(X_Main_control[Strata_Person == k & treat_Person_bool == 1, , drop = FALSE]))
    Beta_full_Sim[[k]][, 1] <- c(
      muY_Sim[k, 1],
      BetaD0_Sim[[k]][, 1],
      Delta0_Sim[k, 1],
      BetaD1_Sim[[k]][, 1],
      Delta1_Sim[k, 1]
    )
  }

  ######################################################################
  ## Gibbs sampling
  ######################################################################
  Y1_use <- matrix(NA_real_, nrow = N, ncol = n.imps)
  Y0_use <- matrix(NA_real_, nrow = N, ncol = n.imps)
  D_Samps_use <- matrix(NA_real_, nrow = N / 2, ncol = n.imps)
  Beta_full_use1 <- matrix(NA_real_, nrow = nrow(Beta_full_Sim[[1]]), ncol = n.imps)
  Beta_full_use2 <- matrix(NA_real_, nrow = nrow(Beta_full_Sim[[2]]), ncol = n.imps)
  ITT_use <- rep(NA_real_, n.imps)
  CACE_use <- rep(NA_real_, n.imps)
  ITT1_use <- rep(NA_real_, n.imps)
  ITT2_use <- rep(NA_real_, n.imps)
  CACE1_use <- rep(NA_real_, n.imps)
  CACE2_use <- rep(NA_real_, n.imps)
  PhiY_use <- matrix(NA_real_, nrow = J, ncol = n.imps)
  PhiD_use <- matrix(NA_real_, nrow = J, ncol = n.imps)
  D_use <- matrix(NA_real_, nrow = length(D), ncol = n.imps)

  # Prior hyperparameters for the latent-strata model, compliance model,
  # outcome model, and cluster-level random effects.
  V_0_inv <- diag(1 / 100, nrow = L + P)
  m_0 <- matrix(0, nrow = L + P, ncol = 1)
  nu_0 <- 7
  S_0 <- diag(1 / 100, nrow = L + P)
  lambda_pk <- 3
  m_alpha <- rep(0, ncol(X_D))
  v_alpha_inv <- diag(1 / 100, ncol(X_D))
  m_beta <- rep(0, ncol(X_Main_obs))
  v_beta_inv <- diag(1 / 100, ncol(X_Main_obs))
  v_tauD_mu <- 25
  m_tauD_mu <- 0
  a_sigmas <- 1
  b_sigmas <- 1
  max_tauPhiD <- 10
  max_tauPhiY <- 225
  alt_phiD <- matrix(0, nrow = J, ncol = K)
  alt_phiY <- matrix(0, nrow = J, ncol = K)

  CZ_obs <- CZ_Matrix[treat_Fac_bool, , drop = FALSE]
  C_obs <- as.matrix(CZ_Matrix[treat_Fac_bool, 1:L, drop = FALSE])
  Z_obs <- as.matrix(CZ_Matrix[treat_Fac_bool, (L + 1):(L + P), drop = FALSE])
  treat_index <- which(treat_Fac_bool)
  J_treat <- sum(treat_Fac_bool)
  control_index <- which(!treat_Fac_bool)
  J_control <- sum(!treat_Fac_bool)

  S_use <- matrix(NA_integer_, nrow = J, ncol = n.imps)
  C_use <- lapply(seq_len(L), function(l) matrix(NA_real_, nrow = J, ncol = n.imps))

  Coeff_Alpha_est <- matrix(NA_real_, nrow = nrow(Alpha_Sim[[1]]), ncol = K)
  Tau2_PhiD_est <- NA_real_
  pi_est <- rep(NA_real_, K)
  Coeff_BetaD1_est <- matrix(NA_real_, nrow = nrow(BetaD1_Sim[[1]]), ncol = K)
  Coeff_D1_est <- rep(NA_real_, K)
  Coeff_D0_est <- rep(NA_real_, K)

  Strata_Fac <- S_Sim[, 1]
  Strata_Person <- S_Sim[person_Fac, 1]
  PhiY_Person <- PhiY_Sim[person_Fac, 1]
  PhiD_Person <- PhiD_Sim[person_Fac, 1]

  CZ_i <- CZ_Matrix
  C_i <- C_Matrix
  Z_i <- Z_Matrix

  for (i in seq_len(Sims - 1)) {
    new <- i + 1

    # Update the current cluster-level latent data object by inserting the most
    # recent imputations for the unobserved control-arm implementation summary.
    for (l in seq_len(L)) {
      CZ_i[treat_Fac_bin == 0, l] <- C_Sim[[l]][treat_Fac_bin == 0, i]
      C_i[treat_Fac_bin == 0, l] <- C_Sim[[l]][treat_Fac_bin == 0, i]
    }

    alpha_current <- lapply(seq_len(K), function(k) Alpha_Sim[[k]][, i])
    beta_current <- lapply(seq_len(K), function(k) Beta_full_Sim[[k]][, i])

    # 1(a) Update the latent-stratum mean vectors for the joint distribution of
    # cluster-level implementation and baseline cluster covariates.
    new_mu <- update_mixture_component_means(
      cz_i = CZ_i,
      strata_fac = Strata_Fac,
      sigma_current = Sigma_Sim[[i]],
      V_0_inv = V_0_inv,
      m_0 = m_0,
      K = K
    )
    for (k in seq_len(K)) {
      mu_Sim[[k]][, new] <- new_mu[, k]
    }

    # 1(b) Update the shared covariance matrix in the latent-strata mixture.
    Sigma_Sim[[new]] <- update_mixture_covariance(
      c_i = C_i,
      z_i = Z_i,
      strata_fac = Strata_Fac,
      mu_draws = new_mu,
      S_0 = S_0,
      nu_0 = nu_0,
      K = K
    )

    # 1(c) Update the stratum membership probabilities.
    pi_Sim[, new] <- update_mixture_weights(
      strata_fac = Strata_Fac,
      lambda_pk = lambda_pk,
      K = K
    )

    # 2(a) Update the latent Gaussian variables for the probit compliance model.
    U_Sim[, new] <- update_latent_probit_variables(
      D_current = D_Sim[, i],
      strata_person = Strata_Person,
      X_D = X_D,
      alpha_current = alpha_current,
      phiD_person = PhiD_Person,
      K = K
    )

    # 2(c) Update the cluster-level random effects in the compliance model.
    PhiD_Sim[, new] <- update_compliance_random_effects(
      U_new = U_Sim[, new],
      X_D = X_D,
      alpha_current = alpha_current,
      strata_person = Strata_Person,
      strata_fac = Strata_Fac,
      provider_fac = provider_Fac,
      person_fac = person_Fac,
      tau2_prev = Tau2_PhiD_Sim[i],
      K = K
    )
    PhiD_Person <- PhiD_Sim[person_Fac, new]

    # 2(d) Update the variance of the cluster-level compliance random effects.
    Tau2_PhiD_Sim[new] <- update_shared_random_effect_variance(
      phi_draw = PhiD_Sim[, new],
      upper_bound_tau2 = max_tauPhiD,
      context = sprintf("%s, replication %d, Gibbs iteration %d, compliance",
                        config$case_study, m, new)
    )

    # 2(b) Update the stratum-specific compliance regression coefficients.
    alpha_new <- update_compliance_coefficients(
      X_D = X_D,
      U_new = U_Sim[, new],
      phiD_person = PhiD_Person,
      strata_person = Strata_Person,
      v_alpha_inv = v_alpha_inv,
      m_alpha = m_alpha,
      K = K
    )
    for (k in seq_len(K)) {
      Alpha_Sim[[k]][, new] <- alpha_new[[k]]
    }

    # 3(c) Update the cluster-level random effects in the outcome model.
    PhiY_Sim[, new] <- update_outcome_random_effects(
      Y_obs = Y_obs,
      X_Main_obs = X_Main_obs,
      beta_current = beta_current,
      strata_person = Strata_Person,
      strata_fac = Strata_Fac,
      provider_fac = provider_Fac,
      person_fac = person_Fac,
      sigma2_prev = Sigma2_Sim[, i],
      tau2_prev = Tau2_PhiY_Sim[i],
      K = K
    )
    PhiY_Person <- PhiY_Sim[person_Fac, new]

    # 3(d) Update the variance of the cluster-level outcome random effects.
    Tau2_PhiY_Sim[new] <- update_shared_random_effect_variance(
      phi_draw = PhiY_Sim[, new],
      upper_bound_tau2 = max_tauPhiY,
      context = sprintf("%s, replication %d, Gibbs iteration %d, outcome",
                        config$case_study, m, new)
    )

    # 3(a) Update the stratum-specific outcome regression coefficients.
    beta_new <- update_outcome_coefficients(
      Y_obs = Y_obs,
      X_Main_obs = X_Main_obs,
      phiY_person = PhiY_Person,
      strata_person = Strata_Person,
      sigma2_current = Sigma2_Sim[, i],
      v_beta_inv = v_beta_inv,
      m_beta = m_beta,
      K = K
    )
    for (k in seq_len(K)) {
      Beta_full_Sim[[k]][, new] <- beta_new[[k]]
      muY_Sim[k, new] <- Beta_full_Sim[[k]][1, new]
      BetaD0_Sim[[k]][, new] <- Beta_full_Sim[[k]][1 + seq_len(YCN), new]
      Delta0_Sim[k, new] <- Beta_full_Sim[[k]][2 + YCN, new]
      BetaD1_Sim[[k]][, new] <- Beta_full_Sim[[k]][(2 + YCN) + seq_len(YCN), new]
      Delta1_Sim[k, new] <- Beta_full_Sim[[k]][ncol(X_Main_obs), new]
    }

    # 3(b) Update the stratum-specific residual variances in the outcome model.
    Sigma2_Sim[, new] <- update_outcome_sigmas(
      Y_obs = Y_obs,
      X_Main_obs = X_Main_obs,
      beta_new = beta_new,
      phiY_person = PhiY_Person,
      strata_person = Strata_Person,
      a_sigmas = a_sigmas,
      b_sigmas = b_sigmas,
      K = K
    )

    # 4. Update latent implementation-stratum membership for each facility.
    strata_update <- update_strata_membership(
      cz_i = CZ_i,
      Y_obs = Y_obs,
      X_Main_obs = X_Main_obs,
      D_current = D_Sim[, i],
      X_D = X_D,
      mu_draws = new_mu,
      sigma_new = Sigma_Sim[[new]],
      beta_new = beta_new,
      sigma2_new = Sigma2_Sim[, new],
      alpha_new = alpha_new,
      phiY_person = PhiY_Person,
      phiD_person = PhiD_Person,
      pi_new = pi_Sim[, new],
      person_fac = person_Fac,
      K = K
    )
    S_Sim[, new] <- strata_update$strata_fac
    Strata_Fac <- S_Sim[, new]
    Strata_Person <- Strata_Fac[Facility]

    # 5. Impute the unobserved cluster-level implementation summaries in the
    # control arm from their conditional Gaussian distribution.
    C_i <- impute_control_compliance_measures(
      C_current = C_i,
      Z_matrix = Z_Matrix,
      mu_draws = new_mu,
      sigma_new = Sigma_Sim[[new]],
      strata_fac_new = Strata_Fac,
      treat_fac_bool = treat_Fac_bool,
      L = L,
      P = P,
      K = K
    )
    for (l in seq_len(L)) {
      C_Sim[[l]][, new] <- C_i[, l]
    }

    # 6. Impute latent compliance indicators for individuals in control
    # facilities from their Bernoulli full conditional.
    D_update <- impute_unobserved_compliance(
      D_current = D_Sim[, i],
      Y_obs = Y_obs,
      X_Main_control_comp = X_Main_control_comp,
      X_Main_control_nocomp = X_Main_control_nocomp,
      X_D = X_D,
      alpha_new = alpha_new,
      beta_new = beta_new,
      sigma2_new = Sigma2_Sim[, new],
      phiY_person = PhiY_Person,
      phiD_person = PhiD_Person,
      strata_person = Strata_Person,
      treat_person_bool = treat_Person_bool,
      K = K
    )
    D_samp_probs <- D_update$probabilities
    D_Sim[, new] <- D_update$D_new

    # Rebuild the outcome-model design matrices under the newly imputed
    # individual compliance indicators.
    outcome_designs <- update_outcome_designs_from_cache(outcome_design_cache, D_Sim[, new])
    X_Main_obs <- outcome_designs$observed
    X_Main_treat <- outcome_designs$treat_all
    X_Main_control <- outcome_designs$control_current
    X_Main_control_comp <- outcome_designs$control_comp
    X_Main_control_nocomp <- outcome_designs$control_nocomp

    if (i %% 500 == 0) {
      print(i)
      print(table(S_Sim[, i], treat_Fac_bin))
    }

    thin_ind <- samp_lookup[new]
    if (thin_ind > 0) {
      S_use[, thin_ind] <- S_Sim[, new]
      for (l in seq_len(L)) {
        C_use[[l]][, thin_ind] <- C_Sim[[l]][, new]
      }

      D_use[, thin_ind] <- D_Sim[, new]
      Y1_use[, thin_ind] <- Y1_Sim[, new]
      Y0_use[, thin_ind] <- Y0_Sim[, new]
      Beta_full_use1[, thin_ind] <- Beta_full_Sim[[1]][, new]
      Beta_full_use2[, thin_ind] <- Beta_full_Sim[[2]][, new]
      PhiY_use[, thin_ind] <- PhiY_Sim[, new]
      PhiD_use[, thin_ind] <- PhiD_Sim[, new]

      # Compute the posterior draw of the super-population ITT and CACE,
      # marginally and within the ordered latent implementation strata.
      estimand_draw <- compute_posterior_draw_estimands(
        X = X,
        X_D = X_D,
        alpha_draws = lapply(seq_len(K), function(k) Alpha_Sim[[k]][, new]),
        beta1_draws = lapply(seq_len(K), function(k) BetaD1_Sim[[k]][, new]),
        delta0_draws = Delta0_Sim[, new],
        delta1_draws = Delta1_Sim[, new],
        pi_draw = pi_Sim[, new],
        tau2_phi_d_draw = Tau2_PhiD_Sim[new],
        mu_draws = sapply(seq_len(K), function(k) mu_Sim[[k]][, new])
      )

      ITT_use[thin_ind] <- estimand_draw["ITT"]
      CACE_use[thin_ind] <- estimand_draw["CACE"]
      ITT1_use[thin_ind] <- estimand_draw["ITT1"]
      ITT2_use[thin_ind] <- estimand_draw["ITT2"]
      CACE1_use[thin_ind] <- estimand_draw["CACE1"]
      CACE2_use[thin_ind] <- estimand_draw["CACE2"]
    }
  }

  ######################################################################
  ## Replication-level posterior summaries
  ######################################################################

  # For each simulated dataset, summarize the retained posterior draws by bias,
  # interval coverage, and interval width relative to the super-population
  # truth values computed above.
  itt_perf <- summarize_estimand_performance(ITT_use, TrueITT)
  ITTPostMean[m] <- itt_perf["post_mean"]
  ITTPostSD[m] <- itt_perf["post_sd"]
  ITTBias[m] <- itt_perf["bias"]
  ITTCover[m] <- itt_perf["cover"]
  ITTIW[m] <- itt_perf["iw"]

  cace_perf <- summarize_estimand_performance(CACE_use, TrueCACE)
  CACEPostMean[m] <- cace_perf["post_mean"]
  CACEPostSD[m] <- cace_perf["post_sd"]
  CACEBias[m] <- cace_perf["bias"]
  CACECover[m] <- cace_perf["cover"]
  CACEIW[m] <- cace_perf["iw"]

  itt1_perf <- summarize_estimand_performance(ITT1_use, TrueITT1)
  ITT1PostMean[m] <- itt1_perf["post_mean"]
  ITT1PostSD[m] <- itt1_perf["post_sd"]
  ITT1Bias[m] <- itt1_perf["bias"]
  ITT1Cover[m] <- itt1_perf["cover"]
  ITT1IW[m] <- itt1_perf["iw"]

  itt2_perf <- summarize_estimand_performance(ITT2_use, TrueITT2)
  ITT2PostMean[m] <- itt2_perf["post_mean"]
  ITT2PostSD[m] <- itt2_perf["post_sd"]
  ITT2Bias[m] <- itt2_perf["bias"]
  ITT2Cover[m] <- itt2_perf["cover"]
  ITT2IW[m] <- itt2_perf["iw"]

  cace1_perf <- summarize_estimand_performance(CACE1_use, TrueCACE1)
  CACE1PostMean[m] <- cace1_perf["post_mean"]
  CACE1PostSD[m] <- cace1_perf["post_sd"]
  CACE1Bias[m] <- cace1_perf["bias"]
  CACE1Cover[m] <- cace1_perf["cover"]
  CACE1IW[m] <- cace1_perf["iw"]

  cace2_perf <- summarize_estimand_performance(CACE2_use, TrueCACE2)
  CACE2PostMean[m] <- cace2_perf["post_mean"]
  CACE2PostSD[m] <- cace2_perf["post_sd"]
  CACE2Bias[m] <- cace2_perf["bias"]
  CACE2Cover[m] <- cace2_perf["cover"]
  CACE2IW[m] <- cace2_perf["iw"]

  credible_intervals[m, ] <- unlist(lapply(
    list(itt_perf, itt1_perf, itt2_perf, cace_perf, cace1_perf, cace2_perf),
    function(x) x[c("lower", "upper")]
  ), use.names = FALSE)

  print(m)
  print(mean(ITTCover, na.rm = TRUE))
  print(mean(CACECover, na.rm = TRUE))
  print(mean(ITT1Cover, na.rm = TRUE))
  print(mean(ITT2Cover, na.rm = TRUE))
  print(mean(CACE1Cover, na.rm = TRUE))
  print(mean(CACE2Cover, na.rm = TRUE))
  print(mean(ITTBias, na.rm = TRUE))
  print(mean(CACEBias, na.rm = TRUE))
  print(mean(ITT1Bias, na.rm = TRUE))
  print(mean(ITT2Bias, na.rm = TRUE))
  print(mean(CACE1Bias, na.rm = TRUE))
  print(mean(CACE2Bias, na.rm = TRUE))

  if (m == 1 || m %% 10 == 0 || m == iterations) {
    EstimandData <- data.frame(
      ITTPostMean,
      ITTPostSD,
      ITTBias,
      ITTCover,
      ITTIW,
      ITT1PostMean,
      ITT1PostSD,
      ITT1Bias,
      ITT1Cover,
      ITT1IW,
      ITT2PostMean,
      ITT2PostSD,
      ITT2Bias,
      ITT2Cover,
      ITT2IW,
      CACEPostMean,
      CACEPostSD,
      CACEBias,
      CACECover,
      CACEIW,
      CACE1PostMean,
      CACE1PostSD,
      CACE1Bias,
      CACE1Cover,
      CACE1IW,
      CACE2PostMean,
      CACE2PostSD,
      CACE2Bias,
      CACE2Cover,
      CACE2IW
    )
  }
}

EstimandData <- cbind(EstimandData, credible_intervals)
EstimandData %>% round(3)

UncertaintyRatioData <- data.frame(
  estimand = c("ITT", "ITT1", "ITT2", "CACE", "CACE1", "CACE2"),
  rbind(
    summarize_uncertainty_ratio(ITTPostMean, ITTPostSD),
    summarize_uncertainty_ratio(ITT1PostMean, ITT1PostSD),
    summarize_uncertainty_ratio(ITT2PostMean, ITT2PostSD),
    summarize_uncertainty_ratio(CACEPostMean, CACEPostSD),
    summarize_uncertainty_ratio(CACE1PostMean, CACE1PostSD),
    summarize_uncertainty_ratio(CACE2PostMean, CACE2PostSD)
  ),
  row.names = NULL
)

UncertaintyRatioDataRounded <- UncertaintyRatioData
numeric_ratio_cols <- vapply(UncertaintyRatioDataRounded, is.numeric, logical(1))
UncertaintyRatioDataRounded[numeric_ratio_cols] <- lapply(
  UncertaintyRatioDataRounded[numeric_ratio_cols],
  round,
  digits = 3
)
UncertaintyRatioDataRounded

############################################################################
## Save reviewer-facing simulation summaries
############################################################################

# The objects saved here are the direct outputs needed to construct the
# simulation tables: the Monte Carlo truth, the optional analytic truth check,
# replicate-level operating characteristics, and posterior-uncertainty
# diagnostics across simulation replications.
output_dir <- file.path(script_dir, "outputs", config$case_study)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  TrueEstimandData,
  file.path(output_dir, paste0(config$case_study, "_true_estimands.csv")),
  row.names = FALSE
)

if (!is.null(TruthFormulaCheck)) {
  write.csv(
    TruthFormulaCheck,
    file.path(output_dir, paste0(config$case_study, "_truth_formula_check.csv")),
    row.names = FALSE
  )
}

write.csv(
  EstimandData,
  file.path(output_dir, paste0(config$case_study, "_estimand_performance_by_replication.csv")),
  row.names = FALSE
)

write.csv(
  UncertaintyRatioData,
  file.path(output_dir, paste0(config$case_study, "_uncertainty_ratios.csv")),
  row.names = FALSE
)

saveRDS(
  list(
    config = config,
    replication_seeds = simulation_seed + seq_len(iterations),
    priors = list(mixture_mean = m_0, mixture_mean_precision = V_0_inv,
                  mixture_df = nu_0, mixture_scale = S_0,
                  mixture_weight_concentration = rep_len(lambda_pk, K),
                  compliance_mean = m_alpha, compliance_precision = v_alpha_inv,
                  outcome_mean = m_beta, outcome_precision = v_beta_inv,
                  residual_shape = a_sigmas, residual_rate = b_sigmas,
                  compliance_sd_upper = sqrt(max_tauPhiD),
                  outcome_sd_upper = sqrt(max_tauPhiY)),
    truth_sums = truth_sums,
    truth_counts = truth_counts,
    TrueEstimandData = TrueEstimandData,
    TruthFormulaCheck = TruthFormulaCheck,
    EstimandData = EstimandData,
    UncertaintyRatioData = UncertaintyRatioData
  ),
  file.path(output_dir, paste0(config$case_study, "_simulation_results.rds"))
)
