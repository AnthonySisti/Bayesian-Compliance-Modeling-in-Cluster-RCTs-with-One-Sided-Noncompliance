############################################################################
## Utility functions for the simulation design and truth calculations.
##
## These functions define the case-study settings, generate complete datasets
## under the data-generating mechanism, and compute the estimands used to
## evaluate posterior operating characteristics.
############################################################################

# Load the packages required for simulation, posterior computation, and the
# parametric forms used in the latent implementation and compliance models.
load_case_study_packages <- function() {
  library("boot", include.only = "inv.logit")
  library(mixtools)
  library(coda)
  library(dplyr)
  library(MCMCpack, include.only = "riwish")
  library(DirichletReg, include.only = "rdirichlet")
  library(mvtnorm, include.only = "dmvnorm")
  library(MASS, include.only = "mvrnorm")
  library(scales, include.only = "alpha")
  library(invgamma)
  library(truncnorm)
  library(sn)
}

# Read the empirical resident-level baseline covariates used to define the
# finite covariate population from which simulated datasets are drawn.
load_simulation_covariates <- function() {
  covariate_path <- "CC_Sim_Covariates.csv"
  covariates <- read.csv(covariate_path, check.names = FALSE)

  expected_names <- c("age", "m3_adl")
  if (!identical(names(covariates), expected_names)) {
    stop(
      "Expected CC_Sim_Covariates.csv columns ",
      paste(expected_names, collapse = ", "),
      " but found ",
      paste(names(covariates), collapse = ", "),
      call. = FALSE
    )
  }

  if (!all(vapply(covariates, is.numeric, logical(1)))) {
    stop("CC_Sim_Covariates.csv must contain only numeric columns.", call. = FALSE)
  }

  covariates
}

# Specify the active simulation scenario. Each case study corresponds to one of
# the data-generating mechanisms described in the manuscript:
# - Case 1: correctly specified model
# - Case 2: misspecified latent implementation distribution
# - Case 3: misspecified individual compliance model
# - Case 4: misspecified outcome model
# - Case 5: correctly specified model with weak latent-stratum separation
# - Case 6: correctly specified model with unequal cluster sizes
make_case_study_spec <- function(case_study = "case1") {
  specs <- list(
    case1 = list(
      label = "Case Study 1: Correctly Specified DGM",
      mixture_mean_gap = 4,
      skew = c(0, 0),
      mix_df = Inf,
      use_burr = FALSE,
      interaction_coeff_treat = c(0, 0),
      interaction_coeff_control = c(0, 0)
    ),
    case2 = list(
      label = "Case Study 2: Misspecified Latent Implementation Strata Model",
      mixture_mean_gap = 4,
      skew = c(2, 2),
      mix_df = 5,
      use_burr = FALSE,
      interaction_coeff_treat = c(0, 0),
      interaction_coeff_control = c(0, 0)
    ),
    case3 = list(
      label = "Case Study 3: Incorrectly Specified Individual Compliance Model",
      mixture_mean_gap = 4,
      skew = c(0, 0),
      mix_df = Inf,
      use_burr = TRUE,
      interaction_coeff_treat = c(0, 0),
      interaction_coeff_control = c(0, 0)
    ),
    case4 = list(
      label = "Case Study 4: Misspecified Outcome Model",
      mixture_mean_gap = 4,
      skew = c(0, 0),
      mix_df = Inf,
      use_burr = FALSE,
      interaction_coeff_treat = c(-2, -2),
      interaction_coeff_control = c(-2, -2)
    ),
    case5 = list(
      label = "Case Study 5: Weak Latent-Stratum Separation",
      mixture_mean_gap = 1.5,
      skew = c(0, 0),
      mix_df = Inf,
      use_burr = FALSE,
      interaction_coeff_treat = c(0, 0),
      interaction_coeff_control = c(0, 0)
    )
  )

  # Keep the Case 1 mechanism; vary only the facility sizes.
  specs$case6 <- specs$case1
  specs$case6$label <- "Case Study 6: Unequal Cluster Sizes"

  if (!case_study %in% names(specs)) {
    stop("Unknown case study: ", case_study, call. = FALSE)
  }

  specs[[case_study]]
}

# Collect the fixed design parameters for the simulation study, including the
# number of facilities, individuals per facility, posterior sampling settings,
# and the case-study specific data-generating mechanism.
make_simulation_config <- function(case_study = "case1") {
  case_spec <- make_case_study_spec(case_study)

  list(
    case_study = case_study,
    case_label = case_spec$label,
    simulation_seed = 1L,
    iterations = 2500L,
    truth_runs = 10000L,
    n_imps = 1000L,
    thin = 5L,
    burn = 200L,
    J = 60L,
    Nj = if (case_study == "case6") rep(c(10L, 30L), 30L) else 20L,
    K = 2L,
    pi_k = c(0.5, 0.5),
    P = 1L,
    L = 1L,
    mixture_mean_gap = case_spec$mixture_mean_gap,
    skew = case_spec$skew,
    mix_df = case_spec$mix_df,
    use_burr = case_spec$use_burr,
    interaction_coeff_treat = case_spec$interaction_coeff_treat,
    interaction_coeff_control = case_spec$interaction_coeff_control
  )
}

# Alternative compliance link used only in the misspecified compliance case
# study.
burr_link <- function(x, c = 0.5) {
  1 - (1 + exp(x))^(-c)
}

# Create the treatment allocation and study indexing structure that is held
# fixed across simulation replications.
simulation_cluster_sizes <- function(config) {
  sizes <- config$Nj
  stopifnot(length(sizes) %in% c(1L, config$J),
            all(is.finite(sizes) & sizes > 0 & sizes == floor(sizes)))
  rep_len(as.integer(sizes), config$J)
}

build_static_simulation_structure <- function(config) {
  sizes <- simulation_cluster_sizes(config)
  n_total <- sum(sizes)
  facility <- rep(seq_len(config$J), times = sizes)
  treat_fac_bin <- ifelse(seq_len(config$J) <= config$J / 2, 0, 1)

  outcome_data <- data.frame(
    Person = seq_len(n_total),
    Facility = facility,
    Treat = treat_fac_bin[facility],
    Y0 = rep(NA_real_, n_total),
    Y1 = rep(NA_real_, n_total),
    D = rep(NA_real_, n_total),
    Strata = rep(NA_integer_, n_total)
  )

  list(
    n_total = n_total,
    facility = facility,
    treat_fac_bin = treat_fac_bin,
    outcome_data = outcome_data
  )
}

# Draw a finite resident-level covariate population from the empirical baseline
# covariate distribution.
sample_covariate_population <- function(covariates, config) {
  n_total <- sum(simulation_cluster_sizes(config))
  sampled_rows <- sample(seq_len(nrow(covariates)), n_total, replace = TRUE)
  x_full <- covariates[sampled_rows, , drop = FALSE]
  rownames(x_full) <- seq_len(nrow(x_full))
  as.matrix(x_full)
}

# Generate one complete dataset under the selected data-generating mechanism.
# The returned object includes both the observed data and the latent
# data-generating quantities needed for posterior evaluation.
generate_case_study_dataset <- function(x_full, static_structure, config) {
  J <- config$J
  K <- config$K
  L <- config$L
  P <- config$P
  N <- static_structure$n_total

  x_matrix <- x_full[sample(seq_len(nrow(x_full)), N, replace = FALSE), , drop = FALSE]
  rownames(x_matrix) <- seq_len(nrow(x_matrix))

  outcome_data <- static_structure$outcome_data
  facility_strata <- sample(seq_len(K), J, prob = config$pi_k, replace = TRUE)
  outcome_data$Strata <- facility_strata[outcome_data$Facility]

  mu_c <- matrix(NA_real_, nrow = L, ncol = K)
  mu_z <- matrix(NA_real_, nrow = P, ncol = K)
  base_mean_sequence <- 1 + config$mixture_mean_gap * seq_len(K)
  for (k in seq_len(K)) {
    mu_c[, k] <- base_mean_sequence[k]
    mu_z[, k] <- base_mean_sequence[k]
  }
  mu_c <- mu_c - mean(mu_c)
  mu_z <- mu_z - mean(mu_z)
  mu <- as.matrix(rbind(mu_c, mu_z))

  rho_sigma <- runif(1, -0.8, 0.8)
  sigma2_c <- runif(1, 0.5, 2)
  sigma2_z <- runif(1, 0.5, 2)
  sigma <- matrix(
    c(
      sigma2_c, rho_sigma * sqrt(sigma2_c * sigma2_z),
      rho_sigma * sqrt(sigma2_c * sigma2_z), sigma2_z
    ),
    nrow = 2,
    byrow = TRUE
  )

  true_cz <- matrix(NA_real_, nrow = J, ncol = L + P)
  for (k in seq_len(K)) {
    fac_in_k <- facility_strata == k
    true_cz[fac_in_k, ] <- sn::rmst(
      n = sum(fac_in_k),
      xi = mu[, k],
      Omega = sigma,
      alpha = config$skew,
      nu = config$mix_df
    )
  }

  x_d <- as.matrix(cbind(1, x_matrix))
  coeff_alpha <- matrix(NA_real_, nrow = ncol(x_d), ncol = K)
  phi_d_fac <- rep(NA_real_, J)
  tau2_phi_d <- 0.25
  for (k in seq_len(K)) {
    coeff_alpha[, k] <- c((k - 1) / 2, rep(-0.25 * k, ncol(x_d) - 1))
    phi_d_fac[facility_strata == k] <- rnorm(sum(facility_strata == k), 0, sqrt(tau2_phi_d))
  }

  phi_d_person <- phi_d_fac[outcome_data$Facility]
  outcome_data$PhiD <- phi_d_person

  compliance_linear_predictor <- diag(x_d %*% coeff_alpha[, outcome_data$Strata, drop = FALSE]) + phi_d_person
  d_probs <- if (config$use_burr) {
    burr_link(compliance_linear_predictor)
  } else {
    pnorm(compliance_linear_predictor)
  }
  outcome_data$D <- rbinom(length(d_probs), 1, d_probs)

  x_main_treat <- cbind(1, as.matrix(x_matrix))
  x_main_treat <- cbind(x_main_treat, 0 * outcome_data$D)
  x_main_treat <- cbind(x_main_treat, outcome_data$D * cbind(as.matrix(x_matrix), 1))

  x_main_control <- cbind(1, as.matrix(x_matrix))
  x_main_control <- cbind(x_main_control, outcome_data$D)
  x_main_control <- cbind(x_main_control, 0 * outcome_data$D * cbind(as.matrix(x_matrix), 1))

  coeff_beta_d0 <- matrix(NA_real_, nrow = ncol(x_matrix), ncol = K)
  coeff_beta_d1 <- matrix(NA_real_, nrow = ncol(x_matrix), ncol = K)
  coeff_d0 <- rep(NA_real_, K)
  coeff_d1 <- rep(NA_real_, K)
  coeff_mu_y <- rep(NA_real_, K)
  phi_y_fac <- rep(NA_real_, J)
  tau2_phi_y <- 9

  for (k in seq_len(K)) {
    coeff_beta_d0[, k] <- rep(k, ncol(x_matrix))
    coeff_beta_d1[, k] <- rep(k, ncol(x_matrix))
    coeff_d0[k] <- k
    coeff_d1[k] <- k + 4.5 + (k - 1)
    coeff_mu_y[k] <- 2 * k
    phi_y_fac[facility_strata == k] <- rnorm(sum(facility_strata == k), 0, sqrt(tau2_phi_y))
  }

  coeff_full_y <- rbind(coeff_mu_y, coeff_beta_d0, coeff_d0, coeff_beta_d1, coeff_d1)
  phi_y_person <- phi_y_fac[outcome_data$Facility]
  outcome_data$PhiY <- phi_y_person
  sigma2_y <- c(16, 16)

  y0_means <- diag(x_main_control %*% coeff_full_y[, outcome_data$Strata, drop = FALSE]) +
    phi_y_person +
    diag(x_main_control[, 2] * x_main_control[, 3]) %*%
      config$interaction_coeff_control[outcome_data$Strata]

  y1_means <- diag(x_main_treat %*% coeff_full_y[, outcome_data$Strata, drop = FALSE]) +
    phi_y_person +
    diag(x_main_control[, 2] * x_main_control[, 3]) %*%
      config$interaction_coeff_control[outcome_data$Strata] +
    outcome_data$D *
      diag(x_main_treat[, 2] * x_main_treat[, 3]) %*%
      config$interaction_coeff_treat[outcome_data$Strata]

  outcome_data$Y0 <- rnorm(N, y0_means, sqrt(sigma2_y[outcome_data$Strata]))
  outcome_data$Y1 <- rnorm(N, y1_means, sqrt(sigma2_y[outcome_data$Strata]))

  list(
    X = x_matrix,
    OutcomeData = outcome_data,
    Facility_Strata = facility_strata,
    True_CZ = true_cz,
    True_C = as.matrix(true_cz[, seq_len(L), drop = FALSE]),
    True_Z = as.matrix(true_cz[, L + seq_len(P), drop = FALSE]),
    Coeff_Alpha = coeff_alpha,
    PhiD_Fac = phi_d_fac,
    Tau2_PhiD = tau2_phi_d,
    Coeff_BetaD0 = coeff_beta_d0,
    Coeff_BetaD1 = coeff_beta_d1,
    Coeff_D0 = coeff_d0,
    Coeff_D1 = coeff_d1,
    Coeff_muY = coeff_mu_y,
    PhiY_Fac = phi_y_fac,
    Tau2_PhiY = tau2_phi_y,
    Sigma2Y = sigma2_y,
    Sigma = sigma
  )
}

# The population contrasts weight eligible residents, not datasets.
# A replication with no eligible residents contributes zero to both totals.
summarize_truth_totals <- function(outcome_data) {
  with(outcome_data, {
    contrast <- Y1 - Y0
    groups <- list(ITT = rep(TRUE, length(contrast)),
                   ITT1 = Strata == 1, ITT2 = Strata == 2,
                   CACE = D == 1,
                   CACE1 = D == 1 & Strata == 1, CACE2 = D == 1 & Strata == 2)
    vapply(groups, function(in_group) {
      c(contrast_sum = sum(contrast[in_group]), count = sum(in_group))
    }, numeric(2))
  })
}

# Construct the fixed parameter values implied by the case-study specification,
# excluding realized random effects and simulated outcomes.
make_baseline_case_study_parameters <- function(x_population, config) {
  n_covariates <- ncol(as.matrix(x_population))
  x_design_cols <- n_covariates + 1L

  coeff_alpha <- matrix(NA_real_, nrow = x_design_cols, ncol = config$K)
  coeff_beta_d0 <- matrix(NA_real_, nrow = n_covariates, ncol = config$K)
  coeff_beta_d1 <- matrix(NA_real_, nrow = n_covariates, ncol = config$K)
  coeff_d0 <- rep(NA_real_, config$K)
  coeff_d1 <- rep(NA_real_, config$K)
  coeff_mu_y <- rep(NA_real_, config$K)

  for (k in seq_len(config$K)) {
    coeff_alpha[, k] <- c((k - 1) / 2, rep(-0.25 * k, x_design_cols - 1L))
    coeff_beta_d0[, k] <- rep(k, n_covariates)
    coeff_beta_d1[, k] <- rep(k, n_covariates)
    coeff_d0[k] <- k
    coeff_d1[k] <- k + 4.5 + (k - 1)
    coeff_mu_y[k] <- 2 * k
  }

  list(
    coeff_alpha = coeff_alpha,
    coeff_beta_d0 = coeff_beta_d0,
    coeff_beta_d1 = coeff_beta_d1,
    coeff_d0 = coeff_d0,
    coeff_d1 = coeff_d1,
    coeff_mu_y = coeff_mu_y,
    tau2_phi_d = 0.25
  )
}

# Indicate whether the active case study satisfies the assumptions under which
# the appendix closed-form estimand formulas are available.
supports_analytic_truth_check <- function(config) {
  !isTRUE(config$use_burr) &&
    all(config$interaction_coeff_treat == 0) &&
    all(config$interaction_coeff_control == 0)
}

# Evaluate the closed-form estimands given the realized data-generating
# coefficients and the empirical covariate population. This provides a direct
# comparison against the Monte Carlo approximation to the truth.
compute_closed_form_truth_from_dgm <- function(x_population, config, dgm_sample) {
  if (!supports_analytic_truth_check(config)) {
    stop(
      "Analytic truth formulas are only available for the probit, no-interaction case studies.",
      call. = FALSE
    )
  }

  x_population <- as.matrix(x_population)
  x_alpha <- cbind(1, x_population) %*% dgm_sample$Coeff_Alpha
  p_ijk <- pnorm(x_alpha / sqrt(dgm_sample$Tau2_PhiD + 1))

  mu_k_plus_x_beta0 <- sweep(x_population %*% dgm_sample$Coeff_BetaD0, 2, dgm_sample$Coeff_muY, "+")
  delta1_k_plus_x_beta1 <- sweep(x_population %*% dgm_sample$Coeff_BetaD1, 2, dgm_sample$Coeff_D1, "+")
  p_ijk_times_delta1_k_plus_x_beta1 <- p_ijk * delta1_k_plus_x_beta1

  itt_y1_components <- mu_k_plus_x_beta0 + p_ijk_times_delta1_k_plus_x_beta1
  true_itt_y1 <- mean(rowSums(sweep(itt_y1_components, 2, config$pi_k, "*")))

  p_ijk_times_delta0_k <- sweep(p_ijk, 2, dgm_sample$Coeff_D0, "*")
  itt_y0_components <- mu_k_plus_x_beta0 + p_ijk_times_delta0_k
  true_itt_y0 <- mean(rowSums(sweep(itt_y0_components, 2, config$pi_k, "*")))

  p_ijk_mean <- colMeans(p_ijk)
  den_cace <- sum(config$pi_k * p_ijk_mean)

  p_ijk_weighted_x_beta_both <- colMeans(p_ijk * (x_population %*% (dgm_sample$Coeff_BetaD0 + dgm_sample$Coeff_BetaD1))) /
    p_ijk_mean
  mu_k_delta1_plus_weighted_x_beta <- dgm_sample$Coeff_muY + p_ijk_weighted_x_beta_both + dgm_sample$Coeff_D1
  num_cace_y1 <- sum(config$pi_k * colMeans(sweep(p_ijk, 2, mu_k_delta1_plus_weighted_x_beta, "*")))
  true_cace_y1 <- num_cace_y1 / den_cace

  p_ijk_weighted_x_beta0 <- colMeans(p_ijk * (x_population %*% dgm_sample$Coeff_BetaD0)) /
    p_ijk_mean
  mu_k_delta0_plus_weighted_x_beta0 <- dgm_sample$Coeff_muY + p_ijk_weighted_x_beta0 + dgm_sample$Coeff_D0
  num_cace_y0 <- sum(config$pi_k * colMeans(sweep(p_ijk, 2, mu_k_delta0_plus_weighted_x_beta0, "*")))
  true_cace_y0 <- num_cace_y0 / den_cace

  itt_y1_by_stratum <- rep(NA_real_, config$K)
  itt_y0_by_stratum <- rep(NA_real_, config$K)
  itt_by_stratum <- rep(NA_real_, config$K)
  cace_y1_by_stratum <- rep(NA_real_, config$K)
  cace_y0_by_stratum <- rep(NA_real_, config$K)
  cace_by_stratum <- rep(NA_real_, config$K)

  for (k in seq_len(config$K)) {
    p_k <- pnorm((cbind(1, x_population) %*% dgm_sample$Coeff_Alpha[, k]) / sqrt(dgm_sample$Tau2_PhiD + 1))
    mean_p_k <- mean(p_k)

    mu_y_plus_x_beta0_k <- dgm_sample$Coeff_muY[k] + mean(x_population %*% dgm_sample$Coeff_BetaD0[, k])
    weighted_x_beta1_delta1_k <- mean(p_k * (x_population %*% dgm_sample$Coeff_BetaD1[, k] + dgm_sample$Coeff_D1[k]))
    weighted_delta0_k <- mean(p_k * dgm_sample$Coeff_D0[k])

    itt_y1_by_stratum[k] <- mu_y_plus_x_beta0_k + weighted_x_beta1_delta1_k
    itt_y0_by_stratum[k] <- mu_y_plus_x_beta0_k + weighted_delta0_k
    itt_by_stratum[k] <- itt_y1_by_stratum[k] - itt_y0_by_stratum[k]

    weighted_x_beta0_beta1_k <- mean(p_k * x_population %*% (dgm_sample$Coeff_BetaD0[, k] + dgm_sample$Coeff_BetaD1[, k]))
    weighted_x_beta0_k <- mean(p_k * x_population %*% dgm_sample$Coeff_BetaD0[, k])
    cace_y1_by_stratum[k] <- dgm_sample$Coeff_muY[k] + weighted_x_beta0_beta1_k / mean_p_k + dgm_sample$Coeff_D1[k]
    cace_y0_by_stratum[k] <- dgm_sample$Coeff_muY[k] + weighted_x_beta0_k / mean_p_k + dgm_sample$Coeff_D0[k]
    cace_by_stratum[k] <- cace_y1_by_stratum[k] - cace_y0_by_stratum[k]
  }

  cace_num <- sum(config$pi_k * cace_by_stratum * p_ijk_mean)
  cace_den <- sum(config$pi_k * p_ijk_mean)

  c(
    ITT = sum(config$pi_k * itt_by_stratum),
    ITT1 = itt_by_stratum[1],
    ITT2 = itt_by_stratum[2],
    CACE = cace_num / cace_den,
    CACE1 = cace_by_stratum[1],
    CACE2 = cace_by_stratum[2],
    ITT_Y1 = true_itt_y1,
    ITT_Y0 = true_itt_y0,
    CACE_Y1 = true_cace_y1,
    CACE_Y0 = true_cace_y0,
    ITT1_Y1 = itt_y1_by_stratum[1],
    ITT1_Y0 = itt_y0_by_stratum[1],
    ITT2_Y1 = itt_y1_by_stratum[2],
    ITT2_Y0 = itt_y0_by_stratum[2],
    CACE1_Y1 = cace_y1_by_stratum[1],
    CACE1_Y0 = cace_y0_by_stratum[1],
    CACE2_Y1 = cace_y1_by_stratum[2],
    CACE2_Y0 = cace_y0_by_stratum[2]
  )
}

# Evaluate the same appendix closed-form estimands using the fixed baseline
# Case Study parameter values rather than parameter values carried inside one
# realized simulated dataset object.
compute_closed_form_truth_from_baseline_parameters <- function(x_population, config) {
  baseline_params <- make_baseline_case_study_parameters(
    x_population = x_population,
    config = config
  )

  x_population <- as.matrix(x_population)
  x_alpha <- cbind(1, x_population) %*% baseline_params$coeff_alpha
  p_ijk <- pnorm(x_alpha / sqrt(baseline_params$tau2_phi_d + 1))

  mu_k_plus_x_beta0 <- sweep(x_population %*% baseline_params$coeff_beta_d0, 2, baseline_params$coeff_mu_y, "+")
  delta1_k_plus_x_beta1 <- sweep(x_population %*% baseline_params$coeff_beta_d1, 2, baseline_params$coeff_d1, "+")
  p_ijk_times_delta1_k_plus_x_beta1 <- p_ijk * delta1_k_plus_x_beta1

  itt_y1_components <- mu_k_plus_x_beta0 + p_ijk_times_delta1_k_plus_x_beta1
  true_itt_y1 <- mean(rowSums(sweep(itt_y1_components, 2, config$pi_k, "*")))

  p_ijk_times_delta0_k <- sweep(p_ijk, 2, baseline_params$coeff_d0, "*")
  itt_y0_components <- mu_k_plus_x_beta0 + p_ijk_times_delta0_k
  true_itt_y0 <- mean(rowSums(sweep(itt_y0_components, 2, config$pi_k, "*")))

  p_ijk_mean <- colMeans(p_ijk)
  den_cace <- sum(config$pi_k * p_ijk_mean)

  p_ijk_weighted_x_beta_both <- colMeans(p_ijk * (x_population %*% (baseline_params$coeff_beta_d0 + baseline_params$coeff_beta_d1))) /
    p_ijk_mean
  mu_k_delta1_plus_weighted_x_beta <- baseline_params$coeff_mu_y + p_ijk_weighted_x_beta_both + baseline_params$coeff_d1
  num_cace_y1 <- sum(config$pi_k * colMeans(sweep(p_ijk, 2, mu_k_delta1_plus_weighted_x_beta, "*")))
  true_cace_y1 <- num_cace_y1 / den_cace

  p_ijk_weighted_x_beta0 <- colMeans(p_ijk * (x_population %*% baseline_params$coeff_beta_d0)) /
    p_ijk_mean
  mu_k_delta0_plus_weighted_x_beta0 <- baseline_params$coeff_mu_y + p_ijk_weighted_x_beta0 + baseline_params$coeff_d0
  num_cace_y0 <- sum(config$pi_k * colMeans(sweep(p_ijk, 2, mu_k_delta0_plus_weighted_x_beta0, "*")))
  true_cace_y0 <- num_cace_y0 / den_cace

  itt_y1_by_stratum <- rep(NA_real_, config$K)
  itt_y0_by_stratum <- rep(NA_real_, config$K)
  itt_by_stratum <- rep(NA_real_, config$K)
  cace_y1_by_stratum <- rep(NA_real_, config$K)
  cace_y0_by_stratum <- rep(NA_real_, config$K)
  cace_by_stratum <- rep(NA_real_, config$K)

  for (k in seq_len(config$K)) {
    p_k <- pnorm((cbind(1, x_population) %*% baseline_params$coeff_alpha[, k]) / sqrt(baseline_params$tau2_phi_d + 1))
    mean_p_k <- mean(p_k)

    mu_y_plus_x_beta0_k <- baseline_params$coeff_mu_y[k] + mean(x_population %*% baseline_params$coeff_beta_d0[, k])
    weighted_x_beta1_delta1_k <- mean(p_k * (x_population %*% baseline_params$coeff_beta_d1[, k] + baseline_params$coeff_d1[k]))
    weighted_delta0_k <- mean(p_k * baseline_params$coeff_d0[k])

    itt_y1_by_stratum[k] <- mu_y_plus_x_beta0_k + weighted_x_beta1_delta1_k
    itt_y0_by_stratum[k] <- mu_y_plus_x_beta0_k + weighted_delta0_k
    itt_by_stratum[k] <- itt_y1_by_stratum[k] - itt_y0_by_stratum[k]

    weighted_x_beta0_beta1_k <- mean(p_k * x_population %*% (baseline_params$coeff_beta_d0[, k] + baseline_params$coeff_beta_d1[, k]))
    weighted_x_beta0_k <- mean(p_k * x_population %*% baseline_params$coeff_beta_d0[, k])
    cace_y1_by_stratum[k] <- baseline_params$coeff_mu_y[k] + weighted_x_beta0_beta1_k / mean_p_k + baseline_params$coeff_d1[k]
    cace_y0_by_stratum[k] <- baseline_params$coeff_mu_y[k] + weighted_x_beta0_k / mean_p_k + baseline_params$coeff_d0[k]
    cace_by_stratum[k] <- cace_y1_by_stratum[k] - cace_y0_by_stratum[k]
  }

  cace_num <- sum(config$pi_k * cace_by_stratum * p_ijk_mean)
  cace_den <- sum(config$pi_k * p_ijk_mean)

  c(
    ITT = sum(config$pi_k * itt_by_stratum),
    ITT1 = itt_by_stratum[1],
    ITT2 = itt_by_stratum[2],
    CACE = cace_num / cace_den,
    CACE1 = cace_by_stratum[1],
    CACE2 = cace_by_stratum[2],
    ITT_Y1 = true_itt_y1,
    ITT_Y0 = true_itt_y0,
    CACE_Y1 = true_cace_y1,
    CACE_Y0 = true_cace_y0,
    ITT1_Y1 = itt_y1_by_stratum[1],
    ITT1_Y0 = itt_y0_by_stratum[1],
    ITT2_Y1 = itt_y1_by_stratum[2],
    ITT2_Y0 = itt_y0_by_stratum[2],
    CACE1_Y1 = cace_y1_by_stratum[1],
    CACE1_Y0 = cace_y0_by_stratum[1],
    CACE2_Y1 = cace_y1_by_stratum[2],
    CACE2_Y0 = cace_y0_by_stratum[2]
  )
}

# Compare the Monte Carlo truth estimates from repeated simulated datasets to
# the corresponding analytic formulas evaluated at the true DGM parameters and
# at the fixed baseline Case Study parameter values. This produces a compact
# table that can be printed, saved, or checked in regression tests.
compare_monte_carlo_and_analytic_truth <- function(x_population, config, dgm_sample, monte_carlo_truth) {
  closed_form_truth <- compute_closed_form_truth_from_dgm(
    x_population = x_population,
    config = config,
    dgm_sample = dgm_sample
  )
  baseline_closed_form_truth <- compute_closed_form_truth_from_baseline_parameters(
    x_population = x_population,
    config = config
  )
  estimand_names <- names(monte_carlo_truth)
  closed_form_truth <- closed_form_truth[estimand_names]
  baseline_closed_form_truth <- baseline_closed_form_truth[estimand_names]

  data.frame(
    estimand = estimand_names,
    monte_carlo = unname(monte_carlo_truth[estimand_names]),
    analytic = unname(closed_form_truth),
    baseline_analytic = unname(baseline_closed_form_truth),
    abs_diff = abs(unname(monte_carlo_truth[estimand_names]) - unname(closed_form_truth)),
    baseline_abs_diff = abs(unname(monte_carlo_truth[estimand_names]) - unname(baseline_closed_form_truth)),
    row.names = NULL
  )
}

# Build the probit compliance-model design matrix used throughout the Gibbs
# sampler. The first column is an intercept and the remaining columns are the
# resident-level covariates from the simulation covariate file.
build_compliance_design_matrix <- function(x_matrix) {
  x_d <- cbind(Intercept = 1, as.matrix(x_matrix))
  colnames(x_d) <- paste0("V", seq_len(ncol(x_d)))
  x_d
}

# Build the reduced outcome-model design matrix used in the simulation study.
# This matches the appendix notation by including baseline terms, the control
# complier indicator, and the treatment-by-compliance terms among treated units.
build_outcome_design_matrix <- function(x_matrix, treat, d_value) {
  x_matrix <- as.matrix(x_matrix)
  d_value <- as.numeric(d_value)
  treat <- as.numeric(treat)

  design <- cbind(1, x_matrix)
  design <- cbind(design, (1 - treat) * d_value)
  design <- cbind(design, treat * d_value * cbind(x_matrix, 1))
  colnames(design) <- paste0("V", seq_len(ncol(design)))
  design
}

# Build the outcome design matrices repeatedly needed by the Gibbs sampler:
# the observed-data matrix, the all-treated matrix, the current all-control
# matrix, and the control matrices under complier/noncomplier assumptions.
build_outcome_design_matrices <- function(x_matrix, treat, d_value) {
  list(
    observed = build_outcome_design_matrix(
      x_matrix = x_matrix,
      treat = treat,
      d_value = d_value
    ),
    treat_all = build_outcome_design_matrix(
      x_matrix = x_matrix,
      treat = rep(1, length(treat)),
      d_value = d_value
    ),
    control_current = build_outcome_design_matrix(
      x_matrix = x_matrix,
      treat = rep(0, length(treat)),
      d_value = d_value
    ),
    control_comp = build_outcome_design_matrix(
      x_matrix = x_matrix,
      treat = rep(0, length(treat)),
      d_value = rep(1, length(treat))
    ),
    control_nocomp = build_outcome_design_matrix(
      x_matrix = x_matrix,
      treat = rep(0, length(treat)),
      d_value = rep(0, length(treat))
    )
  )
}

# Combine the observed potential outcomes into the observed outcome vector used
# in the outcome-model updates.
build_observed_outcomes <- function(y0, y1, treat) {
  y_obs <- rep(NA_real_, length(treat))
  y_obs[treat == 1] <- y1[treat == 1]
  y_obs[treat == 0] <- y0[treat == 0]
  y_obs
}

# Sample the latent implementation-strata mean vectors from the conjugate MVN
# full conditional implied by the Normal mixture model.
update_mixture_component_means <- function(cz_i, strata_fac, sigma_current, V_0_inv, m_0, K) {
  sigma_inv <- solve(sigma_current)
  num_features <- ncol(cz_i)
  mu_draws <- matrix(NA_real_, nrow = num_features, ncol = K)

  for (k in seq_len(K)) {
    fac_in_k <- strata_fac == k
    J_k <- sum(fac_in_k)
    posterior_cov <- solve(V_0_inv + J_k * sigma_inv)

    if (J_k > 0) {
      cz_k_bar <- matrix(colMeans(cz_i[fac_in_k, , drop = FALSE]), ncol = 1)
      posterior_mean <- posterior_cov %*% (V_0_inv %*% m_0 + J_k * sigma_inv %*% cz_k_bar)
    } else {
      posterior_mean <- posterior_cov %*% (V_0_inv %*% m_0)
    }

    mu_draws[, k] <- MASS::mvrnorm(1, posterior_mean, posterior_cov)
  }

  mu_draws
}

# Sample the Albert-Chib latent Gaussian variables for the probit compliance
# model conditional on the current latent compliance indicators.
update_latent_probit_variables <- function(D_current, strata_person, X_D, alpha_current, phiD_person, K) {
  U_new <- rep(NA_real_, length(D_current))

  for (k in seq_len(K)) {
    people_in_k <- strata_person == k
    if (!any(people_in_k)) {
      next
    }

    D1_idx <- which(D_current == 1 & people_in_k)
    D0_idx <- which(D_current == 0 & people_in_k)

    if (length(D1_idx) > 0) {
      U_new[D1_idx] <- truncnorm::rtruncnorm(
        length(D1_idx),
        a = 0,
        b = Inf,
        mean = X_D[D1_idx, , drop = FALSE] %*% alpha_current[[k]] + phiD_person[D1_idx]
      )
    }

    if (length(D0_idx) > 0) {
      U_new[D0_idx] <- truncnorm::rtruncnorm(
        length(D0_idx),
        a = -Inf,
        b = 0,
        mean = X_D[D0_idx, , drop = FALSE] %*% alpha_current[[k]] + phiD_person[D0_idx]
      )
    }
  }

  U_new
}

# Sample a shared random-effect variance under a Uniform(0, U) prior on the
# standard deviation tau, which implies a truncated inverse-gamma full
# conditional for tau^2 with shape (J - 1) / 2 in the fixed-J simulation study.
update_shared_random_effect_variance <- function(phi_draw, upper_bound_tau2, context = NULL) {
  fail <- function(reason) {
    where <- if (is.null(context)) "" else paste0(" (", context, ")")
    stop("Random-effect variance update", where, ": ", reason, call. = FALSE)
  }

  shape <- (length(phi_draw) - 1) / 2
  rate <- 1 / 2 * sum(phi_draw^2)
  if (shape <= 0 || !is.finite(rate) || rate <= 0 ||
      length(upper_bound_tau2) != 1L || !is.finite(upper_bound_tau2) ||
      upper_bound_tau2 <= 0) {
    fail("shape, rate, and variance bound must be finite and positive")
  }

  # Under a uniform prior on the SD, the precision conditional is Gamma
  # restricted to (1 / upper_bound_tau2, Inf). Invert its log survival function.
  log_tail <- stats::pgamma(1 / upper_bound_tau2, shape, rate = rate,
                            lower.tail = FALSE, log.p = TRUE)
  if (!is.finite(log_tail)) fail("the log survival probability at the bound is not finite")
  log_prob <- log_tail + log(stats::runif(1L))
  precision <- stats::qgamma(log_prob, shape, rate = rate, lower.tail = FALSE, log.p = TRUE)
  variance <- 1 / precision
  if (!is.finite(variance) || variance <= 0 || variance >= upper_bound_tau2) {
    fail("the sampled variance is not strictly inside the truncation interval")
  }
  variance
}

# Sample the outcome residual variances from their conjugate inverse-gamma full
# conditional.
update_outcome_sigmas <- function(Y_obs, X_Main_obs, beta_new, phiY_person, strata_person, a_sigmas, b_sigmas, K) {
  sigma2_draws <- rep(NA_real_, K)

  for (k in seq_len(K)) {
    people_in_k <- strata_person == k
    sigma2_shape <- a_sigmas + sum(people_in_k) / 2

    if (sum(people_in_k) > 0) {
      X_in_k <- X_Main_obs[people_in_k, , drop = FALSE]
      sigma2_ss <- sum(
        (
          Y_obs[people_in_k] -
            X_in_k %*% beta_new[[k]] -
            phiY_person[people_in_k]
        )^2
      )
      sigma2_rate <- b_sigmas + 1 / 2 * sigma2_ss
    } else {
      sigma2_rate <- b_sigmas
    }

    sigma2_draws[k] <- invgamma::rinvgamma(1, sigma2_shape, sigma2_rate)
  }

  sigma2_draws
}

# Impute the latent compliance indicators for individuals in control facilities
# using the outcome and probit compliance models.
impute_unobserved_compliance <- function(D_current, Y_obs, X_Main_control_comp, X_Main_control_nocomp, X_D, alpha_new, beta_new, sigma2_new, phiY_person, phiD_person, strata_person, treat_person_bool, K) {
  D_new <- D_current
  D_samp_probs <- rep(NA_real_, length(D_current))

  for (k in seq_len(K)) {
    UnObs_D_k <- strata_person == k & !treat_person_bool

    if (sum(UnObs_D_k) > 0) {
      Y_UnObs_D_k <- Y_obs[UnObs_D_k]
      Y_mean_comp <- X_Main_control_comp[UnObs_D_k, , drop = FALSE] %*% beta_new[[k]] + phiY_person[UnObs_D_k]
      Y_mean_nocomp <- X_Main_control_nocomp[UnObs_D_k, , drop = FALSE] %*% beta_new[[k]] + phiY_person[UnObs_D_k]

      comp_dens <- dnorm(Y_UnObs_D_k, Y_mean_comp, sqrt(sigma2_new[k])) *
        pnorm(X_D[UnObs_D_k, , drop = FALSE] %*% alpha_new[[k]] + phiD_person[UnObs_D_k])

      noncomp_dens <- dnorm(Y_UnObs_D_k, Y_mean_nocomp, sqrt(sigma2_new[k])) *
        (1 - pnorm(X_D[UnObs_D_k, , drop = FALSE] %*% alpha_new[[k]] + phiD_person[UnObs_D_k]))

      D_samp_probs[UnObs_D_k] <- comp_dens / (comp_dens + noncomp_dens)
      D_new[UnObs_D_k] <- rbinom(sum(UnObs_D_k), 1, D_samp_probs[UnObs_D_k])
    }
  }

  list(D_new = D_new, probabilities = D_samp_probs)
}

# Identify the ordered implementation-stratum labels used for reporting the
# stratum-specific estimands. The simulation study treats the stratum with the
# smaller first mixture-mean coordinate as "type 1" to stabilize labels across
# posterior draws.
identify_ordered_strata <- function(mu_draws) {
  ordered_index <- order(mu_draws[1, ])
  list(lower = ordered_index[1], upper = ordered_index[length(ordered_index)])
}

# Summarize operating-characteristic metrics for one estimand across the
# retained posterior draws in a single simulation replication.
summarize_estimand_performance <- function(draw_values, true_value) {
  lower <- unname(quantile(draw_values, 0.025, na.rm = TRUE))
  upper <- unname(quantile(draw_values, 0.975, na.rm = TRUE))
  posterior_mean <- mean(draw_values, na.rm = TRUE)
  posterior_sd <- stats::sd(draw_values, na.rm = TRUE)

  c(
    post_mean = posterior_mean,
    post_sd = posterior_sd,
    bias = posterior_mean - true_value,
    cover = as.numeric(true_value >= lower & true_value <= upper),
    iw = upper - lower,
    lower = lower,
    upper = upper
  )
}

# Assemble the true super-population estimands for the current case study into
# a compact table that can be printed, inspected, or written to disk by the
# main script.
build_true_estimand_table <- function(case_study, case_label, true_itt, true_itt1, true_itt2, true_cace, true_cace1, true_cace2) {
  data.frame(
    case_study = rep(case_study, 6),
    case_label = rep(case_label, 6),
    estimand = c("ITT", "ITT1", "ITT2", "CACE", "CACE1", "CACE2"),
    true_value = c(true_itt, true_itt1, true_itt2, true_cace, true_cace1, true_cace2),
    row.names = NULL
  )
}

# Compare posterior uncertainty within simulated datasets with the empirical
# variability of posterior point estimates across simulated datasets. The median
# posterior-SD ratio is reported on the estimand scale so that the saved CSVs
# match the manuscript tables: median posterior SD across simulation
# replications divided by the empirical SD of posterior means across
# replications. Its square is the corresponding median variance-ratio
# diagnostic requested by the reviewers.
summarize_uncertainty_ratio <- function(post_mean, post_sd) {
  median_post_sd <- stats::median(post_sd, na.rm = TRUE)
  empirical_sd_post_mean <- stats::sd(post_mean, na.rm = TRUE)
  sd_ratio <- median_post_sd / empirical_sd_post_mean

  c(
    median_post_sd = median_post_sd,
    empirical_sd_post_mean = empirical_sd_post_mean,
    sd_ratio = sd_ratio,
    variance_ratio = sd_ratio^2
  )
}
