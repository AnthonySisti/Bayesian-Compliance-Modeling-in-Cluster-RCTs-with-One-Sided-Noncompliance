# Model helper functions for the binary outcome analysis.
#
# These functions implement the conditional updates used in the Gibbs sampler
# for the latent principal-stratification model. The compliance and outcome
# models are both probit models and are sampled using latent Gaussian
# augmentation. Missing cluster-level implementation measures, latent individual
# compliance indicators, missing potential outcomes, and cluster-specific random
# effects are updated as part of the posterior computation.
#
# Super-population estimands are evaluated from integrated probit risks at
# each retained posterior draw. Predictive datasets are still simulated for
# the posterior predictive checks.

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

load_binary_outcome_packages <- function() {
  library("boot", include.only = "inv.logit")
  library("mixtools")
  library(coda)
  library(MCMCpack, include.only = "riwish")
  library(DirichletReg, include.only = "rdirichlet")
  library(mixtools)
  library(mvtnorm, include.only = c("dmvnorm", "rmvnorm"))
  library(MASS, include.only = "mvrnorm")
  library(scales, include.only = "alpha")
  library(invgamma)
  library(dplyr)
  invisible(TRUE)
}

make_binary_outcome_config <- function(
    seed = 20260422L,
    burn = 500L,
    thin = 3L,
    n_imps = 1000L,
    K = 2L,
    plot_initial_clusters = TRUE,
    output_dir = NULL) {
  if (is.null(output_dir)) {
    output_dir <- file.path(get_script_dir(), "outputs")
  }

  burn <- as.integer(burn)
  thin <- as.integer(thin)
  n_imps <- as.integer(n_imps)

  list(
    seed = as.integer(seed),
    burn = burn,
    thin = thin,
    n_imps = n_imps,
    K = as.integer(K),
    Sims = burn + thin * n_imps,
    samps_use = burn + seq_len(n_imps) * thin,
    plot_initial_clusters = isTRUE(plot_initial_clusters),
    output_dir = output_dir
  )
}

update_mixture_means_original <- function(CZ_i, strata_fac, Sigma_prev, V_0_inv, m_0, K) {
  Sigma_inv <- solve(Sigma_prev)
  L_plus_P <- ncol(CZ_i)
  new_mu <- matrix(NA_real_, nrow = L_plus_P, ncol = K)

  for (k in seq_len(K)) {
    fac_in_k <- strata_fac == k
    J_k <- sum(fac_in_k)
    VJSig_inv <- solve(V_0_inv + J_k * Sigma_inv)

    if (J_k > 0) {
      CZ_k_bar <- as.matrix(colMeans(CZ_i[fac_in_k, , drop = FALSE]))
      mu_hat_k <- VJSig_inv %*% (V_0_inv %*% m_0 + J_k * Sigma_inv %*% CZ_k_bar)
    } else {
      mu_hat_k <- VJSig_inv %*% (V_0_inv %*% m_0)
    }

    new_mu[, k] <- MASS::mvrnorm(1, mu_hat_k, VJSig_inv)
  }

  new_mu
}

update_mixture_covariance_original <- function(C_i, Z_i, strata_fac, mu_draws, S_0, nu_0, J, K) {
  S_J <- S_0

  for (k in seq_len(K)) {
    fac_in_k <- strata_fac == k
    J_k <- sum(fac_in_k)
    if (J_k > 0) {
      CZ_k <- rbind(t(C_i[fac_in_k, , drop = FALSE]), t(Z_i[fac_in_k, , drop = FALSE]))
      CZ_k_minus_mu_k <- sweep(CZ_k, 1, mu_draws[, k])
      S_J <- S_J + Reduce(
        "+",
        lapply(seq_len(J_k), function(j) CZ_k_minus_mu_k[, j] %*% t(CZ_k_minus_mu_k[, j]))
      )
    }
  }

  MCMCpack::riwish(nu_0 + J, S_J)
}

update_mixture_weights_original <- function(strata_fac, lambda_pk, K) {
  num_fac_each_k <- numeric(0)
  for (k in seq_len(K)) {
    num_fac_each_k <- c(num_fac_each_k, sum(strata_fac == k))
  }
  DirichletReg::rdirichlet(1, lambda_pk + num_fac_each_k)
}

sample_latent_u_original <- function(D_current, X_D, alpha_list_prev, phiD_person, strata_person, K) {
  U_new <- rep(NA_real_, length(D_current))

  for (k in seq_len(K)) {
    if (sum(strata_person == k) > 0) {
      D1_idx <- which(D_current == 1 & strata_person == k)
      D0_idx <- which(D_current == 0 & strata_person == k)

      if (length(D1_idx) > 0) {
        U_new[D1_idx] <- truncnorm::rtruncnorm(
          length(D1_idx),
          a = 0,
          b = Inf,
          mean = X_D[D1_idx, , drop = FALSE] %*% alpha_list_prev[[k]] + phiD_person[D1_idx]
        )
      }
      if (length(D0_idx) > 0) {
        U_new[D0_idx] <- truncnorm::rtruncnorm(
          length(D0_idx),
          a = -Inf,
          b = 0,
          mean = X_D[D0_idx, , drop = FALSE] %*% alpha_list_prev[[k]] + phiD_person[D0_idx]
        )
      }
    }
  }

  U_new
}

update_compliance_random_effects_original <- function(U_new, X_D, alpha_list_prev, strata_person, strata_fac, provider_fac, person_fac, tau2_prev, Nj, K) {
  phi_draw <- rep(0, length(provider_fac))

  for (k in seq_len(K)) {
    people_in_k <- strata_person == k
    fac_in_k <- provider_fac[strata_fac == k]
    num_fac_in_k <- length(fac_in_k)
    peop_fac <- person_fac[people_in_k]

    if (num_fac_in_k > 0) {
      U_phi <- U_new[people_in_k] - X_D[people_in_k, , drop = FALSE] %*% alpha_list_prev[[k]]
      Njk <- Nj[fac_in_k]
      var_phi <- tau2_prev / (Njk * tau2_prev + 1)
      mean_den_phi <- Njk + (1 / tau2_prev)
      mean_num_phi <- tapply(U_phi, peop_fac, sum)
      mu_phi <- matrix(mean_num_phi / mean_den_phi, nrow = num_fac_in_k, ncol = 1)
      phi_draw[fac_in_k] <- MASS::mvrnorm(1, mu = mu_phi, Sigma = diag(var_phi, num_fac_in_k))
    }
  }

  phi_draw
}

update_truncated_variance_original <- function(phi_draw, upper_bound_tau2, context = NULL) {
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

update_compliance_coefficients_step2 <- function(U_new, X_D, strata_person, phiD_person, v_alpha_inv, m_alpha, K) {
  X_shared <- X_D[, -1, drop = FALSE]
  strata_indicators <- sapply(seq_len(K), function(k) as.numeric(strata_person == k))
  if (!is.matrix(strata_indicators)) {
    strata_indicators <- matrix(strata_indicators, ncol = K)
  }
  design_matrix <- cbind(strata_indicators, X_shared)
  design_matrix_t <- t(design_matrix)
  posterior_cov <- chol2inv(chol(design_matrix_t %*% design_matrix + v_alpha_inv))
  posterior_mean <- posterior_cov %*% (design_matrix_t %*% (U_new - phiD_person) + v_alpha_inv %*% m_alpha)
  coeff_draw <- MASS::mvrnorm(1, mu = posterior_mean, Sigma = posterior_cov)

  mu_draw <- coeff_draw[seq_len(K)]
  alpha_shared_draw <- coeff_draw[K + seq_len(ncol(X_shared))]
  alpha_list_draw <- lapply(seq_len(K), function(k) c(mu_draw[k], alpha_shared_draw))

  list(
    mu_draw = mu_draw,
    alpha_shared_draw = alpha_shared_draw,
    alpha_list_draw = alpha_list_draw,
    coeff_draw = coeff_draw
  )
}

sample_latent_q_original <- function(Y_current, X_Main_obs, beta_full_prev, phiY_person, strata_person, K) {
  Q_new <- rep(NA_real_, length(Y_current))

  for (k in seq_len(K)) {
    if (sum(strata_person == k) > 0) {
      Y1_idx <- which(Y_current == 1 & strata_person == k)
      Y0_idx <- which(Y_current == 0 & strata_person == k)

      if (length(Y1_idx) > 0) {
        Q_new[Y1_idx] <- truncnorm::rtruncnorm(
          length(Y1_idx),
          a = 0,
          b = Inf,
          mean = X_Main_obs[Y1_idx, , drop = FALSE] %*% beta_full_prev[[k]] + phiY_person[Y1_idx]
        )
      }
      if (length(Y0_idx) > 0) {
        Q_new[Y0_idx] <- truncnorm::rtruncnorm(
          length(Y0_idx),
          a = -Inf,
          b = 0,
          mean = X_Main_obs[Y0_idx, , drop = FALSE] %*% beta_full_prev[[k]] + phiY_person[Y0_idx]
        )
      }
    }
  }

  Q_new
}

update_outcome_random_effects_original <- function(Q_new, X_Main_obs, beta_full_prev, strata_person, strata_fac, provider_fac, person_fac, tau2_prev, Nj, K) {
  phi_draw <- rep(0, length(provider_fac))

  for (k in seq_len(K)) {
    people_in_k <- strata_person == k
    fac_in_k <- provider_fac[strata_fac == k]
    num_fac_in_k <- length(fac_in_k)
    peop_fac <- person_fac[people_in_k]

    if (num_fac_in_k > 0) {
      Njk <- Nj[fac_in_k]
      Q_phi <- Q_new[people_in_k] - X_Main_obs[people_in_k, , drop = FALSE] %*% beta_full_prev[[k]]
      var_phi <- tau2_prev / (Njk * tau2_prev + 1)
      mean_den_phi <- Njk + (1 / tau2_prev)
      mean_num_phi <- tapply(Q_phi, peop_fac, sum)
      mu_phi <- matrix(mean_num_phi / mean_den_phi, nrow = num_fac_in_k, ncol = 1)
      phi_draw[fac_in_k] <- MASS::mvrnorm(1, mu = mu_phi, Sigma = diag(var_phi, num_fac_in_k))
    }
  }

  phi_draw
}

update_outcome_noncovariates_original <- function(Q_new, X_Main_obs, phiY_person, strata_person, beta_full_prev, noncovsX, covsX, v_beta_k_inv, m_beta_k, K) {
  muY_new <- rep(NA_real_, K)
  delta0_new <- rep(NA_real_, K)
  delta1_new <- rep(NA_real_, K)

  for (k in seq_len(K)) {
    people_in_k <- strata_person == k
    num_in_k <- sum(people_in_k)

    if (num_in_k > 0) {
      X_in_k <- X_Main_obs[people_in_k, noncovsX, drop = FALSE]
      X_in_k_t <- t(X_in_k)
      vXkXk_inv <- chol2inv(chol(X_in_k_t %*% X_in_k + v_beta_k_inv))
      Xc_in_k <- X_Main_obs[people_in_k, covsX, drop = FALSE]
      Betac <- beta_full_prev[[k]][covsX]
      Betak_hat <- vXkXk_inv %*% (X_in_k_t %*% (Q_new[people_in_k] - Xc_in_k %*% Betac - phiY_person[people_in_k]) +
        v_beta_k_inv %*% m_beta_k)
      Sigmak_hat <- vXkXk_inv
    } else {
      Betak_hat <- m_beta_k
      Sigmak_hat <- solve(v_beta_k_inv)
    }

    muy_deltak_new <- MASS::mvrnorm(1, mu = Betak_hat, Sigma = Sigmak_hat)
    muY_new[k] <- muy_deltak_new[1]
    delta0_new[k] <- muy_deltak_new[2]
    delta1_new[k] <- muy_deltak_new[3]
  }

  list(muY_draw = muY_new, delta0_draw = delta0_new, delta1_draw = delta1_new)
}

update_outcome_shared_covariates_original <- function(Q_new, X_Main_obs, phiY_person, strata_person, muY_draw, delta0_draw, delta1_draw, covsX, noncovsX, v_beta_inv, m_beta) {
  X_in_k <- X_Main_obs[, covsX, drop = FALSE]
  X_in_k_t <- t(X_in_k)
  vXkXk_inv <- chol2inv(chol(X_in_k_t %*% X_in_k + v_beta_inv))
  Xc_in_k <- X_Main_obs[, noncovsX, drop = FALSE]
  Betac <- rbind(muY_draw[strata_person], delta0_draw[strata_person], delta1_draw[strata_person])
  Betak_hat <- vXkXk_inv %*% (X_in_k_t %*% (Q_new - diag(Xc_in_k %*% Betac) - phiY_person) +
    v_beta_inv %*% m_beta)
  Sigmak_hat <- vXkXk_inv
  MASS::mvrnorm(1, mu = Betak_hat, Sigma = Sigmak_hat)
}

assemble_outcome_betas_original <- function(muY_draw, beta_covs_draw, delta0_draw, delta1_draw, YCN, K) {
  beta_full_draw <- vector("list", K)
  for (k in seq_len(K)) {
    beta_full_draw[[k]] <- c(
      muY_draw[k],
      beta_covs_draw[1:YCN],
      delta0_draw[k],
      beta_covs_draw[(YCN + 1):(2 * YCN)],
      delta1_draw[k]
    )
  }
  beta_full_draw
}

update_strata_membership_original <- function(CZ_i, Y_current, X_Main_obs, D_current, X_D, mu_draws, Sigma_new, beta_full_draw, alpha_list_draw, phiY_person, phiD_person, pi_new, person_fac, J, K, facility_keys = seq_len(J), context = NULL) {
  log_weights <- matrix(NA_real_, nrow = J, ncol = K)
  for (k in seq_len(K)) {
    log_cz <- apply(CZ_i, 1, mvtnorm::dmvnorm, mean = mu_draws[, k],
                    sigma = Sigma_new, log = TRUE)
    eta_y <- X_Main_obs %*% beta_full_draw[[k]] + phiY_person
    log_y <- stats::pnorm((2 * Y_current - 1) * eta_y, log.p = TRUE)
    eta_d <- X_D %*% alpha_list_draw[[k]] + phiD_person
    log_d <- stats::pnorm((2 * D_current - 1) * eta_d, log.p = TRUE)

    # The facility likelihood is the product over its residents.
    log_y_fac <- rowsum(as.numeric(log_y), person_fac, reorder = TRUE)[, 1L]
    log_d_fac <- rowsum(as.numeric(log_d), person_fac, reorder = TRUE)[, 1L]
    log_weights[, k] <- log_cz + log_y_fac + log_d_fac + log(pi_new[k])
  }

  row_max <- apply(log_weights, 1, max)
  bad <- which(!is.finite(row_max))
  if (length(bad)) {
    where <- if (is.null(context)) "" else paste0(" (", context, ")")
    stop("Cannot normalize latent-stratum weights", where, " for facility ",
         paste(facility_keys[bad], collapse = ", "), ".", call. = FALSE)
  }

  Sk_probs <- exp(log_weights - row_max)
  Sk_probs <- Sk_probs / rowSums(Sk_probs)
  S_new <- rep(NA_integer_, J)
  for (j in seq_len(J)) {
    S_new[j] <- which(stats::rmultinom(1, 1, Sk_probs[j, ]) == 1)
  }

  list(strata_fac_draw = S_new, strata_prob_draw = Sk_probs)
}

impute_control_c_original <- function(C_Sim_prev, Z_Matrix, mu_draws, Sigma_new, strata_fac_new, treat_fac_bool, L, P, K) {
  SigmaCC <- as.matrix(Sigma_new[1:L, 1:L, drop = FALSE])
  SigmaCZ <- matrix(Sigma_new[1:L, (L + 1):(L + P), drop = FALSE], nrow = L, ncol = P)
  SigmaZC <- t(SigmaCZ)
  SigmaZZ <- as.matrix(Sigma_new[(L + 1):(L + P), (L + 1):(L + P), drop = FALSE])
  SigmaZZ_inv <- solve(SigmaZZ)
  C_draw <- C_Sim_prev

  for (k in seq_len(K)) {
    control_and_k <- strata_fac_new == k & !treat_fac_bool
    J_k_cont <- sum(control_and_k)
    if (J_k_cont > 0) {
      mu_Ck <- mu_draws[1:L, k]
      mu_Zk <- mu_draws[(L + 1):(L + P), k]

      Ck_hat <- mu_Ck + SigmaCZ %*% SigmaZZ_inv %*% t(as.matrix(sweep(as.matrix(Z_Matrix[control_and_k, , drop = FALSE]), 2, mu_Zk)))
      Ck_hat_list <- lapply(seq_len(ncol(Ck_hat)), function(i) Ck_hat[, i])
      Sigma_mu_C_hat <- SigmaCC - SigmaCZ %*% SigmaZZ_inv %*% SigmaZC
      sampled_Cs <- matrix(unlist(lapply(Ck_hat_list, MASS::mvrnorm, n = 1, Sigma = Sigma_mu_C_hat)),
        nrow = J_k_cont, ncol = L, byrow = TRUE)

      for (l in seq_len(L)) {
        C_draw[[l]][control_and_k] <- sampled_Cs[, l]
      }
    }
  }

  C_draw
}

impute_unobserved_d_original <- function(Y_current, X_Main_comp, X_Main_nocomp,
                                         beta_full_draw, phiY_person, phiD_person, X_D,
                                         alpha_list_draw, strata_person, D_missing) {
  D_samp_probs <- D_new <- rep(NA_real_, length(D_missing))
  for (k in seq_along(beta_full_draw)) {
    idx <- which(strata_person == k & D_missing)
    if (!length(idx)) next

    # Evaluate both compliance types under the resident's actual assignment.
    eta_y1 <- X_Main_comp[idx, , drop = FALSE] %*% beta_full_draw[[k]] + phiY_person[idx]
    eta_y0 <- X_Main_nocomp[idx, , drop = FALSE] %*% beta_full_draw[[k]] + phiY_person[idx]
    eta_d <- X_D[idx, , drop = FALSE] %*% alpha_list_draw[[k]] + phiD_person[idx]
    # Use the signed probit predictor to evaluate either binary outcome.
    sign_y <- 2 * Y_current[idx] - 1
    log_y1 <- stats::pnorm(sign_y * eta_y1, log.p = TRUE)
    log_y0 <- stats::pnorm(sign_y * eta_y0, log.p = TRUE)
    log_w1 <- log_y1 + stats::pnorm(eta_d, log.p = TRUE)
    log_w0 <- log_y0 + stats::pnorm(eta_d, lower.tail = FALSE, log.p = TRUE)
    prob <- stats::plogis(log_w1 - log_w0)
    if (any(!is.finite(prob))) {
      stop("Cannot normalize missing-compliance weights in stratum ", k, call. = FALSE)
    }
    D_samp_probs[idx] <- prob
    D_new[idx] <- stats::rbinom(length(idx), 1, prob)
  }
  list(D_draw = D_new, D_prob_draw = D_samp_probs)
}

rebuild_outcome_design_original <- function(X, treat_person_bin, D_draw) {
  X_Main_obs <- cbind(1, as.matrix(X))
  X_Main_obs <- cbind(X_Main_obs, (1 - treat_person_bin) * D_draw * 1)
  X_Main_obs <- cbind(X_Main_obs, treat_person_bin * D_draw * cbind(as.matrix(X), 1))
  colnames(X_Main_obs) <- paste0("V", seq_len(ncol(X_Main_obs)))
  X_Main_obs
}

impute_missing_y_original <- function(Y_current, missingYs, X_Main_obs, beta_full_draw, phiY_person, strata_person, K) {
  Y_updated <- Y_current

  for (k in seq_len(K)) {
    missing_idx <- missingYs & strata_person == k
    Ymiss_mean <- stats::pnorm(X_Main_obs[missing_idx, , drop = FALSE] %*% beta_full_draw[[k]] +
      phiY_person[missing_idx])
    Y_updated[missing_idx] <- stats::rbinom(sum(missing_idx), 1, Ymiss_mean)
  }

  Y_updated
}

# Standardize the integrated probit risks over each draw's empirical strata.
# The compliance and outcome random effects are independent under the model.
compute_binary_estimands <- function(X, X_D, strata_person, pi_new, alpha_list_draw,
                                     beta_full_draw, tau2_phi_d, tau2_phi_y) {
  K <- length(alpha_list_draw)
  sd_d <- sqrt(1 + tau2_phi_d)
  sd_y <- sqrt(1 + tau2_phi_y)
  arm_itt <- arm_cace <- matrix(NA_real_, K, 2L)
  log_complier_prob <- numeric(K)

  for (k in seq_len(K)) {
    idx <- which(strata_person == k)
    if (!length(idx)) {
      stop("Cannot standardize over an empty empirical stratum: ", k, call. = FALSE)
    }
    x_k <- X[idx, , drop = FALSE]
    beta <- beta_full_draw[[k]]
    eta_d <- drop(X_D[idx, , drop = FALSE] %*% alpha_list_draw[[k]]) / sd_d
    log_p <- stats::pnorm(eta_d, log.p = TRUE)
    p <- exp(log_p)
    p0 <- stats::pnorm(eta_d, lower.tail = FALSE)

    eta0 <- drop(rebuild_outcome_design_original(x_k, 0, 0) %*% beta)
    eta1_complier <- drop(rebuild_outcome_design_original(x_k, 1, 1) %*% beta)
    eta0_complier <- drop(rebuild_outcome_design_original(x_k, 0, 1) %*% beta)
    q0 <- stats::pnorm(eta0 / sd_y)
    q_complier <- stats::pnorm(cbind(eta1_complier, eta0_complier) / sd_y)
    arm_itt[k, ] <- colMeans(p * q_complier + p0 * q0)

    # Normalize the resident weights among compliers before exponentiating.
    max_log_p <- max(log_p)
    if (!is.finite(max_log_p)) {
      stop("Cannot compute complier weights in stratum ", k, call. = FALSE)
    }
    relative_p <- exp(log_p - max_log_p)
    arm_cace[k, ] <- colSums(relative_p * q_complier) / sum(relative_p)
    log_complier_prob[k] <- max_log_p + log(mean(relative_p))
  }

  # Marginal complier risks use the mixture weights conditional on compliance.
  log_mass <- log(pi_new) + log_complier_prob
  if (!is.finite(max(log_mass))) {
    stop("Cannot compute marginal complier weights.", call. = FALSE)
  }
  mass <- exp(log_mass - max(log_mass))
  pooled_itt <- colSums(pi_new * arm_itt)
  pooled_cace <- colSums(mass * arm_cace) / sum(mass)

  estimands <- c(ITT = pooled_itt[1] - pooled_itt[2],
                 ITT_Y1 = pooled_itt[1], ITT_Y0 = pooled_itt[2])
  for (k in seq_len(K)) {
    estimands[paste0("ITT", k, c("", "_Y1", "_Y0"))] <-
      c(arm_itt[k, 1] - arm_itt[k, 2], arm_itt[k, ])
  }
  estimands <- c(estimands, CACE = pooled_cace[1] - pooled_cace[2],
                 CACE_Y1 = pooled_cace[1], CACE_Y0 = pooled_cace[2])
  for (k in seq_len(K)) {
    estimands[paste0("CACE", k, c("", "_Y1", "_Y0"))] <-
      c(arm_cace[k, 1] - arm_cace[k, 2], arm_cace[k, ])
  }
  estimands
}

simulate_ppc_original <- function(J, Facility, X, X_D, treat_person_bin, alpha_list_draw, beta_full_draw, tau2_phi_d, tau2_phi_y, pi_draw) {
  K <- length(alpha_list_draw)
  N <- nrow(X)
  PPcheckD <- rep(NA_real_, N)
  PPcheckY <- rep(NA_real_, N)
  if (K == 2L) {
    PPcheckS <- stats::rbinom(J, 1, pi_draw[2]) + 1
  } else {
    PPcheckS <- sample(seq_len(K), J, replace = TRUE, prob = pi_draw)
  }
  PPcheckS_person <- PPcheckS[Facility]

  for (k in seq_len(K)) {
    people_in_k <- PPcheckS_person == k
    PP_phiD <- stats::rnorm(J, mean = 0, sd = sqrt(tau2_phi_d))
    PP_phiY <- stats::rnorm(J, mean = 0, sd = sqrt(tau2_phi_y))
    PP_phiD_person <- PP_phiD[Facility]
    PP_phiY_person <- PP_phiY[Facility]

    PP_Dmean <- stats::pnorm(X_D[people_in_k, , drop = FALSE] %*% alpha_list_draw[[k]] + PP_phiD_person[people_in_k])
    PPcheckD[people_in_k] <- stats::rbinom(sum(people_in_k), 1, PP_Dmean)

    PP_X_Main_obsk <- cbind(1, as.matrix(X[people_in_k, , drop = FALSE]))
    PP_X_Main_obsk <- cbind(PP_X_Main_obsk, (1 - treat_person_bin[people_in_k]) * PPcheckD[people_in_k] * 1)
    PP_X_Main_obsk <- cbind(PP_X_Main_obsk, treat_person_bin[people_in_k] * PPcheckD[people_in_k] * cbind(as.matrix(X[people_in_k, , drop = FALSE]), 1))
    colnames(PP_X_Main_obsk) <- paste0("V", seq_len(ncol(PP_X_Main_obsk)))

    PP_Ymean <- stats::pnorm(PP_X_Main_obsk %*% beta_full_draw[[k]] + PP_phiY_person[people_in_k])
    PPcheckY[people_in_k] <- stats::rbinom(sum(people_in_k), 1, PP_Ymean)
  }

  list(
    PPcheckD = PPcheckD,
    PPcheckY = PPcheckY,
    PPcheckS = PPcheckS,
    PPcheckS_person = PPcheckS_person
  )
}
