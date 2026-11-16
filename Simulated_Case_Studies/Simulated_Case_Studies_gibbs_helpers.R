############################################################################
## Gibbs-sampling updates and posterior estimand calculations.
##
## These functions implement the full conditional updates used in the Bayesian
## model for latent implementation strata, individual compliance, and outcomes.
############################################################################

# Construct facility-level indexing objects used to move between cluster-level
# and individual-level quantities within the Gibbs sampler.
build_facility_index_cache <- function(person_fac, J) {
  person_fac_int <- as.integer(person_fac)
  provider_fac <- seq_len(J)
  people_by_facility <- split(seq_along(person_fac_int), person_fac_int)

  list(
    provider_fac = provider_fac,
    person_fac = person_fac_int,
    person_fac_factor = factor(person_fac_int, levels = provider_fac),
    people_by_facility = people_by_facility
  )
}

# Build the components of the outcome-model design matrix that are fixed within
# a dataset. The only time-varying element is the latent compliance indicator.
build_outcome_design_cache <- function(x_matrix, treat) {
  x_matrix <- as.matrix(x_matrix)
  treat <- as.integer(treat)
  zero_d <- rep(0, nrow(x_matrix))
  x_interaction <- cbind(x_matrix, 1)
  base_columns <- cbind(1, x_matrix)
  colnames(base_columns) <- NULL
  colnames(x_interaction) <- NULL

  list(
    treat = treat,
    base_columns = base_columns,
    x_interaction = x_interaction,
    treat_zero = build_outcome_design_matrix(x_matrix, treat, zero_d),
    control_zero = build_outcome_design_matrix(x_matrix, rep(0L, length(treat)), zero_d)
  )
}

# Reconstruct the observed, treated, and control outcome-model design matrices
# after updating the latent compliance indicators.
update_outcome_designs_from_cache <- function(design_cache, d_value) {
  d_value <- as.numeric(d_value)
  interaction_term <- d_value * design_cache$x_interaction

  treat_all <- cbind(
    design_cache$base_columns,
    0 * d_value,
    interaction_term
  )
  control_current <- cbind(
    design_cache$base_columns,
    d_value,
    0 * interaction_term
  )
  observed <- control_current
  treated_rows <- design_cache$treat == 1
  observed[treated_rows, ] <- treat_all[treated_rows, , drop = FALSE]
  design_names <- paste0("V", seq_len(ncol(observed)))
  colnames(observed) <- design_names
  colnames(treat_all) <- design_names
  colnames(control_current) <- design_names
  control_comp <- cbind(
    design_cache$base_columns,
    1,
    0 * design_cache$x_interaction
  )
  control_nocomp <- design_cache$control_zero
  colnames(control_comp) <- design_names
  colnames(control_nocomp) <- design_names

  list(
    observed = observed,
    treat_all = treat_all,
    control_current = control_current,
    control_comp = control_comp,
    control_nocomp = control_nocomp
  )
}

# Aggregate individual log-likelihood contributions to the facility level,
# which is required for updating latent facility strata.
sum_by_facility <- function(values, person_fac_factor, J) {
  summed <- rowsum(values, group = person_fac_factor, reorder = FALSE)
  out <- numeric(J)
  out[as.integer(rownames(summed))] <- summed[, 1]
  out
}

# Draw from a multivariate Normal distribution when a Cholesky factorization is
# convenient for the full conditional update.
rmvnorm_chol <- function(n, mean, covariance) {
  mean <- as.numeric(mean)
  covariance <- as.matrix(covariance)
  chol_cov <- chol(covariance)
  z <- matrix(rnorm(n * length(mean)), nrow = n, ncol = length(mean))
  sweep(z %*% chol_cov, 2, mean, FUN = "+")
}

# Update the shared covariance matrix in the Gaussian mixture model for the
# latent implementation summaries and baseline cluster covariates.
update_mixture_covariance <- function(c_i, z_i, strata_fac, mu_draws, S_0, nu_0, K) {
  S_J <- S_0

  for (k in seq_len(K)) {
    fac_in_k <- strata_fac == k
    J_k <- sum(fac_in_k)

    if (J_k > 0) {
      cz_k <- rbind(t(c_i[fac_in_k, , drop = FALSE]), t(z_i[fac_in_k, , drop = FALSE]))
      cz_k_minus_mu_k <- sweep(cz_k, 1, mu_draws[, k])
      S_J <- S_J + Reduce(
        "+",
        lapply(
          seq_len(J_k),
          function(j) cz_k_minus_mu_k[, j] %*% t(cz_k_minus_mu_k[, j])
        )
      )
    }
  }

  MCMCpack::riwish(nu_0 + length(strata_fac), S_J)
}

# Update the latent-stratum mixing proportions from their Dirichlet full
# conditional.
update_mixture_weights <- function(strata_fac, lambda_pk, K) {
  counts <- tabulate(strata_fac, nbins = K)
  as.numeric(DirichletReg::rdirichlet(1, lambda_pk + counts))
}

# Update the facility-level random effects in the individual compliance model
# from their Gaussian full conditional.
update_compliance_random_effects <- function(U_new, X_D, alpha_current, strata_person, strata_fac, provider_fac, person_fac, tau2_prev, K) {
  phi_draw <- rep(0, length(provider_fac))

  for (k in seq_len(K)) {
    people_in_k <- strata_person == k
    fac_in_k <- provider_fac[strata_fac == k]
    num_fac_in_k <- length(fac_in_k)
    peop_fac <- person_fac[people_in_k]

    if (num_fac_in_k == 0) {
      next
    }

    U_phi <- U_new[people_in_k] - X_D[people_in_k, , drop = FALSE] %*% alpha_current[[k]]
    # The precision adds one unit for each resident's latent probit response.
    facility_index <- factor(peop_fac, levels = fac_in_k)
    n_k <- tabulate(facility_index, nbins = num_fac_in_k)
    var_phi <- tau2_prev / (n_k * tau2_prev + 1)
    mean_den_phi <- n_k + (1 / tau2_prev)
    mean_num_phi <- tapply(U_phi, facility_index, sum)
    mu_phi <- matrix(mean_num_phi / mean_den_phi, nrow = num_fac_in_k, ncol = 1)
    phi_draw[fac_in_k] <- MASS::mvrnorm(1, mu = mu_phi, Sigma = diag(var_phi, num_fac_in_k))
  }

  phi_draw
}

# Update the facility-level random effects in the outcome model from their
# Gaussian full conditional.
update_outcome_random_effects <- function(Y_obs, X_Main_obs, beta_current, strata_person, strata_fac, provider_fac, person_fac, sigma2_prev, tau2_prev, K) {
  phi_draw <- rep(0, length(provider_fac))

  for (k in seq_len(K)) {
    people_in_k <- strata_person == k
    fac_in_k <- provider_fac[strata_fac == k]
    num_fac_in_k <- length(fac_in_k)
    peop_fac <- person_fac[people_in_k]

    if (num_fac_in_k == 0) {
      next
    }

    Y_phi <- Y_obs[people_in_k] - X_Main_obs[people_in_k, , drop = FALSE] %*% beta_current[[k]]
    # Each resident contributes the stratum-specific outcome precision.
    facility_index <- factor(peop_fac, levels = fac_in_k)
    n_k <- tabulate(facility_index, nbins = num_fac_in_k)
    var_phi <- (sigma2_prev[k] * tau2_prev) / (n_k * tau2_prev + sigma2_prev[k])
    mean_den_phi <- n_k + (sigma2_prev[k] / tau2_prev)
    mean_num_phi <- tapply(Y_phi, facility_index, sum)
    mu_phi <- matrix(mean_num_phi / mean_den_phi, nrow = num_fac_in_k, ncol = 1)
    phi_draw[fac_in_k] <- MASS::mvrnorm(1, mu = mu_phi, Sigma = diag(var_phi, num_fac_in_k))
  }

  phi_draw
}

# Update the stratum-specific regression coefficients in the probit compliance
# model.
update_compliance_coefficients <- function(X_D, U_new, phiD_person, strata_person, v_alpha_inv, m_alpha, K) {
  alpha_draws <- vector("list", K)
  prior_cov <- solve(v_alpha_inv)

  for (k in seq_len(K)) {
    people_in_k <- strata_person == k

    if (any(people_in_k)) {
      X_in_k <- X_D[people_in_k, , drop = FALSE]
      rhs <- U_new[people_in_k] - phiD_person[people_in_k]
      posterior_precision <- crossprod(X_in_k) + v_alpha_inv
      posterior_cov <- chol2inv(chol(posterior_precision))
      posterior_mean <- posterior_cov %*% (crossprod(X_in_k, rhs) + v_alpha_inv %*% m_alpha)
    } else {
      posterior_cov <- prior_cov
      posterior_mean <- m_alpha
    }

    alpha_draws[[k]] <- as.numeric(MASS::mvrnorm(1, mu = posterior_mean, Sigma = posterior_cov))
  }

  alpha_draws
}

# Update the stratum-specific regression coefficients in the outcome model.
update_outcome_coefficients <- function(Y_obs, X_Main_obs, phiY_person, strata_person, sigma2_current, v_beta_inv, m_beta, K) {
  beta_draws <- vector("list", K)

  for (k in seq_len(K)) {
    people_in_k <- strata_person == k

    if (any(people_in_k)) {
      X_in_k <- X_Main_obs[people_in_k, , drop = FALSE]
      X_in_k_t <- t(X_in_k)
      sigma_inv <- 1 / sigma2_current[k]
      posterior_cov <- chol2inv(chol(sigma_inv * X_in_k_t %*% X_in_k + v_beta_inv))
      posterior_mean <- posterior_cov %*% (
        X_in_k_t %*% (sigma_inv * (Y_obs[people_in_k] - phiY_person[people_in_k])) +
          v_beta_inv %*% m_beta
      )
    } else {
      posterior_cov <- solve(v_beta_inv)
      posterior_mean <- m_beta
    }

    beta_draws[[k]] <- as.numeric(MASS::mvrnorm(1, mu = posterior_mean, Sigma = posterior_cov))
  }

  beta_draws
}

# Update the latent facility implementation stratum by combining the
# facility-level likelihood contributions from the mixture model, compliance
# model, and outcome model.
update_strata_membership <- function(cz_i, Y_obs, X_Main_obs, D_current, X_D, mu_draws, sigma_new, beta_new, sigma2_new, alpha_new, phiY_person, phiD_person, pi_new, person_fac, K) {
  J <- nrow(cz_i)
  person_fac_factor <- factor(person_fac, levels = seq_len(J))
  Sk_probs <- matrix(NA_real_, nrow = J, ncol = K)

  for (k in seq_len(K)) {
    CZ_dens <- mvtnorm::dmvnorm(cz_i, mean = mu_draws[, k], sigma = sigma_new, log = TRUE)
    Y_mean_k <- drop(X_Main_obs %*% beta_new[[k]] + phiY_person)
    Y_dens_fac <- sum_by_facility(
      dnorm(Y_obs, Y_mean_k, sqrt(sigma2_new[k]), log = TRUE),
      person_fac_factor = person_fac_factor,
      J = J
    )

    D_dens_fac <- sum_by_facility(
      dbinom(
        D_current,
        1,
        pnorm(drop(X_D %*% alpha_new[[k]] + phiD_person)),
        log = TRUE
      ),
      person_fac_factor = person_fac_factor,
      J = J
    )

    Sk_probs[, k] <- CZ_dens + Y_dens_fac + D_dens_fac + log(pi_new[k])
  }

  row_max <- apply(Sk_probs, 1, max)
  Sk_probs <- exp(Sk_probs - row_max)
  Sk_probs <- Sk_probs / rowSums(Sk_probs)

  strata_fac_new <- rep(NA_integer_, J)
  for (j in seq_len(J)) {
    strata_fac_new[j] <- which(rmultinom(1, 1, Sk_probs[j, ]) == 1)
  }

  list(strata_fac = strata_fac_new, probabilities = Sk_probs)
}

# Impute the unobserved cluster-level implementation summaries in the control
# arm from the conditional Gaussian distribution implied by the mixture model.
impute_control_compliance_measures <- function(C_current, Z_matrix, mu_draws, sigma_new, strata_fac_new, treat_fac_bool, L, P, K) {
  C_updated <- C_current
  SigmaCC <- as.matrix(sigma_new[1:L, 1:L, drop = FALSE])
  SigmaCZ <- as.matrix(sigma_new[1:L, (L + 1):(L + P), drop = FALSE])
  SigmaZC <- t(SigmaCZ)
  SigmaZZ <- as.matrix(sigma_new[(L + 1):(L + P), (L + 1):(L + P), drop = FALSE])
  SigmaZZ_inv <- solve(SigmaZZ)
  conditional_cov <- SigmaCC - SigmaCZ %*% SigmaZZ_inv %*% SigmaZC

  for (k in seq_len(K)) {
    control_and_k <- which(strata_fac_new == k & !treat_fac_bool)
    if (length(control_and_k) == 0) {
      next
    }

    mu_Ck <- as.numeric(mu_draws[1:L, k])
    mu_Zk <- as.numeric(mu_draws[(L + 1):(L + P), k])
    centered_z <- sweep(as.matrix(Z_matrix[control_and_k, , drop = FALSE]), 2, mu_Zk, FUN = "-")
    conditional_mean <- t(mu_Ck + SigmaCZ %*% SigmaZZ_inv %*% t(centered_z))
    conditional_mean_list <- lapply(seq_len(nrow(conditional_mean)), function(i) conditional_mean[i, ])
    sampled_Cs <- matrix(
      unlist(lapply(conditional_mean_list, MASS::mvrnorm, n = 1, Sigma = conditional_cov)),
      nrow = length(control_and_k),
      ncol = L,
      byrow = TRUE
    )
    C_updated[control_and_k, seq_len(L)] <- sampled_Cs
  }

  C_updated
}

# Compute one posterior draw of the super-population ITT and CACE estimands,
# both marginally and within the ordered latent implementation strata.
compute_posterior_draw_estimands <- function(X, X_D, alpha_draws, beta1_draws, delta0_draws, delta1_draws, pi_draw, tau2_phi_d_draw, mu_draws) {
  K <- length(alpha_draws)
  X_effect <- as.matrix(cbind(X, 1))
  p_ijk <- matrix(NA_real_, nrow = nrow(X_D), ncol = K)
  delta_k <- delta1_draws - delta0_draws
  beta_diff_phi <- matrix(NA_real_, nrow = ncol(X_effect), ncol = K)

  for (k in seq_len(K)) {
    p_ijk[, k] <- pnorm(drop(X_D %*% alpha_draws[[k]]) / sqrt(tau2_phi_d_draw + 1))
    beta_diff_phi[, k] <- c(beta1_draws[[k]], delta_k[k])
  }

  treatment_effect_k <- X_effect %*% beta_diff_phi
  weighted_treatment_effect_k <- treatment_effect_k * p_ijk
  itt_by_stratum <- colMeans(weighted_treatment_effect_k)
  itt_overall <- mean(weighted_treatment_effect_k %*% pi_draw)

  sum_p_k <- colSums(p_ijk)
  cace_by_stratum <- colSums(weighted_treatment_effect_k) / sum_p_k

  # Average the stratum-specific contrasts over the joint distribution of X
  # and stratum among compliers. Exclusion makes the numerator equal to ITT.
  pooled_complier_prob <- sum(pi_draw * colMeans(p_ijk))
  cace_overall <- itt_overall / pooled_complier_prob

  ordered_strata <- identify_ordered_strata(mu_draws)

  c(
    ITT = itt_overall,
    CACE = cace_overall,
    ITT1 = itt_by_stratum[ordered_strata$lower],
    ITT2 = itt_by_stratum[ordered_strata$upper],
    CACE1 = cace_by_stratum[ordered_strata$lower],
    CACE2 = cace_by_stratum[ordered_strata$upper]
  )
}
