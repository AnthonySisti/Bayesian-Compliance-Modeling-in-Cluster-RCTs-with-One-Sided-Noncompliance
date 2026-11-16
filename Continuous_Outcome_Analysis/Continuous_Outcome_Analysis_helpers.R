# Supplementary Materials Helper File: multichain continuous-outcome analysis.
#
# This file contains the workflow-level functions used by the main analysis
# script. The Gibbs updating kernels are sourced from the model-helper file;
# this file organizes those updates into a complete Bayesian analysis with
# retained posterior draws, multichain execution, post hoc relabeling of the
# latent strata, MCMC diagnostics, and posterior predictive checks.
#
# Conceptually, the workflow proceeds in five stages:
# 1. define the Gibbs-sampler configuration,
# 2. run a single posterior chain under the continuous outcome model,
# 3. relabel the latent mixture strata in a deterministic reporting convention,
# 4. pool retained draws across chains and compute convergence diagnostics,
# 5. evaluate posterior predictive checks on the combined retained draws.

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
source(file.path(script_dir, "Continuous_Outcome_Analysis_model_helpers.R"))
load_continuous_outcome_packages()

# Analysis configuration.
# `burn`, `thin`, and `n_imps` determine the total number of Gibbs iterations
# and which posterior draws are retained for inference. `K = 2` fixes the
# latent compliance-mixture model to two implementation strata, matching the
# analysis described in the manuscript.
make_continuous_outcome_config <- function(
    seed = 20260422L,
    burn = 500L,
    thin = 3L,
    n_imps = 1000L,
    K = 2L,
    plot_initial_clusters = FALSE,
    output_dir = file.path(script_dir, "outputs")) {
  burn <- as.integer(burn)
  thin <- as.integer(thin)
  n_imps <- as.integer(n_imps)
  Sims <- burn + thin * n_imps
  samps_use <- burn + seq_len(n_imps) * thin

  list(
    seed = as.integer(seed),
    burn = burn,
    thin = thin,
    n_imps = n_imps,
    Sims = Sims,
    samps_use = samps_use,
    K = as.integer(K),
    plot_initial_clusters = isTRUE(plot_initial_clusters),
    output_dir = output_dir
  )
}

# Posterior summaries for the primary causal estimands. These are the pooled
# summaries reported after combining the relabeled retained draws across
# chains.
summarize_primary_estimands <- function(out_df) {
  summary_df <- data.frame(
    estimand = colnames(out_df),
    mean = colMeans(out_df),
    sd = apply(out_df, 2, stats::sd),
    L_95_CrI = apply(out_df, 2, stats::quantile, 0.025),
    U_95_CrI = apply(out_df, 2, stats::quantile, 0.975),
    row.names = NULL
  )

  dplyr::mutate(summary_df, dplyr::across(where(is.numeric), ~ round(.x, 3)))
}

summarize_outcome_levels <- function(out_df) {
  summary_df <- data.frame(
    estimand = colnames(out_df),
    mean = colMeans(out_df),
    sd = apply(out_df, 2, stats::sd),
    L_95_CrI = apply(out_df, 2, stats::quantile, 0.025),
    U_95_CrI = apply(out_df, 2, stats::quantile, 0.975),
    row.names = NULL
  )

  dplyr::mutate(summary_df, dplyr::across(where(is.numeric), ~ round(.x, 3)))
}

# Multichain diagnostics for the primary estimands. Gelman-Rubin statistics
# summarize between-chain mixing, while effective sample size summarizes the
# precision of Monte Carlo estimation after accounting for within-chain
# autocorrelation.
build_estimand_mcmc_list <- function(chain_draws, estimand_name) {
  coda::mcmc.list(
    lapply(chain_draws, function(draw_obj) {
      coda::mcmc(as.matrix(draw_obj$posterior_draws$estimands_relabelled[, estimand_name, drop = FALSE]))
    })
  )
}

compute_primary_estimand_diagnostics <- function(chain_draws, estimand_names = NULL) {
  if (is.null(estimand_names)) {
    estimand_names <- colnames(chain_draws[[1]]$posterior_draws$estimands_relabelled)
  }

  diagnostic_rows <- lapply(estimand_names, function(estimand_name) {
    mcmc_list <- build_estimand_mcmc_list(chain_draws, estimand_name)
    gelman_fit <- coda::gelman.diag(mcmc_list, autoburnin = FALSE, multivariate = FALSE)
    ess_by_chain <- vapply(mcmc_list, function(chain) {
      as.numeric(coda::effectiveSize(chain))
    }, numeric(1))

    data.frame(
      estimand = estimand_name,
      rhat_point = as.numeric(gelman_fit$psrf[1, 1]),
      rhat_upper = as.numeric(gelman_fit$psrf[1, 2]),
      ess_total = sum(ess_by_chain),
      ess_min_chain = min(ess_by_chain),
      row.names = NULL
    )
  })

  diagnostic_df <- do.call(rbind, diagnostic_rows)
  dplyr::mutate(diagnostic_df, dplyr::across(where(is.numeric), ~ round(.x, 3)))
}

combine_primary_estimand_draws <- function(chain_draws) {
  do.call(
    rbind,
    lapply(chain_draws, function(draw_obj) draw_obj$posterior_draws$estimands_relabelled)
  )
}

combine_outcome_level_draws <- function(chain_draws) {
  do.call(
    rbind,
    lapply(chain_draws, function(draw_obj) draw_obj$posterior_draws$outcome_levels_relabelled)
  )
}

# The latent mixture labels are not intrinsically identified, so the retained
# draws are reordered using the first coordinate of the facility-level mixture
# mean. Stratum 1 is defined as the component with the larger first-coordinate
# mean, and stratum 2 is defined as the component with the smaller
# first-coordinate mean. This affects label-dependent quantities such as
# stratum-specific ITT and CACE, but leaves pooled label-invariant estimands
# unchanged.
identify_ordered_strata <- function(mu_draw) {
  order(mu_draw[1, ], decreasing = TRUE)
}

relabel_strata_draw <- function(strata_draw, ordered_index) {
  inverse_map <- integer(length(ordered_index))
  inverse_map[ordered_index] <- seq_along(ordered_index)
  inverse_map[as.integer(strata_draw)]
}

make_primary_estimand_names <- function(K) {
  stratum_ids <- seq_len(K)
  names <- c("ITT", paste0("ITT", stratum_ids))
  if (K >= 2L) {
    pair_index <- utils::combn(stratum_ids, 2)
    names <- c(names, apply(pair_index, 2, function(pair) paste0("ITT", pair[1], "_ITT", pair[2])))
  }
  names <- c(names, "CACE", paste0("CACE", stratum_ids))
  if (K >= 2L) {
    pair_index <- utils::combn(stratum_ids, 2)
    names <- c(names, apply(pair_index, 2, function(pair) paste0("CACE", pair[1], "_CACE", pair[2])))
  }
  names
}

make_outcome_level_names <- function(K) {
  stratum_ids <- seq_len(K)
  c(
    "ITT_Y1",
    "ITT_Y0",
    as.vector(rbind(paste0("ITT", stratum_ids, "_Y1"), paste0("ITT", stratum_ids, "_Y0"))),
    "CACE_Y1",
    "CACE_Y0",
    as.vector(rbind(paste0("CACE", stratum_ids, "_Y1"), paste0("CACE", stratum_ids, "_Y0")))
  )
}

assemble_primary_estimand_draw <- function(ITT, ITT_k, CACE, CACE_k) {
  K <- length(ITT_k)
  stratum_ids <- seq_len(K)
  values <- c(ITT = ITT)
  values <- c(values, stats::setNames(ITT_k, paste0("ITT", stratum_ids)))
  if (K >= 2L) {
    pair_index <- utils::combn(stratum_ids, 2)
    for (pair_col in seq_len(ncol(pair_index))) {
      pair <- pair_index[, pair_col]
      values[paste0("ITT", pair[1], "_ITT", pair[2])] <- ITT_k[pair[1]] - ITT_k[pair[2]]
    }
  }
  values <- c(values, CACE = CACE)
  values <- c(values, stats::setNames(CACE_k, paste0("CACE", stratum_ids)))
  if (K >= 2L) {
    pair_index <- utils::combn(stratum_ids, 2)
    for (pair_col in seq_len(ncol(pair_index))) {
      pair <- pair_index[, pair_col]
      values[paste0("CACE", pair[1], "_CACE", pair[2])] <- CACE_k[pair[1]] - CACE_k[pair[2]]
    }
  }
  values[make_primary_estimand_names(K)]
}

assemble_outcome_level_draw <- function(ITT_Y1, ITT_Y0, Y1_k, Y0_k, CACE_Y1, CACE_Y0, Y1_D_k, Y0_D_k) {
  K <- length(Y1_k)
  stratum_ids <- seq_len(K)
  values <- c(ITT_Y1 = ITT_Y1, ITT_Y0 = ITT_Y0)
  for (k in stratum_ids) {
    values[paste0("ITT", k, "_Y1")] <- Y1_k[k]
    values[paste0("ITT", k, "_Y0")] <- Y0_k[k]
  }
  values <- c(values, CACE_Y1 = CACE_Y1, CACE_Y0 = CACE_Y0)
  for (k in stratum_ids) {
    values[paste0("CACE", k, "_Y1")] <- Y1_D_k[k]
    values[paste0("CACE", k, "_Y0")] <- Y0_D_k[k]
  }
  values[make_outcome_level_names(K)]
}

combine_ppc_draws <- function(chain_draws, draw_name) {
  do.call(
    cbind,
    lapply(chain_draws, function(draw_obj) draw_obj$posterior_draws[[draw_name]])
  )
}

compute_latent_class_diagnostics <- function(chain_draws, ambiguity_thresholds = c(0.70, 0.80)) {
  posterior_draws <- lapply(chain_draws, `[[`, "posterior_draws")
  S_use_raw <- combine_ppc_draws(chain_draws, "S_use")
  ordered_index <- do.call(rbind, lapply(posterior_draws, `[[`, "ordered_index"))
  K <- ncol(ordered_index)

  S_use <- S_use_raw
  for (draw_idx in seq_len(ncol(S_use_raw))) {
    S_use[, draw_idx] <- relabel_strata_draw(S_use_raw[, draw_idx], ordered_index[draw_idx, ])
  }

  membership_prob <- sapply(seq_len(K), function(k) rowMeans(S_use == k, na.rm = TRUE))
  if (K == 1L) {
    membership_prob <- matrix(membership_prob, ncol = 1L)
  }
  colnames(membership_prob) <- paste0("Pr_S", seq_len(K))

  modal_class <- max.col(membership_prob, ties.method = "first")
  max_membership_prob <- apply(membership_prob, 1, max)

  facility_membership <- data.frame(
    facility = seq_len(nrow(membership_prob)),
    modal_stratum = modal_class,
    max_posterior_probability = max_membership_prob,
    membership_prob,
    row.names = NULL,
    check.names = FALSE
  )
  facility_membership <- dplyr::mutate(
    facility_membership,
    dplyr::across(where(is.numeric), ~ round(.x, 3))
  )

  ambiguous_membership <- do.call(
    rbind,
    lapply(ambiguity_thresholds, function(threshold) {
      n_ambiguous <- sum(max_membership_prob < threshold, na.rm = TRUE)
      data.frame(
        threshold = threshold,
        n_ambiguous = n_ambiguous,
        prop_ambiguous = n_ambiguous / length(max_membership_prob),
        row.names = NULL
      )
    })
  )
  ambiguous_membership <- dplyr::mutate(
    ambiguous_membership,
    dplyr::across(where(is.numeric), ~ round(.x, 3))
  )

  class_size_draws <- sapply(seq_len(K), function(k) colSums(S_use == k, na.rm = TRUE))
  if (K == 1L) {
    class_size_draws <- matrix(class_size_draws, ncol = 1L)
  }
  colnames(class_size_draws) <- paste0("S", seq_len(K))

  class_size_summary <- data.frame(
    stratum = paste0("S", seq_len(K)),
    mean = colMeans(class_size_draws),
    sd = apply(class_size_draws, 2, stats::sd),
    L_95_CrI = apply(class_size_draws, 2, stats::quantile, 0.025),
    U_95_CrI = apply(class_size_draws, 2, stats::quantile, 0.975),
    row.names = NULL
  )
  class_size_summary <- dplyr::mutate(
    class_size_summary,
    dplyr::across(where(is.numeric), ~ round(.x, 3))
  )

  list(
    facility_membership = facility_membership,
    ambiguous_membership = ambiguous_membership,
    class_size_summary = class_size_summary,
    class_size_draws = class_size_draws,
    relabelled_strata_draws = S_use
  )
}

compute_relabeling_stability <- function(chain_draws, near_zero_threshold = 0.05) {
  posterior_draws <- lapply(chain_draws, `[[`, "posterior_draws")
  ordered_index <- do.call(rbind, lapply(posterior_draws, `[[`, "ordered_index"))
  K <- ncol(ordered_index)

  if (K == 1L) {
    return(data.frame(
      margin = "not applicable",
      mean = NA_real_,
      sd = NA_real_,
      L_95_CrI = NA_real_,
      median = NA_real_,
      U_95_CrI = NA_real_,
      prop_below_0.05 = NA_real_,
      interpretation = "No label switching diagnostic is needed when K = 1.",
      check.names = FALSE
    ))
  }

  margin_list <- vector("list", length(posterior_draws))
  for (chain_idx in seq_along(posterior_draws)) {
    mu_use <- posterior_draws[[chain_idx]]$mu_use
    chain_order <- posterior_draws[[chain_idx]]$ordered_index
    chain_margins <- matrix(NA_real_, nrow = dim(mu_use)[3], ncol = K - 1L)

    for (draw_idx in seq_len(dim(mu_use)[3])) {
      ordered_first_coordinate <- mu_use[1, chain_order[draw_idx, ], draw_idx]
      chain_margins[draw_idx, ] <- ordered_first_coordinate[-K] - ordered_first_coordinate[-1L]
    }

    margin_list[[chain_idx]] <- chain_margins
  }

  margins <- do.call(rbind, margin_list)
  colnames(margins) <- paste0("ordered_gap_", seq_len(K - 1L))

  result <- data.frame(
    margin = colnames(margins),
    mean = colMeans(margins, na.rm = TRUE),
    sd = apply(margins, 2, stats::sd, na.rm = TRUE),
    L_95_CrI = apply(margins, 2, stats::quantile, probs = 0.025, na.rm = TRUE),
    median = apply(margins, 2, stats::quantile, probs = 0.50, na.rm = TRUE),
    U_95_CrI = apply(margins, 2, stats::quantile, probs = 0.975, na.rm = TRUE),
    prop_below_0.05 = colMeans(margins < near_zero_threshold, na.rm = TRUE),
    row.names = NULL,
    check.names = FALSE
  )

  dplyr::mutate(result, dplyr::across(where(is.numeric), ~ round(.x, 3)))
}

compute_stratified_ppc_summary <- function(
    outcome_label,
    compliance_label,
    PPcheckY,
    PPcheckD,
    PPcheckS_person,
    Y_obs_use,
    D_use,
    S_use,
    treat_Person_bool,
    Facility,
    K) {
  ppc_draw_rows <- vector("list", K)
  ppc_result_rows <- vector("list", K)
  observed_rows <- vector("list", K)

  for (k in seq_len(K)) {
    PP_signal <- PP_noise <- Obs_signal <- Obs_noise <- rep(NA_real_, ncol(PPcheckY))

    for (i in seq_len(ncol(PPcheckY))) {
      pp_treat <- PPcheckY[treat_Person_bool & PPcheckD[, i] == 1 & PPcheckS_person[, i] == k, i]
      pp_control <- PPcheckY[!treat_Person_bool & PPcheckD[, i] == 1 & PPcheckS_person[, i] == k, i]
      obs_treat <- Y_obs_use[treat_Person_bool & D_use[, i] == 1 & S_use[, i][Facility] == k, i]
      obs_control <- Y_obs_use[!treat_Person_bool & D_use[, i] == 1 & S_use[, i][Facility] == k, i]

      PP_signal[i] <- abs(mean(pp_treat) - mean(pp_control))
      Obs_signal[i] <- abs(mean(obs_treat) - mean(obs_control))
      PP_noise[i] <- sqrt(stats::var(pp_treat) / length(pp_treat) + stats::var(pp_control) / length(pp_control))
      Obs_noise[i] <- sqrt(stats::var(obs_treat) / length(obs_treat) + stats::var(obs_control) / length(obs_control))
    }

    stratum_label <- paste0("S=", k)
    ppc_result_rows[[k]] <- data.frame(
      Outcome = outcome_label,
      `Compliance Level` = compliance_label,
      Stratum = stratum_label,
      Signal = mean(PP_signal > Obs_signal, na.rm = TRUE),
      Noise = mean(PP_noise > Obs_noise, na.rm = TRUE),
      `Signal/Noise` = mean(PP_signal / PP_noise > Obs_signal / Obs_noise, na.rm = TRUE),
      check.names = FALSE
    )

    ppc_draw_rows[[k]] <- data.frame(
      Stratum = stratum_label,
      draw = seq_along(PP_signal),
      PP_signal = PP_signal,
      PP_noise = PP_noise,
      Obs_signal = Obs_signal,
      Obs_noise = Obs_noise,
      check.names = FALSE
    )

    observed_rows[[k]] <- data.frame(
      Stratum = stratum_label,
      Obs_signal = mean(Obs_signal, na.rm = TRUE),
      Obs_noise = mean(Obs_noise, na.rm = TRUE),
      check.names = FALSE
    )
  }

  list(
    pp_check_results = do.call(rbind, ppc_result_rows),
    pp_check_draws = do.call(rbind, ppc_draw_rows),
    pp_check_observed_means = do.call(rbind, observed_rows)
  )
}

summarize_ppc_discrepancy <- function(component, measure, stratum, pp_values, observed_values) {
  data.frame(
    component = component,
    measure = measure,
    stratum = stratum,
    pp_mean = mean(pp_values, na.rm = TRUE),
    observed_mean = mean(observed_values, na.rm = TRUE),
    pp_gt_observed = mean(pp_values > observed_values, na.rm = TRUE),
    row.names = NULL,
    check.names = FALSE
  )
}

compute_expanded_ppc_results <- function(chain_draws, outcome_label) {
  posterior_draws <- lapply(chain_draws, `[[`, "posterior_draws")
  fixed_data <- chain_draws[[1]]$fixed_data
  treat_Person_bool <- fixed_data$treat_Person_bool
  Facility <- fixed_data$Facility
  CZ_Matrix <- fixed_data$CZ_Matrix
  treat_Fac_bool <- fixed_data$treat_Fac_bool

  PPcheckD <- combine_ppc_draws(chain_draws, "PPcheckD")
  PPcheckY <- combine_ppc_draws(chain_draws, "PPcheckY")
  PPcheckS_raw <- combine_ppc_draws(chain_draws, "PPcheckS")
  PPcheckS_person_raw <- combine_ppc_draws(chain_draws, "PPcheckS_person")
  D_use <- combine_ppc_draws(chain_draws, "D_use")
  Y_obs_use <- combine_ppc_draws(chain_draws, "Y_obs_use")
  S_use_raw <- combine_ppc_draws(chain_draws, "S_use")
  ordered_index <- do.call(rbind, lapply(posterior_draws, `[[`, "ordered_index"))
  K <- ncol(ordered_index)

  PPcheckS <- PPcheckS_raw
  PPcheckS_person <- PPcheckS_person_raw
  S_use <- S_use_raw
  for (draw_idx in seq_len(ncol(S_use_raw))) {
    PPcheckS[, draw_idx] <- relabel_strata_draw(PPcheckS_raw[, draw_idx], ordered_index[draw_idx, ])
    PPcheckS_person[, draw_idx] <- relabel_strata_draw(PPcheckS_person_raw[, draw_idx], ordered_index[draw_idx, ])
    S_use[, draw_idx] <- relabel_strata_draw(S_use_raw[, draw_idx], ordered_index[draw_idx, ])
  }

  rows <- list()
  add_row <- function(component, measure, stratum, pp_values, observed_values) {
    rows[[length(rows) + 1L]] <<- summarize_ppc_discrepancy(component, measure, stratum, pp_values, observed_values)
  }

  add_row("Compliance", "person-level compliance rate", "All", colMeans(PPcheckD, na.rm = TRUE), colMeans(D_use, na.rm = TRUE))
  add_row("Outcome", paste0(outcome_label, " mean"), "All", colMeans(PPcheckY, na.rm = TRUE), colMeans(Y_obs_use, na.rm = TRUE))
  add_row("Outcome", paste0(outcome_label, " variance"), "All", apply(PPcheckY, 2, stats::var, na.rm = TRUE), apply(Y_obs_use, 2, stats::var, na.rm = TRUE))

  for (k in seq_len(K)) {
    person_stratum_pp <- PPcheckS_person == k
    person_stratum_obs <- S_use[Facility, , drop = FALSE] == k
    add_row(
      "Compliance",
      "person-level compliance rate",
      paste0("S=", k),
      sapply(seq_len(ncol(PPcheckD)), function(i) mean(PPcheckD[person_stratum_pp[, i], i], na.rm = TRUE)),
      sapply(seq_len(ncol(D_use)), function(i) mean(D_use[person_stratum_obs[, i], i], na.rm = TRUE))
    )
    add_row(
      "Outcome",
      paste0(outcome_label, " mean"),
      paste0("S=", k),
      sapply(seq_len(ncol(PPcheckY)), function(i) mean(PPcheckY[person_stratum_pp[, i], i], na.rm = TRUE)),
      sapply(seq_len(ncol(Y_obs_use)), function(i) mean(Y_obs_use[person_stratum_obs[, i], i], na.rm = TRUE))
    )
    add_row(
      "Outcome",
      paste0(outcome_label, " variance"),
      paste0("S=", k),
      sapply(seq_len(ncol(PPcheckY)), function(i) stats::var(PPcheckY[person_stratum_pp[, i], i], na.rm = TRUE)),
      sapply(seq_len(ncol(Y_obs_use)), function(i) stats::var(Y_obs_use[person_stratum_obs[, i], i], na.rm = TRUE))
    )
    add_row(
      "Latent class size",
      "number of facilities",
      paste0("S=", k),
      colSums(PPcheckS == k, na.rm = TRUE),
      colSums(S_use == k, na.rm = TRUE)
    )
  }

  impl_rows <- list()
  draw_offset <- 0L
  for (chain_idx in seq_along(posterior_draws)) {
    mu_use <- posterior_draws[[chain_idx]]$mu_use
    Sigma_use <- posterior_draws[[chain_idx]]$Sigma_use
    chain_order <- posterior_draws[[chain_idx]]$ordered_index
    chain_S_use <- posterior_draws[[chain_idx]]$S_use
    for (draw_idx in seq_len(dim(mu_use)[3])) {
      global_draw_idx <- draw_offset + draw_idx
      replicated_CZ <- matrix(NA_real_, nrow = nrow(CZ_Matrix), ncol = ncol(CZ_Matrix))
      for (k in seq_len(K)) {
        facilities_in_k <- PPcheckS[, global_draw_idx] == k
        if (any(facilities_in_k)) {
          source_k <- chain_order[draw_idx, k]
          replicated_CZ[facilities_in_k, ] <- mvtnorm::rmvnorm(
            sum(facilities_in_k),
            mean = mu_use[, source_k, draw_idx],
            sigma = Sigma_use[[draw_idx]]
          )
        }
      }
      for (k in seq_len(K)) {
        observed_facilities <- relabel_strata_draw(chain_S_use[, draw_idx], chain_order[draw_idx, ]) == k
        replicated_facilities <- PPcheckS[, global_draw_idx] == k
        impl_rows[[length(impl_rows) + 1L]] <- data.frame(
          draw = global_draw_idx,
          stratum = paste0("S=", k),
          PP_C_mean_treated = mean(replicated_CZ[replicated_facilities & treat_Fac_bool, 1], na.rm = TRUE),
          Obs_C_mean_treated = mean(CZ_Matrix[observed_facilities & treat_Fac_bool, 1], na.rm = TRUE),
          PP_Z_mean = mean(replicated_CZ[replicated_facilities, 2], na.rm = TRUE),
          Obs_Z_mean = mean(CZ_Matrix[observed_facilities, 2], na.rm = TRUE),
          check.names = FALSE
        )
      }
    }
    draw_offset <- draw_offset + dim(mu_use)[3]
  }

  impl_draws <- do.call(rbind, impl_rows)
  for (k in seq_len(K)) {
    stratum_draws <- impl_draws[impl_draws$stratum == paste0("S=", k), ]
    add_row("Implementation measure", "treated-facility C mean", paste0("S=", k), stratum_draws$PP_C_mean_treated, stratum_draws$Obs_C_mean_treated)
    add_row("Implementation measure", "facility Z mean", paste0("S=", k), stratum_draws$PP_Z_mean, stratum_draws$Obs_Z_mean)
  }

  list(
    summary = dplyr::mutate(do.call(rbind, rows), dplyr::across(where(is.numeric), ~ round(.x, 3))),
    implementation_draws = impl_draws
  )
}

# Posterior predictive checks are computed after combining retained draws from
# all chains. The summary statistics follow the original analysis script: they
# compare posterior predictive and observed signal, noise, and signal-to-noise
# ratios within each latent compliance stratum.
compute_combined_ppc_results <- function(chain_draws) {
  first_chain <- chain_draws[[1]]
  posterior_draws <- lapply(chain_draws, `[[`, "posterior_draws")
  fixed_data <- first_chain$fixed_data

  PPcheckD <- combine_ppc_draws(chain_draws, "PPcheckD")
  PPcheckY <- combine_ppc_draws(chain_draws, "PPcheckY")
  PPcheckS_person_raw <- combine_ppc_draws(chain_draws, "PPcheckS_person")
  D_use <- combine_ppc_draws(chain_draws, "D_use")
  Y_obs_use <- combine_ppc_draws(chain_draws, "Y_obs_use")
  S_use_raw <- combine_ppc_draws(chain_draws, "S_use")
  ordered_index <- do.call(rbind, lapply(posterior_draws, `[[`, "ordered_index"))

  PPcheckS_person <- PPcheckS_person_raw
  S_use <- S_use_raw
  for (draw_idx in seq_len(ncol(PPcheckS_person_raw))) {
    PPcheckS_person[, draw_idx] <- relabel_strata_draw(PPcheckS_person_raw[, draw_idx], ordered_index[draw_idx, ])
    S_use[, draw_idx] <- relabel_strata_draw(S_use_raw[, draw_idx], ordered_index[draw_idx, ])
  }

  treat_Person_bool <- fixed_data$treat_Person_bool
  Facility <- fixed_data$Facility
  K <- ncol(ordered_index)

  compute_stratified_ppc_summary(
    outcome_label = "Continuous",
    compliance_label = "30 Minutes",
    PPcheckY = PPcheckY,
    PPcheckD = PPcheckD,
    PPcheckS_person = PPcheckS_person,
    Y_obs_use = Y_obs_use,
    D_use = D_use,
    S_use = S_use,
    treat_Person_bool = treat_Person_bool,
    Facility = Facility,
    K = K
  )
}

# Trace plots for the primary estimands provide a direct visual assessment of
# mixing across the independently seeded chains.
make_chain_trace_pdf <- function(chain_results, file) {
  K <- ncol(chain_results[[1]]$posterior_draws$ordered_index)
  estimand_names <- c("ITT", paste0("ITT", seq_len(K)), "CACE", paste0("CACE", seq_len(K)))
  estimand_names <- estimand_names[estimand_names %in% colnames(chain_results[[1]]$posterior_draws$estimands_relabelled)]
  grDevices::pdf(file, width = 10, height = 8)
  on.exit(grDevices::dev.off(), add = TRUE)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(3, 2), mar = c(3.5, 3.5, 2.5, 1), oma = c(0, 0, 0, 0))

  for (estimand_name in estimand_names) {
    y_limits <- range(
      unlist(lapply(chain_results, function(chain_fit) {
        chain_fit$posterior_draws$estimands_relabelled[, estimand_name]
      })),
      finite = TRUE
    )

    graphics::plot(
      NA,
      xlim = c(1, nrow(chain_results[[1]]$posterior_draws$estimands_relabelled)),
      ylim = y_limits,
      xlab = "Retained draw",
      ylab = estimand_name,
      main = paste("Trace plot:", estimand_name)
    )

    for (chain_id in seq_along(chain_results)) {
      graphics::lines(
        chain_results[[chain_id]]$posterior_draws$estimands_relabelled[, estimand_name],
        col = chain_id,
        lwd = 1
      )
    }

    graphics::legend(
      "topright",
      legend = paste("Chain", seq_along(chain_results)),
      col = seq_along(chain_results),
      lty = 1,
      bty = "n",
      cex = 0.8
    )
  }

  invisible(file)
}

# Single-chain analysis function.
#
# Inputs:
# - `config`: Gibbs-sampler settings and output location.
# - `chain_id`: chain label used for organizing saved output.
# - `make_trace_pdf`: whether to save chain-level trace plots.
#
# Output:
# - A list containing retained posterior draws, posterior summaries, posterior
#   predictive summaries, and a compact representation of the terminal Gibbs
#   state.
#
# Statistical role:
# This function carries out one complete Bayesian analysis under the continuous
# outcome model. It initializes the latent-mixture model, cycles through the
# full Gibbs sampler, retains the posterior draws used for inference, and
# stores the quantities required for convergence assessment and posterior
# predictive checking.
# Align facility rows by key; assignment comes from the resident file.
align_facility_inputs <- function(cz, outcome, x) {
  if (!all(c("Facility", "C", "Z") %in% names(cz)) ||
      !all(c("Facility", "Treat", "D", "Y1", "Y0") %in% names(outcome))) {
    stop("Input files are missing required facility or resident columns.", call. = FALSE)
  }
  if (nrow(x) != nrow(outcome) || !nrow(outcome) ||
      !all(vapply(x, is.numeric, logical(1))) || any(!is.finite(as.matrix(x)))) {
    stop("Expected one complete numeric covariate row per outcome row.", call. = FALSE)
  }
  if (anyNA(cz$Facility) || anyDuplicated(cz$Facility) ||
      anyNA(outcome$Facility) || !setequal(cz$Facility, outcome$Facility)) {
    stop("Facility keys must be unique in the facility file and match the resident file.", call. = FALSE)
  }
  if (!all(outcome$Treat %in% 0:1)) {
    stop("Resident assignment must be coded 0 or 1.", call. = FALSE)
  }
  keys <- sort(cz$Facility)
  person_fac <- match(outcome$Facility, keys)
  groups <- split(outcome$Treat, factor(person_fac, levels = seq_along(keys)))
  if (any(vapply(groups, function(w) length(unique(w)) != 1L, logical(1)))) {
    stop("Assignment must be constant within each facility.", call. = FALSE)
  }
  treated <- vapply(groups, function(w) w[1L] == 1L, logical(1))
  cz <- cz[match(keys, cz$Facility), c("C", "Z"), drop = FALSE]
  rownames(cz) <- NULL
  if (!all(is.finite(cz$Z)) || !all(is.finite(cz$C[treated])) ||
      !all(is.na(cz$C[!treated]))) {
    stop("Z must be complete and C must be observed in treated facilities only.", call. = FALSE)
  }
  list(cz = cz, facility_keys = keys, person_fac = person_fac, treated = unname(treated))
}

run_continuous_outcome_single_chain <- function(
    config = make_continuous_outcome_config(),
    chain_id = 1L,
    make_trace_pdf = FALSE) {
  set.seed(config$seed)

  data_dir <- script_dir
  chain_output_dir <- file.path(config$output_dir, paste0("chain_", chain_id))
  dir.create(chain_output_dir, recursive = TRUE, showWarnings = FALSE)

  CZ_Matrix <- read.csv(file.path(data_dir, "CZ_Matrix_Synth_OG.csv"))
  X <- read.csv(file.path(data_dir, "X_Synth_OG.csv"))
  OutcomeData <- read.csv(file.path(data_dir, "Outcome_Synth_Continuous_OG.csv"))
  inputs <- align_facility_inputs(CZ_Matrix, OutcomeData, X)
  CZ_Matrix <- inputs$cz

  X <- data.matrix(X)
  X_shared <- X
  X_D <- data.matrix(cbind(Int = 1, X_shared))

  C_Matrix <- dplyr::select(CZ_Matrix, C)
  Z_Matrix <- dplyr::select(CZ_Matrix, Z)

  treat_Fac_bool <- inputs$treated
  treat_Fac_bin <- as.numeric(treat_Fac_bool)
  treat_Person_bool <- OutcomeData$Treat == 1
  treat_Person_bin <- OutcomeData$Treat
  provider_Fac <- seq_along(inputs$facility_keys)
  person_Fac <- inputs$person_fac

  K <- config$K
  N <- nrow(X)
  D <- OutcomeData$D
  Y1 <- OutcomeData$Y1
  Y0 <- OutcomeData$Y0
  Facility <- person_Fac
  L <- ncol(C_Matrix)
  P <- ncol(Z_Matrix)
  if (is.null(L)) {
    L <- 1
  }
  if (is.null(P)) {
    P <- 1
  }
  J <- nrow(CZ_Matrix)
  treat_index <- which(treat_Fac_bool)
  J_treat <- sum(treat_Fac_bool)
  control_index <- which(!treat_Fac_bool)
  J_control <- sum(!treat_Fac_bool)

  CZ <- CZ_Matrix[treat_Fac_bin == 1, ]
  if (K == 1L) {
    check_mix <- list(
      mu = list(colMeans(CZ)),
      posterior = matrix(1, nrow = nrow(CZ), ncol = 1L)
    )
    init_clust_treat <- rep(1L, nrow(CZ))
    init_clust <- rep(1L, J)
  } else {
    tot_S <- rep(0, K)
    while (min(tot_S) < 3) {
      check_mix <- try(mixtools::mvnormalmixEM(CZ, arbvar = FALSE, k = K), silent = TRUE)
      if ("try-error" %in% class(check_mix)) {
        next
      }
      init_clust_treat <- apply(check_mix$posterior, 1, which.max)
      init_clust <- rep(NA_integer_, J)
      init_clust[treat_Fac_bin == 1] <- init_clust_treat

      for (k in seq_len(K)) {
        tot_S[k] <- sum(init_clust_treat == k)
      }
    }

    cont_fac <- which(treat_Fac_bin == 0)
    for (idx in seq_along(cont_fac)) {
      z_dist <- sapply(seq_len(K), function(k) abs(CZ_Matrix[cont_fac[idx], 2] - check_mix$mu[[k]][2]))
      init_clust[cont_fac[idx]] <- which.min(z_dist)
    }
  }

  if (isTRUE(config$plot_initial_clusters)) {
    graphics::pairs(CZ, col = apply(check_mix$posterior, 1, which.max))
  }

  Sims <- config$Sims
  samps_use <- config$samps_use
  n_imps <- config$n_imps

  S_Sim <- matrix(NA_integer_, J, Sims)
  S_Sim[, 1] <- init_clust
  Strata_Fac <- S_Sim[, 1]
  Strata_Person <- Strata_Fac[person_Fac]

  pi_Sim <- matrix(NA_real_, nrow = K, ncol = Sims)
  pi_Sim[, 1] <- 1 / K

  Sigma_Sim <- vector("list", Sims)
  Sigma_Sim[[1]] <- diag(0, nrow = L + P)
  for (k in seq_len(K)) {
    Sigma_Sim[[1]] <- Sigma_Sim[[1]] + stats::cov(CZ_Matrix[treat_Fac_bin == 1 & S_Sim[, 1] == k, ])
  }

  mu_Sim <- vector("list", K)
  for (k in seq_len(K)) {
    mu_Sim[[k]] <- matrix(NA_real_, nrow = L + P, ncol = Sims)
    mu_Sim[[k]][, 1] <- colMeans(CZ_Matrix[treat_Fac_bin == 1 & S_Sim[, 1] == k, ])
  }

  C_Sim <- vector("list", L)
  for (l in seq_len(L)) {
    C_Sim[[l]] <- matrix(NA_real_, nrow = J, ncol = Sims)
    C_Sim[[l]][treat_Fac_bin == 1, 1:Sims] <- C_Matrix[treat_Fac_bin == 1, l]
    for (k in seq_len(K)) {
      C_Sim[[l]][treat_Fac_bin == 0 & init_clust == k, 1] <- mean(C_Matrix[treat_Fac_bin == 1 & init_clust == k, l])
    }
  }

  D_Sim <- matrix(NA_real_, nrow = N, ncol = Sims)
  D_Sim[treat_Person_bool, 1] <- D[treat_Person_bool]
  for (s in seq_len(Sims)) {
    D_Sim[treat_Person_bool, s] <- D[treat_Person_bool]
  }

  Y0_Sim <- matrix(NA_real_, nrow = N, ncol = Sims)
  Y1_Sim <- matrix(NA_real_, nrow = N, ncol = Sims)
  Y1_Sim[treat_Person_bool, 1] <- Y1[treat_Person_bool]
  Y0_Sim[!treat_Person_bool, 1] <- Y0[!treat_Person_bool]
  for (s in seq_len(Sims)) {
    Y1_Sim[treat_Person_bool, s] <- Y1[treat_Person_bool]
    Y0_Sim[!treat_Person_bool, s] <- Y0[!treat_Person_bool]
  }

  muY_Sim <- matrix(NA_real_, nrow = K, ncol = Sims)
  PhiY_Sim <- matrix(NA_real_, nrow = J, ncol = Sims)
  PhiY_Sim[, 1] <- 0
  Sigma2_Sim <- rep(NA_real_, Sims)
  Sigma2_Sim[1] <- 1
  Tau2_PhiY_Sim <- rep(NA_real_, Sims)
  Tau2_PhiY_Sim[1] <- 1
  PhiY_Person <- rep(0, N)

  Tau2_PhiD_Sim <- rep(NA_real_, Sims)
  Tau2_PhiD_Sim[1] <- 0.1
  PhiD_Person <- rep(0, N)
  PhiD_Sim <- matrix(0, J, Sims)

  Alpha_Sim <- vector("list", K)
  for (k in seq_len(K)) {
    Alpha_Sim[[k]] <- matrix(NA_real_, nrow = ncol(X_D), ncol = Sims)
  }

  strata_design <- sapply(seq_len(K), function(k) as.numeric(Strata_Person == k))
  if (K == 1L) {
    strata_design <- matrix(strata_design, ncol = 1L)
  }
  colnames(strata_design) <- paste0("strata", seq_len(K))
  CheckData_comp <- as.data.frame(cbind(D = D_Sim[, 1], strata_design, X_shared))
  CheckData_comp$D[!treat_Person_bool] <- NA
  prob_D <- rep(NA_real_, length(D))
  reg_comp <- stats::glm(D ~ . - 1, data = CheckData_comp, family = stats::binomial("probit"))
  comp_coef <- summary(reg_comp)$coefficients[, 1]
  muD_init <- comp_coef[1:K]
  alpha_shared_init <- comp_coef[-(1:K)]
  for (k in seq_len(K)) {
    Alpha_Sim[[k]][, 1] <- c(muD_init[k], alpha_shared_init)
    people_in_k <- Strata_Person == k
    prob_D[people_in_k] <- stats::pnorm(drop(X_D[people_in_k, , drop = FALSE] %*% Alpha_Sim[[k]][, 1]))
  }
  D_Sim[is.na(OutcomeData$D), 1] <- stats::rbinom(sum(is.na(OutcomeData$D)), 1, prob_D[is.na(OutcomeData$D)])

  U_Sim <- matrix(NA_real_, nrow = nrow(X_D), ncol = Sims)
  U_Sim[D_Sim[, 1] == 1, 1] <- truncnorm::rtruncnorm(sum(D_Sim[, 1]), a = 0, b = Inf)
  U_Sim[D_Sim[, 1] == 0, 1] <- truncnorm::rtruncnorm(sum(1 - D_Sim[, 1]), a = -Inf, b = 0)

  X_Main_obs <- cbind(1, as.matrix(X))
  X_Main_obs <- cbind(X_Main_obs, (1 - OutcomeData$Treat) * D_Sim[, 1] * 1)
  X_Main_obs <- cbind(X_Main_obs, OutcomeData$Treat * D_Sim[, 1] * cbind(as.matrix(X), 1))
  colnames(X_Main_obs) <- paste0("V", seq_len(ncol(X_Main_obs)))

  X_Main_treat <- cbind(1, as.matrix(X))
  X_Main_treat <- cbind(X_Main_treat, 0 * D_Sim[, 1] * 1)
  X_Main_treat <- cbind(X_Main_treat, 1 * D_Sim[, 1] * cbind(as.matrix(X), 1))
  colnames(X_Main_treat) <- paste0("V", seq_len(ncol(X_Main_treat)))

  X_Main_control <- cbind(1, as.matrix(X))
  X_Main_control <- cbind(X_Main_control, 1 * D_Sim[, 1] * 1)
  X_Main_control <- cbind(X_Main_control, 0 * D_Sim[, 1] * cbind(as.matrix(X), 1))
  colnames(X_Main_control) <- paste0("V", seq_len(ncol(X_Main_control)))

  YCN <- ncol(X)
  BetaD0_Sim <- vector("list", K)
  BetaD1_Sim <- vector("list", K)
  Delta0_Sim <- matrix(NA_real_, nrow = K, ncol = Sims)
  Delta1_Sim <- matrix(NA_real_, nrow = K, ncol = Sims)
  Beta_full_Sim <- vector("list", K)

  for (k in seq_len(K)) {
    BetaD0_Sim[[k]] <- matrix(NA_real_, nrow = YCN, ncol = Sims)
    BetaD1_Sim[[k]] <- matrix(NA_real_, nrow = YCN, ncol = Sims)
    Beta_full_Sim[[k]] <- matrix(NA_real_, nrow = ncol(X_Main_obs), ncol = Sims)
  }

  Y_obs <- matrix(NA_real_, nrow = nrow(OutcomeData), ncol = Sims)
  for (s in seq_len(Sims)) {
    Y_obs[treat_Person_bool, s] <- Y1[treat_Person_bool]
    Y_obs[!treat_Person_bool, s] <- Y0[!treat_Person_bool]
  }

  CheckData <- as.data.frame(cbind(Y_obs[, 1], X_Main_obs))
  colnames(CheckData)[1] <- "Y"

  for (k in seq_len(K)) {
    reg1k <- stats::lm(Y ~ . - 1, data = CheckData[Strata_Person == k, ])
    muY_Sim[k, 1] <- summary(reg1k)$coefficients[1, 1]
    BetaD0_Sim[[k]][, 1] <- summary(reg1k)$coefficients[1 + seq_len(nrow(BetaD0_Sim[[1]])), 1]
    Delta0_Sim[k, 1] <- summary(reg1k)$coefficients[2 + nrow(BetaD0_Sim[[1]]), 1]
    BetaD1_Sim[[k]][, 1] <- summary(reg1k)$coefficients[2 + nrow(BetaD0_Sim[[1]]) + seq_len(nrow(BetaD1_Sim[[1]])), 1]
    Delta1_Sim[k, 1] <- summary(reg1k)$coefficients[ncol(X_Main_obs), 1]

    Y_obs[treat_Person_bool & is.na(Y1), 1] <- stats::predict(reg1k, newdata = as.data.frame(X_Main_treat[treat_Person_bool & is.na(Y1), ]))
    Y_obs[!treat_Person_bool & is.na(Y0), 1] <- stats::predict(reg1k, newdata = as.data.frame(X_Main_control[!treat_Person_bool & is.na(Y0), ]))
    Y1_Sim[Strata_Person == k & !treat_Person_bool, 1] <- stats::predict(reg1k, newdata = as.data.frame(X_Main_treat[Strata_Person == k & !treat_Person_bool, ]))
    Y0_Sim[Strata_Person == k & treat_Person_bool, 1] <- stats::predict(reg1k, newdata = as.data.frame(X_Main_control[Strata_Person == k & treat_Person_bool, ]))
    Beta_full_Sim[[k]][, 1] <- c(muY_Sim[k, 1], BetaD0_Sim[[k]][, 1], Delta0_Sim[k, 1], BetaD1_Sim[[k]][, 1], Delta1_Sim[k, 1])
  }

  Y1_Sim[treat_Person_bool, 1] <- Y_obs[treat_Person_bool, 1]
  Y0_Sim[!treat_Person_bool, 1] <- Y_obs[!treat_Person_bool, 1]

  covsX <- c(2:(YCN + 1), (YCN + 3):(ncol(X_Main_obs) - 1))
  noncovsX <- c(1, YCN + 2, ncol(X_Main_obs))
  nu_0 <- 7
  lambda_pk <- 3
  V_0_inv <- diag(c(1 / 1000, 1 / 100000))
  m_0 <- matrix(0, nrow = L + P, ncol = 1)
  S_0 <- diag(1 / 10000, nrow = L + P)
  m_alpha <- rep(0, K + ncol(X_shared))
  v_alpha_inv <- diag(1 / 25, K + ncol(X_shared))
  m_beta <- rep(0, length(covsX))
  m_beta_k <- rep(0, length(noncovsX))
  v_beta_inv <- diag(1 / 100, length(covsX))
  v_beta_k_inv <- diag(c(1 / 10000, 1 / 16, 1 / 100))
  a_sigmas <- 1
  b_sigmas <- 1
  max_tauPhiD <- 10
  max_tauPhiY <- 225

  estimand_names <- make_primary_estimand_names(K)
  outcome_level_names <- make_outcome_level_names(K)
  estimands_raw <- matrix(NA_real_, nrow = n_imps, ncol = length(estimand_names), dimnames = list(NULL, estimand_names))
  estimands_relabelled <- matrix(NA_real_, nrow = n_imps, ncol = length(estimand_names), dimnames = list(NULL, estimand_names))
  outcome_levels_raw <- matrix(NA_real_, nrow = n_imps, ncol = length(outcome_level_names), dimnames = list(NULL, outcome_level_names))
  outcome_levels_relabelled <- matrix(NA_real_, nrow = n_imps, ncol = length(outcome_level_names), dimnames = list(NULL, outcome_level_names))
  mu_use <- array(NA_real_, dim = c(L + P, K, n_imps))
  Sigma_use <- vector("list", n_imps)
  pi_use <- matrix(NA_real_, nrow = K, ncol = n_imps)
  ordered_index_use <- matrix(NA_integer_, nrow = n_imps, ncol = K)
  raw_strata_use <- matrix(NA_integer_, nrow = J, ncol = n_imps)
  D_use <- matrix(NA_real_, nrow = N, ncol = n_imps)
  Y_obs_use <- matrix(NA_real_, nrow = N, ncol = n_imps)
  PPcheckD <- matrix(NA_real_, nrow = N, ncol = n_imps)
  PPcheckY <- matrix(NA_real_, nrow = N, ncol = n_imps)
  PPcheckS <- matrix(NA_real_, nrow = J, ncol = n_imps)
  PPcheckS_person <- matrix(NA_real_, nrow = N, ncol = n_imps)

  missingYs <- (treat_Person_bool & is.na(OutcomeData$Y1)) | (!treat_Person_bool & is.na(OutcomeData$Y0))
  Nj <- table(person_Fac)
  CZ_i <- CZ_Matrix
  C_i <- C_Matrix
  Z_i <- Z_Matrix

  X_Main_comp <- rebuild_outcome_design_original(X, treat_Person_bin, 1)
  X_Main_nocomp <- rebuild_outcome_design_original(X, treat_Person_bin, 0)

  for (i in seq_len(Sims - 1)) {
    new <- i + 1

    for (l in seq_len(L)) {
      CZ_i[treat_Fac_bin == 0, l] <- C_Sim[[l]][treat_Fac_bin == 0, i]
      C_i[treat_Fac_bin == 0, l] <- C_Sim[[l]][treat_Fac_bin == 0, i]
    }

    new_mu <- update_mixture_means_original(
      CZ_i = CZ_i,
      strata_fac = S_Sim[, i],
      Sigma_prev = Sigma_Sim[[i]],
      V_0_inv = V_0_inv,
      m_0 = m_0,
      K = K
    )
    for (k in seq_len(K)) {
      mu_Sim[[k]][, new] <- new_mu[, k]
    }

    Sigma_Sim[[new]] <- update_mixture_covariance_original(
      C_i = C_i,
      Z_i = Z_i,
      strata_fac = S_Sim[, i],
      mu_draws = new_mu,
      S_0 = S_0,
      nu_0 = nu_0,
      J = J,
      K = K
    )

    pi_Sim[, new] <- as.numeric(update_mixture_weights_original(
      strata_fac = Strata_Fac,
      lambda_pk = lambda_pk,
      K = K
    ))

    U_Sim[, new] <- sample_latent_u_original(
      D_current = D_Sim[, i],
      X_D = X_D,
      alpha_list_prev = lapply(seq_len(K), function(k) Alpha_Sim[[k]][, i]),
      phiD_person = PhiD_Person,
      strata_person = Strata_Person,
      K = K
    )

    PhiD_Sim[, new] <- update_compliance_random_effects_original(
      U_new = U_Sim[, new],
      X_D = X_D,
      alpha_list_prev = lapply(seq_len(K), function(k) Alpha_Sim[[k]][, i]),
      strata_person = Strata_Person,
      strata_fac = Strata_Fac,
      provider_fac = provider_Fac,
      person_fac = person_Fac,
      tau2_prev = Tau2_PhiD_Sim[i],
      Nj = Nj,
      K = K
    )
    PhiD_Person <- PhiD_Sim[person_Fac, new]

    Tau2_PhiD_Sim[new] <- update_truncated_variance_original(
      PhiD_Sim[, new], max_tauPhiD,
      context = sprintf("chain %d, Gibbs iteration %d, compliance", chain_id, new)
    )

    alpha_update <- update_compliance_coefficients_step1(
      U_new = U_Sim[, new],
      X_D = X_D,
      strata_person = Strata_Person,
      phiD_person = PhiD_Person,
      v_alpha_inv = v_alpha_inv,
      m_alpha = m_alpha,
      K = K
    )
    for (k in seq_len(K)) {
      Alpha_Sim[[k]][, new] <- alpha_update$alpha_list_draw[[k]]
    }

    PhiY_Sim[, new] <- update_outcome_random_effects_original(
      Y_obs = Y_obs[, i],
      X_Main_obs = X_Main_obs,
      beta_full_prev = lapply(seq_len(K), function(k) Beta_full_Sim[[k]][, i]),
      strata_person = Strata_Person,
      strata_fac = Strata_Fac,
      provider_fac = provider_Fac,
      person_fac = person_Fac,
      sigma2_prev = Sigma2_Sim[i],
      tau2_prev = Tau2_PhiY_Sim[i],
      Nj = Nj,
      K = K
    )
    PhiY_Person <- PhiY_Sim[person_Fac, new]

    Tau2_PhiY_Sim[new] <- update_truncated_variance_original(
      PhiY_Sim[, new], max_tauPhiY,
      context = sprintf("chain %d, Gibbs iteration %d, outcome", chain_id, new)
    )

    outcome_noncov_update <- update_outcome_noncovariates_original(
      Y_obs = Y_obs[, i],
      X_Main_obs = X_Main_obs,
      phiY_person = PhiY_Person,
      strata_person = Strata_Person,
      beta_full_prev = lapply(seq_len(K), function(k) Beta_full_Sim[[k]][, i]),
      sigma2_prev = Sigma2_Sim[i],
      noncovsX = noncovsX,
      covsX = covsX,
      v_beta_k_inv = v_beta_k_inv,
      m_beta_k = m_beta_k,
      K = K
    )
    muY_Sim[, new] <- outcome_noncov_update$muY_draw
    Delta0_Sim[, new] <- outcome_noncov_update$delta0_draw
    Delta1_Sim[, new] <- outcome_noncov_update$delta1_draw

    Beta_covs_new <- update_outcome_shared_covariates_original(
      Y_obs = Y_obs[, i],
      X_Main_obs = X_Main_obs,
      phiY_person = PhiY_Person,
      strata_person = Strata_Person,
      muY_draw = muY_Sim[, new],
      delta0_draw = Delta0_Sim[, new],
      delta1_draw = Delta1_Sim[, new],
      sigma2_prev = Sigma2_Sim[i],
      covsX = covsX,
      noncovsX = noncovsX,
      v_beta_inv = v_beta_inv,
      m_beta = m_beta
    )
    beta_full_draw <- assemble_outcome_betas_original(
      muY_draw = muY_Sim[, new],
      beta_covs_draw = Beta_covs_new,
      delta0_draw = Delta0_Sim[, new],
      delta1_draw = Delta1_Sim[, new],
      YCN = YCN,
      K = K
    )
    for (k in seq_len(K)) {
      Beta_full_Sim[[k]][, new] <- beta_full_draw[[k]]
    }

    Sigma2_Sim[new] <- update_sigma2_original(
      Y_obs = Y_obs[, i],
      X_Main_obs = X_Main_obs,
      beta_full_draw = beta_full_draw,
      strata_person = Strata_Person,
      phiY_person = PhiY_Person,
      a_sigmas = a_sigmas,
      b_sigmas = b_sigmas,
      N = N,
      K = K
    )

    strata_update <- update_strata_membership_original(
      CZ_i = CZ_i,
      Y_obs = Y_obs[, i],
      X_Main_obs = X_Main_obs,
      D_current = D_Sim[, i],
      X_D = X_D,
      mu_draws = new_mu,
      Sigma_new = Sigma_Sim[[new]],
      beta_full_draw = beta_full_draw,
      Sigma2_new = Sigma2_Sim[new],
      alpha_list_draw = lapply(seq_len(K), function(k) Alpha_Sim[[k]][, new]),
      phiY_person = PhiY_Person,
      phiD_person = PhiD_Person,
      pi_new = pi_Sim[, new],
      person_fac = person_Fac,
      J = J,
      K = K,
      facility_keys = inputs$facility_keys,
      context = sprintf("chain %d, Gibbs iteration %d", chain_id, new)
    )
    S_Sim[, new] <- strata_update$strata_fac_draw
    Strata_Fac <- S_Sim[, new]
    Strata_Person <- Strata_Fac[Facility]

    C_draw_new <- impute_control_c_original(
      C_Sim_prev = lapply(seq_len(L), function(l) C_Sim[[l]][, i]),
      Z_Matrix = Z_Matrix,
      mu_draws = new_mu,
      Sigma_new = Sigma_Sim[[new]],
      strata_fac_new = S_Sim[, new],
      treat_fac_bool = treat_Fac_bool,
      control_index = control_index,
      L = L,
      P = P,
      K = K
    )
    for (l in seq_len(L)) {
      C_Sim[[l]][, new] <- C_draw_new[[l]]
    }

    D_update <- impute_unobserved_d_original(
      Y_obs = Y_obs[, i],
      X_Main_comp = X_Main_comp,
      X_Main_nocomp = X_Main_nocomp,
      beta_full_draw = beta_full_draw,
      phiY_person = PhiY_Person,
      phiD_person = PhiD_Person,
      X_D = X_D,
      alpha_list_draw = lapply(seq_len(K), function(k) Alpha_Sim[[k]][, new]),
      strata_person = Strata_Person,
      D_missing = is.na(OutcomeData$D),
      Sigma2_new = Sigma2_Sim[new]
    )
    D_Sim[!is.na(D_update$D_draw), new] <- D_update$D_draw[!is.na(D_update$D_draw)]

    X_Main_obs <- rebuild_outcome_design_original(
      X = X,
      treat_person_bin = OutcomeData$Treat,
      D_draw = D_Sim[, new]
    )

    Y_obs[, new] <- Y_obs[, i]
    Y_obs[, new] <- impute_missing_y_original(
      Y_obs_prev = Y_obs[, new],
      missingYs = missingYs,
      X_Main_obs = X_Main_obs,
      beta_full_draw = beta_full_draw,
      phiY_person = PhiY_Person,
      strata_person = Strata_Person,
      Sigma2_new = Sigma2_Sim[new],
      K = K
    )

    Y1_Sim[treat_Person_bool, new] <- Y_obs[treat_Person_bool, new]
    Y0_Sim[!treat_Person_bool, new] <- Y_obs[!treat_Person_bool, new]

    if (new %in% samps_use) {
      thin_ind <- match(new, samps_use)
      mu_use[, , thin_ind] <- new_mu
      Sigma_use[[thin_ind]] <- Sigma_Sim[[new]]
      pi_use[, thin_ind] <- pi_Sim[, new]
      raw_strata_use[, thin_ind] <- S_Sim[, new]
      D_use[, thin_ind] <- D_Sim[, new]
      Y_obs_use[, thin_ind] <- Y_obs[, new]

      Coeff_Alpha_est <- matrix(NA_real_, nrow = nrow(Alpha_Sim[[1]]), ncol = K)
      Coeff_BetaD1_est <- matrix(NA_real_, nrow = YCN, ncol = K)
      Coeff_BetaD0_est <- matrix(NA_real_, nrow = YCN, ncol = K)
      Coeff_D1_est <- rep(NA_real_, K)
      Coeff_D0_est <- rep(NA_real_, K)
      muY_est <- rep(NA_real_, K)
      Tau2_PhiD_est <- rep(Tau2_PhiD_Sim[new], K)
      pi_est <- pi_Sim[, new]
      X_k <- vector("list", K)
      XD_k <- vector("list", K)

      for (k in seq_len(K)) {
        Coeff_Alpha_est[, k] <- Alpha_Sim[[k]][, new]
        muY_est[k] <- Beta_full_Sim[[k]][1, new]
        Coeff_BetaD1_est[, k] <- Beta_full_Sim[[k]][(YCN + 3):(2 * YCN + 2), new]
        Coeff_BetaD0_est[, k] <- Beta_full_Sim[[k]][2:(YCN + 1), new]
        Coeff_D1_est[k] <- Delta1_Sim[k, new]
        Coeff_D0_est[k] <- Delta0_Sim[k, new]
        X_k[[k]] <- X[Strata_Person == k, , drop = FALSE]
        XD_k[[k]] <- X_D[Strata_Person == k, , drop = FALSE]
      }

      Y1_k <- rep(NA_real_, K)
      Y0_k <- rep(NA_real_, K)
      ITT_k <- rep(NA_real_, K)
      for (k in seq_len(K)) {
        XD_Alpha <- as.matrix(XD_k[[k]]) %*% Coeff_Alpha_est[, k]
        p_ijk <- stats::pnorm(XD_Alpha / sqrt(Tau2_PhiD_est[k] + 1))
        muY_weighted_XBeta0_k <- muY_est[k] + mean(X_k[[k]] %*% Coeff_BetaD0_est[, k])
        weighted_XBeta1_delta1_k <- mean(p_ijk * (X_k[[k]] %*% Coeff_BetaD1_est[, k] + Coeff_D1_est[k]))
        weighted_delta0_k <- mean(p_ijk * Coeff_D0_est[k])
        Y1_k[k] <- muY_weighted_XBeta0_k + weighted_XBeta1_delta1_k
        Y0_k[k] <- muY_weighted_XBeta0_k + weighted_delta0_k
        ITT_k[k] <- Y1_k[k] - Y0_k[k]
      }

      p_ijk_mean <- rep(NA_real_, K)
      Y1_D_k <- rep(NA_real_, K)
      Y0_D_k <- rep(NA_real_, K)
      CACE_k <- rep(NA_real_, K)
      for (k in seq_len(K)) {
        XD_Alpha <- as.matrix(XD_k[[k]]) %*% Coeff_Alpha_est[, k]
        p_ijk <- stats::pnorm(XD_Alpha / sqrt(Tau2_PhiD_est[k] + 1))
        p_ijk_mean[k] <- mean(p_ijk)
        weighted_XBeta0_kBeta1_k <- mean(p_ijk * X_k[[k]] %*% (Coeff_BetaD0_est[, k] + Coeff_BetaD1_est[, k]))
        weighted_XBeta0_k <- mean(p_ijk * X_k[[k]] %*% Coeff_BetaD0_est[, k])
        Y1_D_k[k] <- muY_est[k] + weighted_XBeta0_kBeta1_k / p_ijk_mean[k] + Coeff_D1_est[k]
        Y0_D_k[k] <- muY_est[k] + weighted_XBeta0_k / p_ijk_mean[k] + Coeff_D0_est[k]
        CACE_k[k] <- Y1_D_k[k] - Y0_D_k[k]
      }

      ITT <- sum(pi_est * ITT_k)
      ITT_Y1 <- sum(pi_est * Y1_k)
      ITT_Y0 <- sum(pi_est * Y0_k)
      CACE_num <- sum(pi_est * CACE_k * p_ijk_mean)
      CACE_den <- sum(pi_est * p_ijk_mean)
      CACE <- CACE_num / CACE_den
      CACE_Y1 <- sum(pi_est * p_ijk_mean * Y1_D_k) / CACE_den
      CACE_Y0 <- sum(pi_est * p_ijk_mean * Y0_D_k) / CACE_den

      ordered_index <- identify_ordered_strata(new_mu)
      ordered_index_use[thin_ind, ] <- ordered_index
      ITT_k_ordered <- ITT_k[ordered_index]
      CACE_k_ordered <- CACE_k[ordered_index]
      Y1_k_ordered <- Y1_k[ordered_index]
      Y0_k_ordered <- Y0_k[ordered_index]
      Y1_D_k_ordered <- Y1_D_k[ordered_index]
      Y0_D_k_ordered <- Y0_D_k[ordered_index]

      estimands_raw[thin_ind, ] <- assemble_primary_estimand_draw(ITT, ITT_k, CACE, CACE_k)
      outcome_levels_raw[thin_ind, ] <- assemble_outcome_level_draw(
        ITT_Y1, ITT_Y0, Y1_k, Y0_k, CACE_Y1, CACE_Y0, Y1_D_k, Y0_D_k
      )

      estimands_relabelled[thin_ind, ] <- assemble_primary_estimand_draw(ITT, ITT_k_ordered, CACE, CACE_k_ordered)
      outcome_levels_relabelled[thin_ind, ] <- assemble_outcome_level_draw(
        ITT_Y1, ITT_Y0, Y1_k_ordered, Y0_k_ordered, CACE_Y1, CACE_Y0, Y1_D_k_ordered, Y0_D_k_ordered
      )

      if (K == 2L) {
        PPcheckS[, thin_ind] <- stats::rbinom(J, 1, pi_Sim[2, new]) + 1
      } else {
        PPcheckS[, thin_ind] <- sample(seq_len(K), J, replace = TRUE, prob = pi_Sim[, new])
      }
      PPcheckS_person[, thin_ind] <- PPcheckS[, thin_ind][Facility]

      for (k in seq_len(K)) {
        people_in_k <- PPcheckS_person[, thin_ind] == k
        PP_phiD <- stats::rnorm(J, mean = 0, sd = sqrt(Tau2_PhiD_Sim[new]))
        PP_phiY <- stats::rnorm(J, mean = 0, sd = sqrt(Tau2_PhiY_Sim[new]))
        PP_phiD_person <- PP_phiD[Facility]
        PP_phiY_person <- PP_phiY[Facility]

        PP_Dmean <- stats::pnorm(X_D[people_in_k, , drop = FALSE] %*% Alpha_Sim[[k]][, new] + PP_phiD_person[people_in_k])
        PPcheckD[people_in_k, thin_ind] <- stats::rbinom(sum(people_in_k), 1, PP_Dmean)

        PP_X_Main_obsk <- cbind(1, as.matrix(X[people_in_k, , drop = FALSE]))
        PP_X_Main_obsk <- cbind(PP_X_Main_obsk, (1 - treat_Person_bin[people_in_k]) * PPcheckD[people_in_k, thin_ind] * 1)
        PP_X_Main_obsk <- cbind(PP_X_Main_obsk, treat_Person_bin[people_in_k] * PPcheckD[people_in_k, thin_ind] * cbind(as.matrix(X[people_in_k, , drop = FALSE]), 1))
        colnames(PP_X_Main_obsk) <- paste0("V", seq_len(ncol(PP_X_Main_obsk)))

        PP_Ymean <- PP_X_Main_obsk %*% Beta_full_Sim[[k]][, new] + PP_phiY_person[people_in_k]
        PPcheckY[people_in_k, thin_ind] <- stats::rnorm(sum(people_in_k), PP_Ymean, sqrt(Sigma2_Sim[new]))
      }
    }
  }

  if (isTRUE(make_trace_pdf)) {
    grDevices::pdf(file.path(chain_output_dir, "Continuous_Outcome_single_chain_traces.pdf"))
    on.exit(grDevices::dev.off(), add = TRUE)
    for (estimand_name in colnames(estimands_relabelled)) {
      graphics::plot(seq_len(n_imps), estimands_relabelled[, estimand_name], type = "l", xlab = "Retained draw", ylab = estimand_name, main = estimand_name)
    }
  }

  list(
    config = config,
    chain_id = as.integer(chain_id),
    posterior_draws = list(
      estimands_raw = estimands_raw,
      estimands_relabelled = estimands_relabelled,
      outcome_levels_raw = outcome_levels_raw,
      outcome_levels_relabelled = outcome_levels_relabelled,
      mu_use = mu_use,
      Sigma_use = Sigma_use,
      pi_use = pi_use,
      ordered_index = ordered_index_use,
      S_use = raw_strata_use,
      D_use = D_use,
      Y_obs_use = Y_obs_use,
      PPcheckD = PPcheckD,
      PPcheckY = PPcheckY,
      PPcheckS = PPcheckS,
      PPcheckS_person = PPcheckS_person
    ),
    posterior_summary_raw = summarize_primary_estimands(as.data.frame(estimands_raw)),
    posterior_summary_relabelled = summarize_primary_estimands(as.data.frame(estimands_relabelled)),
    outcome_level_summary_raw = summarize_outcome_levels(as.data.frame(outcome_levels_raw)),
    outcome_level_summary_relabelled = summarize_outcome_levels(as.data.frame(outcome_levels_relabelled)),
    pp_check_results = compute_combined_ppc_results(list(list(
      posterior_draws = list(
        ordered_index = ordered_index_use,
        mu_use = mu_use,
        Sigma_use = Sigma_use,
        D_use = D_use,
        Y_obs_use = Y_obs_use,
        S_use = raw_strata_use,
        PPcheckD = PPcheckD,
        PPcheckY = PPcheckY,
        PPcheckS_person = PPcheckS_person
      ),
      fixed_data = list(
        treat_Person_bool = treat_Person_bool,
        Facility = Facility,
        CZ_Matrix = CZ_Matrix
      )
    )))$pp_check_results,
    final_state = list(
      S_last = S_Sim[, Sims],
      pi_last = pi_Sim[, Sims],
      mu_last = do.call(cbind, lapply(seq_len(K), function(k) mu_Sim[[k]][, Sims])),
      Sigma_last = Sigma_Sim[[Sims]]
    ),
    fixed_data = list(
      treat_Person_bool = treat_Person_bool,
      Facility = Facility,
      facility_keys = inputs$facility_keys,
      treat_Fac_bool = treat_Fac_bool,
      CZ_Matrix = CZ_Matrix
    )
  )
}

# Multichain analysis function.
#
# Inputs:
# - `config`: common analysis configuration shared across chains.
# - `n_chains`: number of independent chains.
# - `chain_seeds`: reproducible chain-specific random seeds.
# - `make_trace_pdf`: whether to save combined trace plots for the primary
#   estimands.
#
# Output:
# - A list containing the chain-specific fits, pooled estimand draws, pooled
#   posterior summaries, MCMC diagnostics, and combined-chain posterior
#   predictive checks.
#
# Statistical role:
# This function is the main inferential wrapper for the application analysis.
# It combines independent chains to support convergence assessment, chain
# pooling, and reporting of posterior predictive summaries on the common
# relabeled estimand scale.
run_continuous_outcome_multichain <- function(
    config = make_continuous_outcome_config(),
    n_chains = 4L,
    chain_seeds = config$seed + seq_len(n_chains) - 1L,
    make_trace_pdf = TRUE) {
  n_chains <- as.integer(n_chains)
  if (length(chain_seeds) != n_chains) {
    stop("`chain_seeds` must have length equal to `n_chains`.", call. = FALSE)
  }

  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

  chain_results <- vector("list", n_chains)
  for (chain_id in seq_len(n_chains)) {
    chain_config <- config
    chain_config$seed <- as.integer(chain_seeds[chain_id])
    chain_results[[chain_id]] <- run_continuous_outcome_single_chain(
      config = chain_config,
      chain_id = chain_id,
      make_trace_pdf = FALSE
    )
  }
  names(chain_results) <- paste0("chain_", seq_len(n_chains))

  pooled_draws <- combine_primary_estimand_draws(chain_results)
  pooled_summary <- summarize_primary_estimands(as.data.frame(pooled_draws))
  pooled_outcome_level_draws <- combine_outcome_level_draws(chain_results)
  pooled_outcome_level_summary <- summarize_outcome_levels(as.data.frame(pooled_outcome_level_draws))
  diagnostics <- compute_primary_estimand_diagnostics(chain_results)
  latent_class_diagnostics <- compute_latent_class_diagnostics(chain_results)
  relabeling_stability <- compute_relabeling_stability(chain_results)
  ppc_results <- compute_combined_ppc_results(chain_results)
  expanded_ppc_results <- compute_expanded_ppc_results(chain_results, outcome_label = "Continuous outcome")

  if (isTRUE(make_trace_pdf)) {
    make_chain_trace_pdf(
      chain_results = chain_results,
      file = file.path(config$output_dir, "Continuous_Outcome_multichain_traces.pdf")
    )
  }

  list(
    config = config,
    chain_seeds = as.integer(chain_seeds),
    chain_results = chain_results,
    pooled_estimand_draws = pooled_draws,
    pooled_posterior_summary = pooled_summary,
    pooled_outcome_level_draws = pooled_outcome_level_draws,
    pooled_outcome_level_summary = pooled_outcome_level_summary,
    diagnostics = diagnostics,
    latent_class_diagnostics = latent_class_diagnostics,
    relabeling_stability = relabeling_stability,
    pp_check_results = ppc_results$pp_check_results,
    pp_check_draws = ppc_results$pp_check_draws,
    pp_check_observed_means = ppc_results$pp_check_observed_means,
    expanded_ppc_results = expanded_ppc_results$summary,
    expanded_ppc_implementation_draws = expanded_ppc_results$implementation_draws
  )
}

export_continuous_outcome_tables <- function(fit, output_dir = fit$config$output_dir, prefix = "Continuous_Outcome") {
  table_dir <- file.path(output_dir, "reviewer_tables")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  write_table <- function(object, suffix) {
    utils::write.csv(
      object,
      file = file.path(table_dir, paste0(prefix, "_", suffix, ".csv")),
      row.names = FALSE
    )
  }

  write_table(fit$pooled_posterior_summary, "posterior_summary")
  write_table(fit$pooled_outcome_level_summary, "outcome_level_summary")
  write_table(fit$diagnostics, "mcmc_diagnostics")
  write_table(fit$latent_class_diagnostics$facility_membership, "facility_membership")
  write_table(fit$latent_class_diagnostics$ambiguous_membership, "ambiguous_membership")
  write_table(fit$latent_class_diagnostics$class_size_summary, "class_size_summary")
  write_table(fit$relabeling_stability, "relabeling_stability")
  write_table(fit$pp_check_results, "original_ppc")
  write_table(fit$expanded_ppc_results, "expanded_ppc")

  run_settings <- data.frame(
    setting = c("seed", "K", "burn", "thin", "n_imps", "n_chains", "chain_seeds"),
    value = c(
      fit$config$seed,
      fit$config$K,
      fit$config$burn,
      fit$config$thin,
      fit$config$n_imps,
      length(fit$chain_seeds),
      paste(fit$chain_seeds, collapse = ", ")
    ),
    row.names = NULL
  )
  write_table(run_settings, "run_settings")

  invisible(table_dir)
}
