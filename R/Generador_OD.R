# ============================================================
# Generador OD para RLM con predictores mixtos, dependencia
# Normalidad y homocedasticidad en errores
# No linealidad leve (cuadrática), con interacciones
# ============================================================

gen_od_rlm <- function(
  n,
  n_cont_norm,
  n_cont_skew,
  n_bin,
  n_fact,
  fact_levels,
  rho,
  n_latent,
  latent_strength,
  R2_target,
  nonlinear_frac,
  nonlinear_strength,
  beta_sd,
  seed = NULL,
  include_interactions = FALSE,
  interaction_types = c("xN:xF"),
  interaction_frac = 0.3,
  interaction_sd = NULL,
  fixed_interactions = NULL,
  fit_with_interactions = TRUE
) {
  if (!is.null(seed)) set.seed(seed)

  # ---- Chequeos ----
  stopifnot(n > 10)
  stopifnot(n_cont_norm >= 0, n_cont_skew >= 0, n_bin >= 0, n_fact >= 0)
  if (n_fact > 0) stopifnot(length(fact_levels) == n_fact)

  if (is.null(interaction_sd)) {
    interaction_sd <- beta_sd / 2
  }

  # ---- 1) Latentes para dependencia global ----
  F_lat <- matrix(rnorm(n * n_latent), nrow = n, ncol = n_latent)
  colnames(F_lat) <- paste0("F", seq_len(n_latent))

  # ---- 2) Continuas normales correlacionadas ----
  X_cont_norm <- NULL
  if (n_cont_norm > 0) {
    idx <- seq_len(n_cont_norm)
    Sigma <- outer(idx, idx, function(i, j) rho^abs(i - j))

    Z <- MASS::mvrnorm(n = n, mu = rep(0, n_cont_norm), Sigma = Sigma)

    L_load <- matrix(rnorm(n_cont_norm * n_latent), nrow = n_cont_norm, ncol = n_latent)
    L_load <- latent_strength * L_load / sqrt(rowSums(L_load^2) + 1e-12)

    X_cont_norm <- Z + F_lat %*% t(L_load)
    colnames(X_cont_norm) <- paste0("xN", seq_len(n_cont_norm))
  }

  # ---- 3) Continuas sesgadas ----
  X_cont_skew <- NULL
  if (n_cont_skew > 0) {
    Zs <- matrix(rnorm(n * n_cont_skew), nrow = n, ncol = n_cont_skew)

    Ls <- matrix(rnorm(n_cont_skew * n_latent), nrow = n_cont_skew, ncol = n_latent)
    Ls <- latent_strength * Ls / sqrt(rowSums(Ls^2) + 1e-12)

    U <- Zs + F_lat %*% t(Ls)

    X_cont_skew <- exp(0.6 * U)
    X_cont_skew <- scale(X_cont_skew)
    colnames(X_cont_skew) <- paste0("xS", seq_len(n_cont_skew))
  }

  # ---- 4) Binarias ----
  X_bin <- NULL
  if (n_bin > 0) {
    Wb <- matrix(rnorm(n_bin * n_latent), nrow = n_bin, ncol = n_latent)
    Wb <- latent_strength * Wb / sqrt(rowSums(Wb^2) + 1e-12)

    linp <- F_lat %*% t(Wb)

    if (!is.null(X_cont_norm)) {
      take <- min(2, n_cont_norm)
      linp <- linp + 0.3 * X_cont_norm[, seq_len(take), drop = FALSE] %*%
        matrix(rnorm(take * n_bin), nrow = take)
    }

    P <- 1 / (1 + exp(-linp))
    X_bin <- matrix(
      rbinom(n * n_bin, size = 1, prob = as.vector(P)),
      nrow = n, ncol = n_bin
    )
    colnames(X_bin) <- paste0("xB", seq_len(n_bin))
  }

  # ---- 5) Factores ----
  X_fact <- list()
  if (n_fact > 0) {
    for (k in seq_len(n_fact)) {
      K <- fact_levels[k]

      Wf <- matrix(rnorm(K * n_latent), nrow = K, ncol = n_latent)
      Wf <- latent_strength * Wf / sqrt(rowSums(Wf^2) + 1e-12)

      logits <- F_lat %*% t(Wf)
      logits <- sweep(logits, 1, apply(logits, 1, max), FUN = "-")
      probs <- exp(logits)
      probs <- probs / rowSums(probs)

      lev <- paste0("L", seq_len(K))
      draw_one <- function(p) sample(lev, size = 1, prob = p)
      fk <- apply(probs, 1, draw_one)

      X_fact[[k]] <- factor(fk, levels = lev)
      names(X_fact)[k] <- paste0("xF", k)
    }
  }

  # ---- 6) Data frame de predictores ----
  X_df <- as.data.frame(matrix(nrow = n, ncol = 0))

  if (!is.null(X_cont_norm)) X_df <- cbind(X_df, as.data.frame(scale(X_cont_norm)))
  if (!is.null(X_cont_skew)) X_df <- cbind(X_df, as.data.frame(X_cont_skew))
  if (!is.null(X_bin))       X_df <- cbind(X_df, as.data.frame(X_bin))
  if (length(X_fact) > 0)    X_df <- cbind(X_df, as.data.frame(X_fact))

  # ============================================================
    # ============================================================
    # 7) Selección de interacciones a nivel de variable
    # ============================================================
    get_vars_by_prefix <- function(df_names, prefix) {
      grep(paste0("^", prefix), df_names, value = TRUE)
    }
    
    make_pairs <- function(vars1, vars2, same_group = FALSE) {
      if (length(vars1) == 0 || length(vars2) == 0) return(character(0))
    
      if (same_group) {
        cmb <- utils::combn(vars1, 2)
        apply(cmb, 2, function(z) paste(z, collapse = ":"))
      } else {
        as.vector(outer(vars1, vars2, function(a, b) paste(a, b, sep = ":")))
      }
    }
    
    df_names <- names(X_df)
    
    vars_xN <- get_vars_by_prefix(df_names, "xN")
    vars_xS <- get_vars_by_prefix(df_names, "xS")
    vars_xB <- get_vars_by_prefix(df_names, "xB")
    vars_xF <- get_vars_by_prefix(df_names, "xF")
    
    candidate_interactions <- character(0)
    
    if (include_interactions) {
      for (tp in interaction_types) {
        if (tp == "xN:xF") candidate_interactions <- c(candidate_interactions, make_pairs(vars_xN, vars_xF))
        if (tp == "xN:xB") candidate_interactions <- c(candidate_interactions, make_pairs(vars_xN, vars_xB))
        if (tp == "xN:xN") candidate_interactions <- c(candidate_interactions, make_pairs(vars_xN, vars_xN, same_group = TRUE))
        if (tp == "xS:xF") candidate_interactions <- c(candidate_interactions, make_pairs(vars_xS, vars_xF))
        if (tp == "xS:xB") candidate_interactions <- c(candidate_interactions, make_pairs(vars_xS, vars_xB))
        if (tp == "xS:xS") candidate_interactions <- c(candidate_interactions, make_pairs(vars_xS, vars_xS, same_group = TRUE))
        if (tp == "xN:xS") candidate_interactions <- c(candidate_interactions, make_pairs(vars_xN, vars_xS))
        if (tp == "xB:xF") candidate_interactions <- c(candidate_interactions, make_pairs(vars_xB, vars_xF))
        if (tp == "xB:xB") candidate_interactions <- c(candidate_interactions, make_pairs(vars_xB, vars_xB, same_group = TRUE))
      }
    
      candidate_interactions <- unique(candidate_interactions)
    }
    
    if (!include_interactions) {
      selected_interactions <- character(0)
    
    } else if (!is.null(fixed_interactions)) {
      invalid_int <- setdiff(fixed_interactions, candidate_interactions)
    
      if (length(invalid_int) > 0) {
        stop(
          "Estas interacciones fijas no son válidas para este escenario: ",
          paste(invalid_int, collapse = ", ")
        )
      }
    
      selected_interactions <- fixed_interactions
    
    } else {
      n_int <- if (length(candidate_interactions) > 0) {
        max(1, floor(length(candidate_interactions) * interaction_frac))
      } else {
        0
      }
    
      selected_interactions <- if (n_int > 0) {
        sample(candidate_interactions,
               size = min(n_int, length(candidate_interactions)),
               replace = FALSE)
      } else {
        character(0)
      }
    }

  # ============================================================
  # 8) Matriz de diseño verdadera
  # ============================================================
  rhs_true <- "."
  if (length(selected_interactions) > 0) {
    rhs_true <- paste(rhs_true, "+", paste(selected_interactions, collapse = " + "))
  }

  formula_true <- as.formula(paste("~", rhs_true))
  X_mm <- model.matrix(formula_true, data = X_df)

  p_mm <- ncol(X_mm)
  beta_true <- rep(0, p_mm)
  names(beta_true) <- colnames(X_mm)
  beta_true["(Intercept)"] <- 0

  is_interaction_col <- grepl(":", colnames(X_mm), fixed = TRUE)
  is_main_col <- !is_interaction_col
  is_main_col[1] <- FALSE

  if (sum(is_main_col) > 0) {
    beta_true[is_main_col] <- rnorm(sum(is_main_col), mean = 0, sd = beta_sd)
  }
  if (sum(is_interaction_col) > 0) {
    beta_true[is_interaction_col] <- rnorm(sum(is_interaction_col), mean = 0, sd = interaction_sd)
  }

  # ============================================================
  # 9) No linealidad leve
  # ============================================================
  cont_names <- c(
    if (!is.null(X_cont_norm)) colnames(X_cont_norm) else character(0),
    if (!is.null(X_cont_skew)) colnames(X_cont_skew) else character(0)
  )
  cont_names <- intersect(cont_names, names(X_df))

  n_nl <- floor(length(cont_names) * nonlinear_frac)
  nl_vars <- if (n_nl > 0) sample(cont_names, n_nl) else character(0)

  gamma <- setNames(rep(0, length(nl_vars)), nl_vars)
  if (length(nl_vars) > 0) {
    gamma[] <- rnorm(length(nl_vars), mean = 0, sd = nonlinear_strength)
  }

  mu_lin <- as.vector(X_mm %*% beta_true)
  mu_nl <- mu_lin

  if (length(nl_vars) > 0) {
    for (v in nl_vars) {
      z <- as.numeric(scale(X_df[[v]]))
      mu_nl <- mu_nl + gamma[v] * (z^2)
    }
  }

  # ============================================================
  # 10) Ajuste de sigma para R2 objetivo
  # ============================================================
  signal_var <- var(mu_nl)
  sigma2 <- signal_var * (1 - R2_target) / max(R2_target, 1e-6)
  sigma <- sqrt(max(sigma2, 1e-8))

  y <- mu_nl + rnorm(n, mean = 0, sd = sigma)

  # ============================================================
  # 11) Data frame final y ajuste observado
  # ============================================================
  df <- cbind.data.frame(y = y, X_df)

  rhs_fit <- if (fit_with_interactions && length(selected_interactions) > 0) {
    paste(". +", paste(selected_interactions, collapse = " + "))
  } else {
    "."
  }

  formula_fit <- as.formula(paste("y ~", rhs_fit))
  fit_lm <- lm(formula_fit, data = df)
  R2_achieved <- summary(fit_lm)$r.squared

  # ============================================================
  # 12) Salida
  # ============================================================
  out <- list(
    data = df,
    fit_lm = fit_lm,
    truth = list(
      X_mm = X_mm,
      beta_true = beta_true,
      sigma = sigma,
      R2_target = R2_target,
      R2_achieved = R2_achieved,
      nonlinear_vars = nl_vars,
      gamma = gamma,
      signal_var = signal_var,
      interaction_terms = selected_interactions,
      formula_true = formula_true,
      formula_fit = formula_fit
    ),
    settings = list(
      n = n,
      n_cont_norm = n_cont_norm,
      n_cont_skew = n_cont_skew,
      n_bin = n_bin,
      n_fact = n_fact,
      fact_levels = fact_levels,
      rho = rho,
      n_latent = n_latent,
      latent_strength = latent_strength,
      include_interactions = include_interactions,
      interaction_types = interaction_types,
      interaction_frac = interaction_frac,
      interaction_sd = interaction_sd,
      fixed_interactions = fixed_interactions,
      fit_with_interactions = fit_with_interactions
    )
  )

  return(out)
}
