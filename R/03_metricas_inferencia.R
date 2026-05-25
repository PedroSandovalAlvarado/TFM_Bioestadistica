infer_metrics_mSD <- function(fits_obj, beta_true = NULL, alpha = 0.05) {
  od_tab  <- fits_obj$od
  sd_tabs <- fits_obj$sds
  m <- length(sd_tabs)

  # Apilar SD en formato largo
  sd_long <- do.call(rbind, Map(function(tab, k) {
    tab$k <- k
    tab
  }, sd_tabs, seq_len(m)))

  # Join SD con OD por término (quedan estimate_sd/od, p_sd/od, etc.)
  merged <- merge(
    sd_long,
    od_tab[, c("term", "estimate", "p", "sig_od")],
    by = "term",
    suffixes = c("_sd", "_od")
  )

  # Renombrar a nombres consistentes
  names(merged)[names(merged) == "estimate_sd"] <- "beta_sd"
  names(merged)[names(merged) == "estimate_od"] <- "beta_od"
  names(merged)[names(merged) == "p_sd"] <- "p_sd"
  names(merged)[names(merged) == "p_od"] <- "p_od"

  # Métricas por (term, k)
  merged$delta_beta <- merged$beta_sd - merged$beta_od
  merged$cover_od   <- (merged$beta_od >= merged$lwr) & (merged$beta_od <= merged$upr)
  merged$rej_sd     <- merged$p_sd < alpha

  # ---------- Resumen por término ----------
  agg_mean <- aggregate(cbind(
    beta_sd = merged$beta_sd,
    delta   = merged$delta_beta,
    cover   = merged$cover_od,
    ci_len  = merged$ci_len,
    rej     = merged$rej_sd
  ) ~ term, data = merged, FUN = mean)

  agg_var <- aggregate(cbind(
    beta_sd = merged$beta_sd,
    delta   = merged$delta_beta
  ) ~ term, data = merged, FUN = var)

  # Renombrar columnas para evitar ambigüedad
  names(agg_mean)[names(agg_mean) == "beta_sd"] <- "mean_beta_sd"
  names(agg_mean)[names(agg_mean) == "delta"]   <- "bias_beta"
  names(agg_mean)[names(agg_mean) == "cover"]   <- "cover_95"
  names(agg_mean)[names(agg_mean) == "ci_len"]  <- "mean_ci_len"
  names(agg_mean)[names(agg_mean) == "rej"]     <- "reject_rate"

  names(agg_var)[names(agg_var) == "beta_sd"] <- "var_beta"
  names(agg_var)[names(agg_var) == "delta"]   <- "var_delta"

  out <- merge(agg_mean, agg_var, by = "term")

  # MSE = Var(delta) + Bias^2
  out$mse_beta <- out$var_delta + out$bias_beta^2

  # Añadir beta_od, p_od, sig_od (del OD)
  aux_od <- unique(merged[, c("term", "beta_od", "p_od", "sig_od")])
  out <- merge(out, aux_od, by = "term")

  # ---------- Tipo I y Potencia ----------
  if (!is.null(beta_true)) {
    # beta_true: named numeric vector (names = términos, e.g. xN1..xN5)
    bt <- data.frame(term = names(beta_true), beta_true = as.numeric(beta_true))
    out <- merge(out, bt, by = "term", all.x = TRUE)

    out$type1 <- ifelse(abs(out$beta_true) < 1e-12, out$reject_rate, NA)
    out$power <- ifelse(abs(out$beta_true) >= 1e-12, out$reject_rate, NA)
    type_label <- "real (beta_true)"
  } else {
    # proxies según OD
    out$type1 <- ifelse(out$sig_od == FALSE, out$reject_rate, NA)
    out$power <- ifelse(out$sig_od == TRUE,  out$reject_rate, NA)
    type_label <- "proxy (según OD)"
  }

  # ---------- Ranking stability ----------
  ks <- sort(unique(merged$k))
  rho <- sapply(ks, function(k) {
    dd <- merged[merged$k == k, ]
    suppressWarnings(cor(abs(dd$beta_sd), abs(dd$beta_od), method = "spearman"))
  })

  rank_summary <- data.frame(
    metric = c("rank_spearman_mean", "rank_spearman_sd"),
    value  = c(mean(rho, na.rm = TRUE), sd(rho, na.rm = TRUE)),
    note   = type_label
  )

  # ---------- Orden final ----------
  out_final <- out[, c(
    "term",
    "bias_beta",
    "var_beta", "mse_beta",
    "cover_95",
    "mean_ci_len",
    "type1",
    "power"
  )]
  
  out_final <- out_final[order(out_final$term), ]

  list(
    by_term = out_final,
    long = merged,
    rank_rho = rho,
    rank_summary = rank_summary,
    meta = list(m = m, alpha = alpha, typeI_power = type_label)
  )
}
