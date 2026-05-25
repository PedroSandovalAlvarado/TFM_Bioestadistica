# Rubin
rubin_combine <- function(fits_obj, conf_level = 0.95) {
  od_tab  <- fits_obj$od
  sd_tabs <- fits_obj$sds
  m <- length(sd_tabs)

  # términos (ya sin intercepto en tu run_m_fits)
  terms <- od_tab$term

  get_vec <- function(trm, col) {
    sapply(sd_tabs, function(tab) tab[tab$term == trm, col])
  }

  alpha <- 1 - conf_level

  rub_list <- lapply(terms, function(trm) {
    q <- get_vec(trm, "estimate")            # betas por SD
    u <- get_vec(trm, "se")^2                # var dentro (within)

    qbar <- mean(q, na.rm = TRUE)
    ubar <- mean(u, na.rm = TRUE)
    b    <- var(q, na.rm = TRUE)             # between
    tvar <- ubar + (1 + 1/m) * b
    se_tot <- sqrt(tvar)

    # df Rubin (Barnard-Rubin approx simplificada)
    if (is.na(b) || b < 1e-12) {
      df <- Inf
    } else {
      r  <- (1 + 1/m) * b / ubar
      df <- (m - 1) * (1 + 1/r)^2
    }

    tstat <- qbar / se_tot
    pval  <- 2 * pt(abs(tstat), df = df, lower.tail = FALSE)

    crit <- qt(1 - alpha/2, df = df)
    lwr <- qbar - crit * se_tot
    upr <- qbar + crit * se_tot

    data.frame(
      term = trm,
      beta_rubin = qbar,
      se_rubin = se_tot,
      df_rubin = df,
      p_rubin = pval,
      lwr_rubin = lwr,
      upr_rubin = upr
    )
  })

  rub <- do.call(rbind, rub_list)

  # unir con OD (para comparar OD vs Rubin)
  out <- merge(
    od_tab[, c("term","estimate","se","p","lwr","upr")],
    rub,
    by = "term"
  )
  names(out)[names(out)=="estimate"] <- "beta_od"
  names(out)[names(out)=="se"] <- "se_od"
  names(out)[names(out)=="p"] <- "p_od"
  names(out)[names(out)=="lwr"] <- "lwr_od"
  names(out)[names(out)=="upr"] <- "upr_od"

  out
}
