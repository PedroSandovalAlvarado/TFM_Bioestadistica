fit_inference <- function(data, formula, conf_level = 0.95) {
  fit <- lm(formula, data = data)
  s <- summary(fit)

  coefs <- as.data.frame(s$coefficients)
  coefs$term <- rownames(coefs); rownames(coefs) <- NULL
  names(coefs) <- c("estimate", "se", "t", "p", "term")

  ci <- tryCatch({
    ci0 <- as.data.frame(confint(fit, level = conf_level))
    ci0$term <- rownames(ci0); rownames(ci0) <- NULL
    names(ci0)[1:2] <- c("lwr", "upr")
    ci0
  }, error = function(e) {
    data.frame(term = coefs$term, lwr = NA_real_, upr = NA_real_)
  })

  out <- merge(coefs, ci, by = "term", all.x = TRUE)
  out$ci_len <- out$upr - out$lwr

  list(table = out, fit = fit)
}

run_m_fits <- function(od_df, sd_list, formula, alpha = 0.05, conf_level = 0.95) {
  stopifnot(is.data.frame(od_df))
  stopifnot(is.list(sd_list), length(sd_list) >= 1)

  # OD
  od_obj <- fit_inference(od_df, formula, conf_level)
  od_fit <- od_obj$table
  od_fit <- subset(od_fit, term != "(Intercept)")
  od_fit$sig_od <- od_fit$p < alpha

  # SDs
  fits_sd_obj <- lapply(sd_list, function(sd) {
    obj <- fit_inference(sd, formula, conf_level)
    tab <- obj$table
    tab <- subset(tab, term != "(Intercept)")
    tab$sig_sd <- tab$p < alpha
    list(table = tab, fit = obj$fit)
  })

  list(
    od = od_fit,
    od_fit = od_obj$fit,
    sds = lapply(fits_sd_obj, `[[`, "table"),
    sds_fit = lapply(fits_sd_obj, `[[`, "fit"),
    alpha = alpha,
    conf_level = conf_level
  )
}
