## Normalizando nombres
normaliza_terminos <- function(term) {
  term2 <- term
  
  # binarios codificados como factor en SD: xB11 -> xB1, xB21 -> xB2
  term2 <- gsub("^xB([0-9]+)1$", "xB\\1", term2)
  
  term2
}

## Para ordenar términos
orden_terminos <- function(term) {
  grupo <- dplyr::case_when(
    term == "(Intercept)" ~ 0,
    str_detect(term, ":") ~ 5,
    str_detect(term, "^xN\\d+$") ~ 1,
    str_detect(term, "^xS\\d+$") ~ 2,
    str_detect(term, "^xB\\d+$") ~ 3,
    str_detect(term, "^xF\\d+L\\d+$") ~ 4,
    TRUE ~ 99
  )

  num_base <- dplyr::case_when(
    str_detect(term, "^xN\\d+$") ~ as.numeric(str_extract(term, "\\d+")),
    str_detect(term, "^xS\\d+$") ~ as.numeric(str_extract(term, "\\d+")),
    str_detect(term, "^xB\\d+$") ~ as.numeric(str_extract(term, "\\d+")),
    str_detect(term, "^xF\\d+L\\d+$") ~ as.numeric(str_extract(term, "(?<=xF)\\d+")),
    TRUE ~ NA_real_
  )

  num_nivel <- case_when(
    str_detect(term, "^xF\\d+L\\d+$") ~ as.numeric(str_extract(term, "(?<=L)\\d+")),
    TRUE ~ NA_real_
  )

  inter_rank <- seq_along(term)
  inter_rank[!str_detect(term, ":")] <- Inf

  data.frame(
    term = term,
    grupo = grupo,
    num_base = ifelse(is.na(num_base), Inf, num_base),
    num_nivel = ifelse(is.na(num_nivel), Inf, num_nivel),
    inter_rank = inter_rank
  ) %>%
    arrange(grupo, num_base, num_nivel, inter_rank)
}

## Preparando tablas para plots
prep_infer_plot_data <- function(tab) {
  ord <- orden_terminos(normaliza_terminos(tab$term))

  tab2 <- tab %>%
    mutate(
      term_original = term,
      term_std = normaliza_terminos(term)
    ) %>%
    left_join(
      ord %>% rename(term_std = term),
      by = "term_std"
    ) %>%
    arrange(grupo, num_base, num_nivel, inter_rank) %>%
    mutate(
      term_plot = factor(term_original, levels = term_original),
      reject_sd = coalesce(power, type1),
      reject_type = ifelse(is.na(power), "Type I", "Power")
    )

  tab2
}

## PLOT sesgo
plot_bias_terms <- function(tab) {
  df <- prep_infer_plot_data(tab)

  ggplot(df, aes(x = term_plot, y = bias_beta)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, color = "gray40") +
    geom_point(size = 2.8, color = "black") +
    geom_line(aes(group = 1), linewidth = 0.4, color = "black") +
    labs(x = NULL, y = "Bias") +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )
}

## PLOT cover
plot_coverage_terms <- function(tab, nominal = 0.95) {
  df <- prep_infer_plot_data(tab)

  ggplot(df, aes(x = term_plot, y = cover_95)) +
    geom_hline(yintercept = nominal, linetype = "dashed", linewidth = 0.4, color = "gray40") +
    geom_point(size = 2.8, color = "black") +
    geom_line(aes(group = 1), linewidth = 0.4, color = "black") +
    labs(x = NULL, y = "Coberage") +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )
}

## PLOT mse beta
plot_mse_terms <- function(tab) {
  df <- prep_infer_plot_data(tab)

  ggplot(df, aes(x = term_plot, y = mse_beta)) +
    geom_point(size = 2.8, color = "black") +
    geom_line(aes(group = 1), linewidth = 0.4, color = "black") +
    labs(x = NULL, y = "MSE") +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )
}

## PLOT power compare 
plot_power_od_vs_sd <- function(tab_sd, tab_od_power) {
  df_sd <- prep_infer_plot_data(tab_sd) %>%
    select(term_original, term_plot, term_std, reject_sd, reject_type)

  df_od <- tab_od_power %>%
    mutate(term_std = normaliza_terminos(term)) %>%
    select(term, term_std, power) %>%
    rename(power_od = power)

  df_plot <- df_sd %>%
    left_join(df_od %>% dplyr::select(term_std, power_od), by = "term_std") %>%
    pivot_longer(
      cols = c(power_od, reject_sd),
      names_to = "source",
      values_to = "value"
    ) %>%
    mutate(
      source = ifelse(source == "power_od", "OD", "SD")
    )

  ggplot(df_plot, aes(x = term_plot, y = value, shape = source, linetype = source, group = source)) +
    geom_point(size = 2.8, color = "black", position = position_dodge(width = 0.2)) +
    geom_line(color = "black", linewidth = 0.4, position = position_dodge(width = 0.2)) +
    scale_shape_manual(values = c(16, 1)) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    labs(
      x = NULL,
      y = "Power",
      shape = NULL,
      linetype = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )
}

todo <- prep_infer_plot_data(res$by_term) %>%
  select(term = term_original, bias_beta, var_beta, mse_beta, cover_95,
                mean_ci_len, type1, power)

p_bias  <- plot_bias_terms(todo)
p_cov   <- plot_coverage_terms(todo)
p_mse   <- plot_mse_terms(todo)
p_power <- plot_power_od_vs_sd(todo, pow_tbl)

p_bias
p_cov
p_mse
p_power
