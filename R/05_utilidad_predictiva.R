# ============================================================
# 1. Utilidad predictiva OD vs SD
# ============================================================

predictive_utility <- function(od_data, syn_list, formula, prop_train = 0.7, seed = 123) {
  
  set.seed(seed)
  
  # Armonizar OD
  od_data <- prepare_regression_data(od_data)
  
  # Armonizar todos los SD
  syn_list <- purrr::map(syn_list, prepare_regression_data)
  
  n <- nrow(od_data)
  idx_train <- sample(seq_len(n), size = floor(prop_train * n))
  idx_test  <- setdiff(seq_len(n), idx_train)
  
  od_train <- od_data[idx_train, ]
  od_test  <- od_data[idx_test, ]
  
  # Modelo en OD
  fit_od <- lm(formula, data = od_train)
  pred_od <- predict(fit_od, newdata = od_test)
  
  rmse_od <- sqrt(mean((od_test$y - pred_od)^2))
  mae_od  <- mean(abs(od_test$y - pred_od))
  
  # Modelos en cada SD, evaluados sobre el mismo test OD
  pred_sd <- purrr::map_dfr(seq_along(syn_list), function(i) {
    
    sd_i <- syn_list[[i]]
    
    fit_sd <- lm(formula, data = sd_i)
    pred_i <- predict(fit_sd, newdata = od_test)
    
    tibble(
      m_id = i,
      rmse = sqrt(mean((od_test$y - pred_i)^2)),
      mae  = mean(abs(od_test$y - pred_i))
    )
  })
  
  list(
    od = tibble(source = "OD", rmse = rmse_od, mae = mae_od),
    sd = pred_sd %>% mutate(source = "SD"),
    test_index = idx_test
  )
}

# ============================================================
# 2. Función auxiliar para armonizar tipos de variables
# ============================================================

prepare_regression_data <- function(df) {
  
  df <- as.data.frame(df)
  
  # Convertir binarias a factor con niveles consistentes
  if ("xB1" %in% names(df)) {
    df$xB1 <- factor(df$xB1, levels = c(0, 1))
  }
  
  if ("xB2" %in% names(df)) {
    df$xB2 <- factor(df$xB2, levels = c(0, 1))
  }
  
  # Asegurar factores categóricos con niveles correctos
  if ("xF1" %in% names(df)) {
    df$xF1 <- factor(df$xF1, levels = c("L1", "L2", "L3"))
  }
  
  if ("xF2" %in% names(df)) {
    df$xF2 <- factor(df$xF2, levels = c("L1", "L2", "L3", "L4"))
  }
  
  return(df)
}
