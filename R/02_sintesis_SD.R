# Función para generar SD desde OD

make_syn_from_od <- function(od, m = 1, method = "normrank", seed = 123) {
  stopifnot("y" %in% names(od))
  set.seed(seed)

  vars <- names(od)
  vs <- c(setdiff(vars, "y"), "y")

  pm <- matrix(1, nrow = length(vars), ncol = length(vars),
               dimnames = list(vars, vars))
  diag(pm) <- 0

  # Evitar que y prediga a los X
  pm[setdiff(vars, "y"), "y"] <- 0

  # Permitir que X prediga a y
  pm["y", setdiff(vars, "y")] <- 1

  # method puede ser un solo string o un vector nombrado
  if (length(method) == 1) {
    meth <- rep(method, length(vars))
    names(meth) <- vars
  } else {
    meth <- method[vars]
  }

  syn_obj <- syn(
    data = od,
    m = m,
    method = meth,
    predictor.matrix = pm,
    visit.sequence = vs,
    seed = seed
  )

  sd_list <- syn_obj$syn
  if (is.data.frame(sd_list)) sd_list <- list(sd_list)

  list(
    syn_obj = syn_obj,
    sd_list = sd_list,
    predictor.matrix = pm,
    visit.sequence = vs,
    method = meth
  )
}
