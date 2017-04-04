TSMNgenerator <- function(n, mu, sigma2, nu = NULL, lower = -Inf, upper = Inf, dist = "Normal"){
  if (length(mu) > 1)
    stop("mu parameter must be a scalar.")

  if (length(sigma2) > 1)
    stop("sigma2 parameter must be a scalar.")
  if (sigma2 <= 0)
    stop("Scale parameter must be positive.")

  if (lower > upper)
    stop("lower must be lower than upper.")
  if (length(lower) > 1 || length(upper) > 1)
    stop("Range values must be scalar.")

  if ((dist != "T") && (dist != "Normal") && (dist != "PearsonVII") && (dist != "Slash") && (dist !=
                                                                                             "CNormal"))
    stop("Distribution family not supported. Check documentation!")

  if (dist == "CNormal") {
    if (length(nu) != 2)
      stop("nu must be a bidimensional vector in case of Contaminated Normal distribution.")
    if (nu[1] <= 0 || nu[1] >= 1)
      stop("nu[1] must lies in (0,1).")
    if (nu[2] <= 0 || nu[2] >= 1)
      stop("nu[2] must lies in (0,1).")
  }

  if (lower <= -1e+05) {
    lower = -1e+05
  }
  if (upper >= 1e+05) {
    upper = 1e+05
  }

  if (dist != "CNormal" && dist != "Normal"){
    if (length(nu) > 1)
      stop("nu parameter must be a scalar.")
    if (length(nu) == 0)
      stop("initial value for nu parameter must be provided in case of Pearson VII, T and Slash.")
    if (nu <= 1)
      stop("nu must be greater than 1 for PearsonVII, T and Slash")
    if(nu >= 30){
      nu = 30
    }
  }

  if (mu > upper || mu < lower){
    stop("mu must be between upper and lower limits")
  }

  generator(n = n, mu = mu, sigma2 = sigma2, nu = nu, lower = lower, upper = upper, dist = dist)
}
