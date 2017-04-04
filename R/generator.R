generator <- function(n, mu, sigma2, nu, lower, upper, dist) {
  if (dist == "Normal") {
    amostra_normal <- inv_TN(n, mu, sigma2, lower, upper)
    return(amostra_normal)
  }
  if (dist == "T") {
    amostra_t <- inv_TT(n, mu, sigma2, nu, lower, upper)
    return(amostra_t)
  }
  if (dist == "Slash") {
    amostra_slash <- inv_TSL(n, mu, sigma2, nu, lower, upper)
    return(amostra_slash)
  }
  if (dist == "CNormal") {
    amostra_CNormal <- inv_TCN(n, mu, sigma2, nu, lower, upper)
    return(amostra_CNormal)
  }
}
