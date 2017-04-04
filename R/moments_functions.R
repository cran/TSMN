AcumPearsonVII <- function(y, mu, sigma2, nu, delta) {
  Acum <- z <- vector(mode = "numeric", length = length(y))
  sigma2a <- sigma2 * (delta/nu)
  for (i in 1:length(y)) {
    z[i] <- (y[i] - mu)/sqrt(sigma2a)
    Acum[i] <- pt(z[i], df = nu)
  }
  return(Acum)
}

AcumSlash <- function(y, mu, sigma2, nu) {
  Acum <- z <- vector(mode = "numeric", length = length(y))
  for (i in 1:length(y)) {
    z[i] <- (y[i] - mu)/sqrt(sigma2)
    f1 <- function(u) nu * u^(nu - 1) * pnorm(z[i] * sqrt(u))
    Acum[i] <- integrate(f1, 0, 1)$value
  }
  return(Acum)
}

AcumCNormal <- function(y, mu, sigma2, nu) {
  Acum <- vector(mode = "numeric", length = length(y))
  for (i in 1:length(y)) {
    eta <- nu[1]
    gama <- nu[2]
    Acum[i] <- eta * pnorm(y[i], mu, sqrt(sigma2/gama)) + (1 - eta) * pnorm(y[i], mu, sqrt(sigma2))
  }
  return(Acum)
}

dPearsonVII <- function(y, mu, sigma2, nu, delta) {
  f <- z <- vector(mode = "numeric", length = length(y))
  sigma2a <- sigma2 * (delta/nu)
  for (i in 1:length(y)) {
    z[i] <- (y[i] - mu)/sqrt(sigma2a)
    f[i] <- dt(z[i], df = nu)/sqrt(sigma2a)
  }
  return(f)
}

dSlash <- function(y, mu, sigma2, nu) {
  resp <- z <- vector(mode = "numeric", length = length(y))
  for (i in 1:length(y)) {
    z[i] <- (y[i] - mu)/sqrt(sigma2)
    f2 <- function(u) nu * u^(nu - 0.5) * dnorm(z[i] * sqrt(u))/sqrt(sigma2)
    resp[i] <- integrate(f2, 0, 1)$value
  }
  return(resp)
}

dCNormal <- function(y, mu, sigma2, nu) {
  Acum <- vector(mode = "numeric", length = length(y))
  for (i in 1:length(y)) {
    eta <- nu[1]
    gama <- nu[2]
    Acum[i] <- eta * dnorm(y[i], mu, sqrt(sigma2/gama)) + (1 - eta) * dnorm(y[i], mu, sqrt(sigma2))
  }
  return(Acum)
}

GamaInc <- function(b, x) {
  res <- vector(mode = "numeric", length = length(x))
  f <- function(t) {
    exp(-t) * t^(b - 1)
  }
  res <- integrate(f, 0, x)$value
  return(res)
}

E_phi <- function(r, a, nu, delta, dist) {
  if (dist == "Normal") {
    resp <- dnorm(a)
  }
  if (dist == "T") {
    Aux0 <- gamma(0.5 * (nu + 2 * r))
    Aux1 <- gamma(nu/2) * sqrt(2 * pi)
    Aux2 <- Aux0/Aux1
    Aux3 <- (0.5 * nu)^(nu/2)
    Aux4 <- (0.5 * (a^2 + nu))^(-0.5 * (nu + 2 * r))
    resp <- Aux2 * Aux3 * Aux4
  }
  if (dist == "PearsonVII") {
    Aux0 <- gamma(0.5 * (nu + 2 * r))
    Aux1 <- gamma(nu/2) * sqrt(2 * pi)
    Aux2 <- Aux0/Aux1
    Aux3 <- (0.5 * delta)^(nu/2)
    Aux4 <- (0.5 * (a^2 + delta))^(-0.5 * (nu + 2 * r))
    resp <- Aux2 * Aux3 * Aux4
  }
  if (dist == "Slash") {
    Aux0 <- nu/sqrt(2 * pi)
    Aux1 <- (0.5 * a^2)^(-(nu + r))
    Aux2 <- GamaInc(nu + r, 0.5 * a^2)
    resp <- Aux0 * Aux1 * Aux2
  }
  if (dist == "CNormal") {
    Aux0 <- nu[1] * nu[2]^(r) * dnorm(a * sqrt(nu[2]))
    Aux1 <- (1 - nu[1]) * dnorm(a)
    resp <- Aux0 + Aux1
  }
  return(resp)
}

E_Phi <- function(r, a, nu, delta, dist) {
  n <- length(a)
  if (dist == "Normal") {
    resp <- pnorm(a)
  }
  if (dist == "T") {
    Aux0 <- gamma(0.5 * (nu + (2 * r)))
    Aux1 <- gamma(nu/2)
    Aux2 <- Aux0/Aux1
    Aux3 <- (0.5 * nu)^(-r)
    Aux4 <- AcumPearsonVII(a, 0, 1, nu + (2 * r), nu)
    resp <- Aux2 * Aux3 * Aux4
  }
  if (dist == "PearsonVII") {
    Aux0 <- gamma(0.5 * (nu + (2 * r)))
    Aux1 <- gamma(nu/2)
    Aux2 <- Aux0/Aux1
    Aux3 <- (0.5 * delta)^(-r)
    Aux4 <- AcumPearsonVII(a, 0, 1, nu + (2 * r), delta)
    resp <- Aux2 * Aux3 * Aux4
  }
  if (dist == "Slash") {
    Aux0 <- nu/(nu + r)
    Aux1 <- AcumSlash(a, 0, 1, nu + r)
    resp <- Aux0 * Aux1
  }
  if (dist == "CNormal") {
    Aux0 <- nu[2]^(r) * AcumCNormal(a, 0, 1, nu)
    Aux1 <- (1 - nu[2]^(r)) * (1 - nu[1]) * pnorm(a)
    resp <- Aux0 + Aux1
  }
  return(resp)
}
