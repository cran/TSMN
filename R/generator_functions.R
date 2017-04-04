pdfNI <- function(y, mu, sigma2, nu, dist) {
  resp <- matrix(0, length(y), 1)
  if (dist == "Normal") {
    resp <- dnorm((y - mu)/sqrt(sigma2))/sqrt(sigma2)
  }
  if (dist == "T") {
    resp <- dt((y - mu)/sqrt(sigma2), df = nu)/sqrt(sigma2)
  }
  if (dist == "Slash") {
    f <- function(u) {
      2 * nu * u^(nu - 1) * dnorm(y, mu, sqrt(u^(-1) * sigma2)) * pnorm(u^(1/2) * (y - mu)/sqrt(sigma2))
    }
    resp <- integrate(Vectorize(f), 0, 1)$value
  }
  if (dist == "CNormal") {
    resp <- 2 * (nu[1] * dnorm(y, mu, sqrt(sigma2/nu[2])) * pnorm(sqrt(nu[2]) * sigma2^(-1/2) *
                                                                    (y - mu)) + (1 - nu[1]) * dnorm(y, mu, sqrt(sigma2)) * pnorm(sigma2^(-1/2) * (y - mu)))
  }
  return(resp)
}

cdfNI <- function(x, mu, sigma2, nu, dist) {
  if (x == -Inf) {
    x = -500
  }
  resp <- matrix(0, length(x), 1)
  if (dist == "T") {
    cdf <- function(y) {
      z = (y - mu)/sqrt(sigma2)
      cdf = 2 * dt(z, df = nu) * pt(sqrt(nu + 1) * z/sqrt(nu + z^2), df = nu + 1)/sqrt(sigma2)
    }
  }
  if (dist == "Slash") {
    cdf <- function(y) {
      f <- function(u) 2 * nu * u^(nu - 1) * dnorm(y, mu, sqrt(u^(-1) * sigma2)) * pnorm(u^(1/2) *
                                                                                           (y - mu)/sqrt(sigma2))
      cdf <- integrate(Vectorize(f), 0, 1)$value
    }
  }
  if (dist == "CNormal") {
    cdf <- function(y) {
      z = (y - mu)/sqrt(sigma2)
      z2 = z * sqrt(nu[2])
      cdf = (2 * nu[1] * dnorm(z2) * pnorm(z2)/sqrt(sigma2/nu[2])) + (2 * (1 - nu[1]) * dnorm(z) *
                                                                        pnorm(z)/sqrt(sigma2))
    }
  }
  for (i in 1:length(x)) {
    resp[i] <- integrate(Vectorize(cdf), -Inf, x[i])$value
  }
  return(resp)
}

inv_TN <- function(n, mu, sigma2, lower, upper) {
  u <- runif(n)
  sigma <- sqrt(sigma2)
  aux <- u * (pnorm((upper - mu)/sqrt(sigma2)) - pnorm((lower - mu)/sqrt(sigma2))) + pnorm((lower - mu)/sqrt(sigma2))
  amostra_x <- mu + sqrt(sigma2) * qnorm(aux)
  return(amostra_x)
}

inv_TT <- function(n, mu, sigma2, nu, lower, upper) {
  u <- runif(n)
  zb <- (upper - mu)/sqrt(sigma2)
  za <- (lower - mu)/sqrt(sigma2)
  aux <- u*(pt(zb,df = nu)-pt(za,df=nu))+pt(za,df=nu)
  amostra.x <- mu + sqrt(sigma2) * qt(aux, df = nu)
  return(amostra.x)
}

inv_TSL <- function(n, mu, sigma2, nu, lower, upper) {
  u <- rbeta(n, shape1 = nu, shape2 = 1)
  a1 <- (lower - mu) * sqrt(u)
  b1 <- (upper - mu) * sqrt(u)
  aux0 <- am.x <- c()
  for (i in 1:n) {
    aux0[i] <- inv_TN(1, 0, sigma2, lower = a1[i], upper = b1[i])
    am.x[i] <- mu + (u[i])^(-1/2) * aux0[i]
  }
  return(am.x)
}

inv_TCN <- function(n, mu, sigma2, nu, lower, upper) {
  p <- runif(n)
  u <- rep(1, n)
  u[p < nu[1]] <- nu[2]
  a1 <- (lower - mu) * sqrt(u)
  b1 <- (upper - mu) * sqrt(u)
  aux0 <- c()
  am.x <- c()
  for (i in 1:n) {
    aux0[i] <- inv_TN(1, 0, sigma2, lower = a1[i], upper = b1[i])
    am.x[i] <- mu + (u[i])^(-1/2) * aux0[i]
  }
  return(am.x)
}
