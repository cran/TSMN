pre_moments <- function(mu, sigma2, nu, delta, lower, upper, dist) {
  if ((dist == "PearsonVII") || (dist == "T")) {
    moments_T_PVII(mu, sigma2, nu, delta, lower, upper, dist)
  } else {
    if (dist == "Slash") {
      moments_Slash(mu, sigma2, nu, delta, lower, upper, dist)
    } else {
      moments(mu, sigma2, nu, delta, lower, upper, dist)
    }
  }
}

moments <- function(mu, sigma2, nu, delta, lower, upper, dist) {
  Lim1 <- (lower - mu)/sqrt(sigma2)
  Lim2 <- (upper - mu)/sqrt(sigma2)

  if (dist == "Normal") {
    FNIb <- pnorm(Lim2)
    FNIa <- pnorm(Lim1)
  }
  if (dist == "CNormal") {
    FNIb <- AcumCNormal(Lim2, 0, 1, nu)
    FNIa <- AcumCNormal(Lim1, 0, 1, nu)
  }

  if (dist == "T") {
    FNIb <- pt(Lim2, nu)
    FNIa <- pt(Lim1, nu)
  }
  if (dist == "PearsonVII") {
    FNIb <- AcumPearsonVII(Lim2, 0, 1, nu, delta)
    FNIa <- AcumPearsonVII(Lim1, 0, 1, nu, delta)
  }
  if (dist == "Slash") {
    FNIb <- AcumSlash(Lim2, 0, 1, nu)
    FNIa <- AcumSlash(Lim1, 0, 1, nu)
  }

  K <- 1/(FNIb - FNIa)

  EU0X0 <- K * (E_Phi(1, Lim2, nu, delta, dist) - E_Phi(1, Lim1, nu, delta, dist))
  EU0X1 <- K * (E_phi(-0.5, Lim1, nu, delta, dist) - E_phi(-0.5, Lim2, nu, delta, dist))
  EU0X2 <- K * (E_Phi(-1, Lim2, nu, delta, dist) - E_Phi(-1, Lim1, nu, delta, dist) +
                  Lim1 * E_phi(-0.5,Lim1, nu, delta, dist) - Lim2 * E_phi(-0.5, Lim2, nu, delta, dist))
  EU0X3 <- 2 * K * (E_phi(-1.5, Lim1, nu, delta, dist) - E_phi(-1.5, Lim2, nu, delta, dist)) +
    K * ((Lim1 * Lim1) * E_phi(-0.5, Lim1, nu, delta, dist) - (Lim2 * Lim2) * E_phi(-0.5, Lim2,
                                                                                    nu, delta, dist))
  EU0X4 <- (3 * K * (E_Phi(-2, Lim2, nu, delta, dist) - E_Phi(-2, Lim1, nu, delta, dist))) +
   3 * K * (Lim1 * E_phi(-1.5, Lim1, nu, delta, dist) - Lim2 * E_phi(-1.5, Lim2, nu, delta, dist)) +
  K * ( (Lim1^3 * E_phi(-0.5, Lim1, nu, delta, dist)) - (Lim2^3 * E_phi(-0.5,Lim2, nu, delta, dist)))

  EU0Y0 <- EU0X0
  EU0Y1 <- mu * EU0X0 + sqrt(sigma2) * EU0X1
  EU0Y2 <- EU0X0 * (mu^2) + 2 * mu * sqrt(sigma2) * EU0X1 + sigma2 * EU0X2
  EU0Y3 <- EU0X0 * (mu^3) + 3 * (sqrt(sigma2) * (mu^2) * EU0X1 + (sqrt(sigma2)^2) * mu * EU0X2) +
    (sqrt(sigma2)^3) * EU0X3
  EU0Y4 <- (mu^4) * EU0X0 + 4 * (sqrt(sigma2) * (mu^3) * EU0X1 + (sqrt(sigma2)^3) * mu * EU0X3) +
    6 * (sqrt(sigma2)^2) * (mu^2) * EU0X2 + (sqrt(sigma2)^4) * EU0X4

  Retur <- list(EY1 = EU0Y1, EY2 = EU0Y2, EY3 = EU0Y3, EY4 = EU0Y4)

  return(Retur)
}

moments_T_PVII <- function(mu, sigma2, nu, delta, lower, upper, dist) {
  Lim1 <- (lower - mu)/sqrt(sigma2)
  Lim2 <- (upper - mu)/sqrt(sigma2)

  if (dist == "T") {
    FNIb <- pt(Lim2, nu)
    FNIa <- pt(Lim1, nu)
  }
  if (dist == "PearsonVII") {
    FNIb <- AcumPearsonVII(Lim2, 0, 1, nu, delta)
    FNIa <- AcumPearsonVII(Lim1, 0, 1, nu, delta)
  }

  K <- 1/(FNIb - FNIa)

  if (nu > 1 && nu <= 2) {
    EU0X0 <- K * (E_Phi(1, Lim2, nu, delta, dist) - E_Phi(1, Lim1, nu, delta, dist))
    EU0X1 <- K * (E_phi(-0.5, Lim1, nu, delta, dist) - E_phi(-0.5, Lim2, nu, delta, dist))
    EU0Y0 <- EU0X0
    EU0Y1 <- mu * EU0X0 + sqrt(sigma2) * EU0X1
    print("For 2nd Moment nu must be greater than 2")
    print("For 3rd Moment nu must be greater than 3")
    print("For 4th Moment nu must be greater than 4")
    Retur <- list(EY1 = EU0Y1)
  }

  if (nu > 2 && nu <= 3) {
    EU0X0 <- K * (E_Phi(1, Lim2, nu, delta, dist) - E_Phi(1, Lim1, nu, delta, dist))
    EU0X1 <- K * (E_phi(-0.5, Lim1, nu, delta, dist) - E_phi(-0.5, Lim2, nu, delta, dist))
    EU0X2 <- K * (E_Phi(-1, Lim2, nu, delta, dist) - E_Phi(-1, Lim1, nu, delta, dist) + Lim1 *
                    E_phi(-0.5, Lim1, nu, delta, dist) - Lim2 * E_phi(-0.5, Lim2, nu, delta, dist))
    EU0Y0 <- EU0X0
    EU0Y1 <- mu * EU0X0 + sqrt(sigma2) * EU0X1
    EU0Y2 <- EU0X0 * (mu^2) + 2 * mu * sqrt(sigma2) * EU0X1 + sigma2 * EU0X2
    print("For 3rd Moment nu must be greater than 3")
    print("For 4th Moment nu must be greater than 4")
    Retur <- list(EY1 = EU0Y1, EY2 = EU0Y2)
  }

  if (nu > 3 && nu <= 4) {
    EU0X0 <- K * (E_Phi(1, Lim2, nu, delta, dist) - E_Phi(1, Lim1, nu, delta, dist))
    EU0X1 <- K * (E_phi(-0.5, Lim1, nu, delta, dist) - E_phi(-0.5, Lim2, nu, delta, dist))
    EU0X2 <- K * (E_Phi(-1, Lim2, nu, delta, dist) - E_Phi(-1, Lim1, nu, delta, dist) + Lim1 *
                    E_phi(-0.5, Lim1, nu, delta, dist) - Lim2 * E_phi(-0.5, Lim2, nu, delta, dist))
    EU0X3 <- 2 * K * (E_phi(-1.5, Lim1, nu, delta, dist) - E_phi(-1.5, Lim2, nu, delta, dist)) +
      K * ((Lim1 * Lim1) * E_phi(-0.5, Lim1, nu, delta, dist) - (Lim2 * Lim2) * E_phi(-0.5,
                                                                                      Lim2, nu, delta, dist))
    EU0Y0 <- EU0X0
    EU0Y1 <- mu * EU0X0 + sqrt(sigma2) * EU0X1
    EU0Y2 <- EU0X0 * (mu^2) + 2 * mu * sqrt(sigma2) * EU0X1 + sigma2 * EU0X2
    EU0Y3 <- EU0X0 * (mu^3) + 3 * (sqrt(sigma2) * (mu^2) * EU0X1 + (sqrt(sigma2)^2) * mu * EU0X2) +
      (sqrt(sigma2)^3) * EU0X3
    print("For 4th Moment nu must be greater than 4")
    Retur <- list(EY1 = EU0Y1, EY2 = EU0Y2, EY3 = EU0Y3)
  }

  if (nu > 4) {
    Retur <- moments(mu = mu, sigma2 = sigma2, nu = nu, delta = delta, lower = lower, upper = upper, dist = dist)
  }

  return(Retur)
}

moments_Slash <- function(mu, sigma2, nu, delta, lower, upper, dist) {
  Lim1 <- (lower - mu)/sqrt(sigma2)
  Lim2 <- (upper - mu)/sqrt(sigma2)

  if (dist == "Slash") {
    FNIb <- AcumSlash(Lim2, 0, 1, nu)
    FNIa <- AcumSlash(Lim1, 0, 1, nu)
  }

  K <- 1/(FNIb - FNIa)

  if (nu > 1 && nu <= 1.5) {
    EU0X0 <- K * (E_Phi(1, Lim2, nu, delta, dist) - E_Phi(1, Lim1, nu, delta, dist))
    EU0X1 <- K * (E_phi(-0.5, Lim1, nu, delta, dist) - E_phi(-0.5, Lim2, nu, delta, dist))
    EU0X2 <- K * (E_Phi(-1, Lim2, nu, delta, dist) - E_Phi(-1, Lim1, nu, delta, dist) + Lim1 *
                    E_phi(-0.5, Lim1, nu, delta, dist) - Lim2 * E_phi(-0.5, Lim2, nu, delta, dist))

    EU0Y0 <- EU0X0
    EU0Y1 <- mu * EU0X0 + sqrt(sigma2) * EU0X1
    EU0Y2 <- EU0X0 * (mu^2) + 2 * mu * sqrt(sigma2) * EU0X1 + sigma2 * EU0X2

    print("For 3th Moment nu must be greater than 1.5")
    print("For 4th Moment nu must be greater than 2")
    Retur <- list(EY1 = EU0Y1, EY2 = EU0Y2)
  }

  if (nu > 1.5 && nu <= 2) {
    EU0X0 <- K * (E_Phi(1, Lim2, nu, delta, dist) - E_Phi(1, Lim1, nu, delta, dist))
    EU0X1 <- K * (E_phi(-0.5, Lim1, nu, delta, dist) - E_phi(-0.5, Lim2, nu, delta, dist))
    EU0X2 <- K * (E_Phi(-1, Lim2, nu, delta, dist) - E_Phi(-1, Lim1, nu, delta, dist) + Lim1 *
                    E_phi(-0.5, Lim1, nu, delta, dist) - Lim2 * E_phi(-0.5, Lim2, nu, delta, dist))
    EU0X3 <- 2 * K * (E_phi(-1.5, Lim1, nu, delta, dist) - E_phi(-1.5, Lim2, nu, delta, dist)) +
      K * ((Lim1 * Lim1) * E_phi(-0.5, Lim1, nu, delta, dist) - (Lim2 * Lim2) * E_phi(-0.5,
                                                                                      Lim2, nu, delta, dist))
    EU0Y0 <- EU0X0
    EU0Y1 <- mu * EU0X0 + sqrt(sigma2) * EU0X1
    EU0Y2 <- EU0X0 * (mu^2) + 2 * mu * sqrt(sigma2) * EU0X1 + sigma2 * EU0X2
    EU0Y3 <- EU0X0 * (mu^3) + 3 * (sqrt(sigma2) * (mu^2) * EU0X1 + (sqrt(sigma2)^2) * mu * EU0X2) +
      (sqrt(sigma2)^3) * EU0X3
    print("For 4th Moment nu must be greater than 2")
    Retur <- list(EY1 = EU0Y1, EY2 = EU0Y2, EY3 = EU0Y3)
  }

  if (nu > 2) {
    Retur <- moments(mu = mu, sigma2 = sigma2, nu = nu, delta = delta, lower = lower, upper = upper, dist = dist)
  }

  return(Retur)
}
