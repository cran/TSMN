\name{TSMNgenerator}
\alias{TSMNgenerator}

\title{Generator of Truncated Scale Mixtures of Normal Distributions}

\description{This function generate random TSMN samples from Normal, Student-t, Slash and Contaminated Normal, using the the inverse method.}

\usage{
TSMNgenerator(n, mu, sigma2, nu = NULL, lower = -Inf, upper = Inf, dist = "Normal")
}

\arguments{
\item{n}{Number of observations.}
\item{mu}{Location parameter.}
\item{sigma2}{Scale parameter.}
\item{nu}{Parameter of the scale variable of the SMN family. Must be NULL in case of Normal distribution. Must be a bidimensional vector in case of contaminated normal distribution (CNormal).}
\item{lower}{Lower bounds.}
\item{upper}{Upper bounds.}
\item{dist}{Distribution to be used: "Normal" for Normal model, "T" for Student-t model, "Slash" for slash model and "CNormal" for contaminated Normal model.}
}

\details{For the contaminated Normal distribution, each component of the bidimensional vector "nu" must lie on (0,1).}

\value{
Returns the sample generated according to the parameters.
}

\references{
 Aldo M. Garay, Victor H. Lachos, Heleno Bolfarine, Celso R. Cabral. "Linear censored regression models with scale mixtures of normal distributions".
 Statistical Papers, 2017, vol. 58, issue 1, pages 247-278
}

\author{Eraldo B. dos Anjos Filho \email{ebdaf1@de.ufpe.br} and Aldo M. Garay \email{agaray@de.ufpe.br}

Maintainer: Eraldo B. dos Anjos Filho \email{ebdaf1@de.ufpe.br}
}

\seealso{
\code{\link{TSMNmoments}}
}

\examples{
## A test sample to compare theoretical and empirical moments, considering the following parameters:

mu = 2
sigma2 = 4
nu = 5
lower = -3
upper = 10
dist = "T"
n = 10000

## Theoretical moments with TSMNmoments
theor<-TSMNmoments(mu=mu, sigma2=sigma2, nu=nu, lower=lower, upper=upper, dist=dist)

## Generate the sample with TSMNgenerator to compute the empirical moments
empir<-TSMNgenerator(n=n, mu=mu, sigma2=sigma2, nu=nu, lower=lower, upper=upper, dist=dist)

## Compare the results
data.frame("1st" = c("Theoretic" = theor$EY1, "Empirical" = mean(empir)),
           "2nd" = c("Theoretic" = theor$EY2, "Empirical" = mean(empir^2)),
           "3rd" = c("Theoretic" = theor$EY3, "Empirical" = mean(empir^3)),
           "4th" = c("Theoretic" = theor$EY4, "Empirical" = mean(empir^4)))

}

\keyword{ truncated scale mixtures of normal distribution }
\keyword{ tsmn }
\keyword{ generator of truncated scale mixtures of normal distributions }
