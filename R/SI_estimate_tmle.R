#' Main TMLE Calculations for the Single Intervention Context-Specific ATE Paramater
#'
#' This function performs all the main TMLE calculations for the single intervention
#' Context-Specific ATE Paramater.
#'
#' @param tmledata \code{data.frame} containing all observed values for the A and Y node,
#' estimate of g, Q, as well as Q evaluated when A=1 and A=0.
#' @param maxIter maximum number of iterations for the iterative TMLE.
#' @param Qbounds bounds for the Q estimates.
#'
#' @return An object of class \code{tstmle}.
#' \describe{
#' \item{tmledata}{Final updates estimates of g, Q and clever covariates.}
#' \item{psi}{Average treatment effect estimated using TMLE.}
#' \item{steps}{Number of steps until convergence of the iterative TMLE.}
#' \item{IC}{Influence function.}
#' \item{ED}{Mean of the final influence function.}
#' \item{ED2}{Mean of the squared final influence function.}
#' }
#'
#' @export
#

tmleSI <- function(tmledata, maxIter=1000, Qbounds=c(1e-4,1-1e-4)){

  order <- 1/nrow(tmledata)

  #Initial estimate:
  eststep <- SI_estimate(tmledata)

  for (iter in seq_len(maxIter)){

    updatestep <- SI_update(eststep$tmledata, Qbounds)
    eststep <- SI_estimate(updatestep$tmledata)

    ED <- sapply(eststep$Dstar, mean)

    if (all(abs(ED) < order)) {
      converge <- T
      break
    }
  }

  ED2 <- sapply(eststep$Dstar, function(x) mean(x^2))

  return(list(tmledata = eststep$tmledata, psi = eststep$ests, steps = iter,
       IC = eststep$Dstar, ED = ED, ED2 = ED2))
}

#' Fluctuate function
#'
#' Function to estimate logistic parametric submodel and get updated estmate logistic fluctuation
#'
#' @param tmledata \code{data.frame} containing all observed values for the A and Y node,
#' estimate of g, Q, as well as Q evaluated when A=1 and A=0.
#' @param flucmod fluctuation model used.
#' @param subset subset of the data used for fluctuation.
#'
#' @importFrom stats glm
#'

fluctuate <- function(tmledata, flucmod, subset = seq_len(nrow(tmledata))) {
  suppressWarnings({
    fluc <- stats::glm(flucmod, data = tmledata[subset, ], family = "binomial")
  })
  list(eps = coef(fluc)[1])
}

#' Update function
#'
#' Function to update the Q part of the likelihood using the specified fluctuation model/
#'
#' @param tmledata \code{data.frame} containing all observed values for the A and Y node,
#' estimate of g, Q, as well as Q evaluated when A=1 and A=0.
#' @param Qbounds bounds for the Q estimates.
#'
#' @importFrom stats plogis qlogis
#'

SI_update <- function(data, Qbounds) {

  subset <- with(data, which(0 < Qk & Qk < 1))
  eps <- 0

  if (length(subset) > 0) {
    #fluctuate Q
    data$Qktrunc <- bound(with(data, Qk), Qbounds)
    qfluc <- fluctuate(data, Y ~ -1 + offset(stats::qlogis(Qktrunc)) + HA)
    eps <- qfluc$eps
    data$Qk <- with(data, stats::plogis(stats::qlogis(Qktrunc) + HA * eps))
    data$Q1k <- with(data, stats::plogis(stats::qlogis(Q1k) + H1 * eps))
    data$Q0k <- with(data, stats::plogis(stats::qlogis(Q0k) + H0 * eps))
  }

  list(tmledata = data, coefs = c(eps))

}

#' Estimate function
#'
#' Function to estimate the ATE parameter and coresponding relevant parts of the influence curve
#' for the single intervention contex-specific parameter.
#'
#' @param tmledata \code{data.frame} containing all observed values for the A and Y node,
#' estimate of g, Q, as well as Q evaluated when A=1 and A=0.
#'

SI_estimate <- function(tmledata) {

  psi <- mean(tmledata$Q1k - tmledata$Q0k)

  tmledata$H1 <- with(tmledata, (1/gk))
  tmledata$H0 <- with(tmledata, -1/(1 - gk))
  tmledata$HA <- with(tmledata, (A * H1 + (1 - A) * H0))

  # influence curves
  Dstar_psi <- with(tmledata, HA * (Y - Qk) + Q1k - Q0k - psi)

  list(tmledata = tmledata, ests = c(psi = psi), Dstar = list(Dstar_psi = Dstar_psi))

}
