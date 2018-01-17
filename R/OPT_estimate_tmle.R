#' Main TMLE Calculations for the Optimal Individualized Treatment Regime Parameter
#'
#' This function performs all the main TMLE calculations for the context-specific
#' optimal individualized treatment Regime Parameter.
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

tmleOPT <- function(tmledata, maxIter=1000, Qbounds=c(1e-4,1-1e-4)){

  order <- 1/nrow(tmledata)

  #Initial estimate:
  eststep <- OPT_estimate(tmledata)

  for (iter in seq_len(maxIter)){

    updatestep <- OPT_update(eststep$tmledata, Qbounds)
    eststep <- OPT_estimate(updatestep$tmledata)

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

#' Prepare for Optimal Individualized Treatment Regime TMLE Procedure
#'
#' This function performs all the necessary steps in order to call the main TMLE function
#' for the optimal individualized treatment regime, \code{tmleOPT}. It also performs all the
#' calculations necessary to obtain ineference on the parameter.
#'
#' @param obsA Observed exposure (A).
#' @param obsY Observed outcome (Y).
#' @param pA1 Estimate of g, or P(A=1) obtained by calling \code{estSplit}.
#' @param Q0W Estimate of Q(0,W), or E(Y|A=0,W) obtained by calling \code{estSplit}.
#' @param Q1W Estimate of Q(1,W), or E(Y|A=1,W) obtained by calling \code{estSplit}.
#' @param ruleA user-specified rule for exposure under which we want to estimate the context-specific
#' mean of the outcome.
#' @param maxIter maximum number of iterations for the iterative TMLE.
#' @param Qbounds bounds for the Q estimates.
#'
#' @return An object of class \code{tstmle}.
#' \describe{
#' \item{tmlePsi}{Context-specific mean of the outcome under user-specified \code{ruleA} estimated
#' using TMLE methodology.}
#' \item{tmleSD}{Standard deviation of the context-specific mean of the outcome under
#' user-specified \code{ruleA} estimated using TMLE methodology.}
#' \item{tmleCI}{Confidence interval of the context-specific mean of the outcome under
#' user-specified \code{ruleA} estimated using TMLE methodology.}
#' \item{IC}{Influence curve for the context-specific parameter under user-specified \code{ruleA}.}
#' \item{steps}{Number of steps until convergence of the iterative TMLE.}
#' \item{initialData}{Initial estimates of g and Q, and observed A and Y.}
#' \item{tmleData}{Final updates estimates of g, Q and clever covariates.}
#' \item{all}{Full results of \code{tmleOPT}.}
#' }
#'
#' @importFrom stats qnorm
#'
#' @export
#

ruletmle<-function(obsA,obsY,pA1,Q0W,Q1W,ruleA,maxIter=1000,Qbounds=c(1e-4,1-1e-4)){

  #Get samples that follow the rule
  #(1 for samples that follow the rule in the observed data)
  A<-as.numeric(obsA==ruleA)

  #If rule is A=0, get p(A=0)
  pA<-pA1*ruleA + (1-pA1)*(1-ruleA)

  #Need E(Y|A=0,W) if rule is 0
  Q<-Q1W*ruleA + Q0W*(1-ruleA)

  #Create appropriate data.frame and run tmle
  tmledata<-data.frame(A=A, Y=obsY, gk=pA, Qk=Q)
  res<-tmleOPT(tmledata, maxIter, Qbounds)

  #Get inference:
  psi<-res$psi
  sd<-sqrt(res$ED2)/sqrt(length(obsA))
  lower <- psi - stats::qnorm(0.975) * sd
  upper <- psi + stats::qnorm(0.975) * sd

  return(list(tmlePsi=psi, tmleSD=sd, tmleCI=c(lower,upper),
              IC=res$IC, steps=res$steps, initialData=tmledata, tmleData=res$tmledata,all=res))
}

#' Estimate function for the context-specific parameter under user-specified rule
#'
#' Function to estimate the context-specific parameter under user-specified rule and
#' coresponding relevant parts of the influence curve for the single intervention.
#'
#' @param tmledata \code{data.frame} containing all observed values for the A and Y node,
#' and either initial or updated estimates for g and Q.
#'

OPT_estimate <- function(tmledata) {

  psi <- mean(tmledata$Qk)

  tmledata$H1 <- with(tmledata, (1/gk))
  tmledata$HA <- with(tmledata, (A*H1))

  # influence curves
  Dstar_psi <- with(tmledata, HA * (Y - Qk) + Qk - psi)

  list(tmledata = tmledata, ests = c(psi = psi), Dstar = list(Dstar_psi = Dstar_psi))

}

#' Update function for the context-specific parameter under user-specified rule
#'
#' Function to update the Q part of the likelihood using the specified fluctuation model for the
#' context-specific parameter under user-specified rule.
#'
#' @param data \code{data.frame} containing all observed values for the A and Y node,
#' and either initial or updated estimates for g and Q.
#' @param Qbounds bounds for the Q estimates.
#'
#' @importFrom stats plogis qlogis
#'

OPT_update <- function(data, Qbounds) {

  subset <- with(data, which(0 < Qk & Qk < 1 & A == 1))
  eps <- 0

  if (length(subset) > 0) {
    #fluctuate Q
    data$Qktrunc <- bound(with(data, Qk), Qbounds)
    qfluc <- fluctuate(data, Y ~ -1 + offset(stats::qlogis(Qktrunc)) + HA)
    eps <- qfluc$eps
    data$Qk <- with(data, stats::plogis(stats::qlogis(Qktrunc) + HA * eps))
  }

  list(tmledata = data, coefs = c(eps))

}

#' Fluctuate function
#'
#' Function to estimate logistic parametric submodel and get updated estmate logistic fluctuation.
#'
#' @param tmledata \code{data.frame} containing all observed values for the A and Y node,
#' and either initial or updated estimates for g and Q.
#' @param flucmod fluctuation model used.
#' @param subset subset of the data used for fluctuation.
#'
#' @importFrom stats glm
#'

fluctuate <- function(tmledata, flucmod, subset = seq_len(nrow(tmledata))) {
  suppressWarnings({
    fluc <- stats::glm(flucmod, data = tmledata[subset, ], family = "binomial")
  })
  list(eps = stats::coef(fluc)[1])
}
