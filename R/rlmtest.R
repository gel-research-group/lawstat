#' Robust Lagrange multiplier test for detecting ARCH/GARCH effect
#'
#' The function performs two resampling techniques to find critical values of the Lagrange multiplier (LM) test,
#' namely permutation and bootstrap, as well as the asymptotic parametric test.
#'
#' @param y  A vector that contains univariate time series observations. Missing values are
#' not allowed.
#' @param x  A matrix that contains control variable time series observations. Missing values are
#' not allowed.
#' @param a.order Order of the autoregressive model which must be a nonnegative integer number. Default value 1.
#' @param sig.level Significance level for testing hypothesis of no ARCH effect. Default value is 0.05.
#' @param crit.type Method of obtaining critical values: "asymptotic" (default) or "nonparametric" options.
#' @param nonp.method Type of resampling if \code{crit.type}="nonparametric": "bootstrap"
#' (default) or "permutation".
#' @param num.resamp Number of resampling replications if \code{crit.type}="nonparametric". Default number is 1,000.
#' @param num.discard Number of bootstrap sample discarded. Default value is 0.
#'
#'
#' @return
#' \item{stats}{Test statistic.}
#' \item{crit.val}{Ressampling based critical value, with a given significance level \code{sig.level}.}
#' \item{p.values}{\code{p-value} of the Lagrange multiplier test.}
#'
#'
#' @references
#' Engle R F. Autoregressive conditional heteroscedasticity with estimates of the variance of United Kingdom inflation[J].
#' Econometrica: Journal of the econometric society, 1982: 987-1007.
#'
#' Gel Y R, Chen B. Robust Lagrange multiplier test for detecting ARCH/GARCH effect using permutation and bootstrap[J].
#' Canadian Journal of Statistics, 2012, 40(3): 405-426.
#'
#' @keywords time series, Lagrange multiplier test, bootstrap, resampling
#'
#' @author Meichen Huang, Tian Jiang, Yulia R. Gel
#'
#' @export
#'
#' @examples
#'
#' Collect daily high price (in USD) of Bitcoin and gold from June 1st, 2016 to May 31st, 2021,
#' and create a dataset with log-returned high price of Bitcoin and gold, named as "GoldBTC".
#' Test for ARCH effects with log-return of high price of Bitcoin and gold with conditions below:
#' H0: No ARCH effect.
#' Ha: There exist ARCH effects.
#'
#' # Let's import dataset "GoldBTC":
#' data(GoldBTC)
#'
#' # Fix seed for reproducible simulations:
#' set.seed(2021)
#'
#' # We want to test ARCH effect in Gold.
#' rlmtest(GoldBTC$Gold, crit.type = "nonparametric",
#'         nonp.method = "bootstrap")
#'
#' ## The list of returned results are below:
#' ## $stat
#' ## [1] 65.84728
#' ## $crit.val
#' ## 95%
#' ## 1.914901
#' ## $p.value
#' ## [1] 0.003
#' ## Since p-value is 0.003, reject null hypothesis at 0.05,
#' ## that is, ARCH effect exists in gold.
#'
#' # We want to test for ARCH effects in Bitcoin.
#' rlmtest(GoldBTC$BTC, crit.type = "nonparametric",
#'         nonp.method = "bootstrap")
#'
#' ## The list of returned results are below:
#' ## $stat
#' ## [1] 8.189278
#' ## $crit.val
#' ## 95%
#' ## 3.541198
#' ## $p.value
#' ## [1] 0.013
#' ## Since p-value is 0.013, reject null hypothesis at 0.05,
#' ## that is, ARCH effects exist in Bitcoin.
#'
#' ## We want to test ARCH effect in gold and Bitcoin.
#' rlmtest(GoldBTC$BTC, GoldBTC$Gold, crit.type = "nonparametric",
#'         nonp.method = "bootstrap")
#'
#' ## The list of returned results are below:
#' ## $stat
#' ## [1] 9.246483
#' ## $crit.val
#' ## 95%
#' ## 3.478484
#' ## $p.value
#' ## [1] 0.015
#' ## Since p-value is 0.015, reject null hypothesis at 0.05,
#' ## that is, ARCH effects exist in Bitcoin controling gold.
#'

rlmtest <- function(y,
                    x = NA,
                    a.order = 1,
                    sig.level = 0.05,
                    crit.type = c("asymptotic", "nonparametric"),
                    nonp.method = c("bootstrap","permutation"),
                    num.resamp = 1000,
                    num.discard = 0)
{
  # the length of time series y
  n <- length(y)

  # regression
  if (is.na(x[1])) {
    ls.fit <- lm(y ~ 1)
  } else {
    ls.fit <- lm(y ~ x)
  }

  # extract residuals
  res <- ls.fit$residuals

  # demean residuals from the regression above ?
  res.demean <- res - mean(res)

  # square the residuals
  res.sq <- res^2
  res.sq.demean <- res.sq-mean(res.sq)

  # run autoregressive ls.fit of a.order on res.sq
  ar.fit <- lm(res.sq.demean[(a.order+1):n]~res.sq.demean[1:(n-a.order)])

  # calculate residuals based on AR(a.order)
  ee <- ar.fit$residuals

  # calculate tau_hat (see the CJS paper)
  SSR0.hat <- sum((res.sq.demean[(a.order+1):n])^2)
  # SSR1.hat <- sum(ee[(a.order+1):n]^2)
  SSR1.hat <- sum(ee^2)
  tau.hat <- (n-a.order)*(1-SSR1.hat/SSR0.hat) # this the observed statistic

  if (crit.type == "asymptotic") {
    pValue <- pchisq(tau.hat, df = a.order, lower.tail = FALSE)
    crit <- qchisq(sig.level, df = a.order, lower.tail = FALSE)
    return(list(stat = tau.hat, crit.val = crit, p.value = pValue))
  } else {
    # resampling test
    tau.bs <- numeric(num.resamp)

    for (j in 1:num.resamp) {

      ### resample residuals
      if(nonp.method == "permutation") {
        # permutation
        e.bs <- sample(res.demean,n,replace=FALSE)
      } else {
        # bootstrap
        e.bs <- sample(res.demean,n+num.discard,replace=TRUE)
      }

      # create a new time series y.bs that mimicks y
      y.bs <- fitted(ls.fit) + e.bs[(1+num.discard):(n+num.discard)]
      # regression ls.fit
      if (is.na(x[1])) {
        ls.fit.bs <- lm(y.bs ~ 1)
      } else {
        ls.fit.bs <- lm(y.bs ~ x)
      }

      # extract residuals
      res.bs <- ls.fit.bs$residuals

      # square the residuals
      res.sq.bs <- (res.bs-mean(res.bs))^2

      res.sq.demean.bs <- res.sq.bs-mean(res.sq.bs)

      # run autoregressive of a.order on res.sq
      ar.fit.bs <- lm(res.sq.demean.bs[(a.order+1):n]~res.sq.demean.bs[1:(n-a.order)])

      # AR(a.order) residuals
      ee.bs <- ar.fit.bs$residuals

      # calculate tau_hat (see the CJS paper)
      SSR0.bs <- sum(res.sq.demean.bs[(a.order+1):n]^2)
      # SSR1.bs <- sum(ee.bs[(a.order+1):n]^2)
      SSR1.bs <- sum(ee.bs^2)
      # this the observed statistic
      tau.bs[j] <- (n-a.order)*(1-SSR1.bs/SSR0.bs) # this the observed statistic
    }

    # associated p-value to report
    pValue.bs = mean(tau.bs > tau.hat)
    # associated critical value to report
    crit.bs = quantile(tau.bs, 1-sig.level)

    return(list(stat = tau.hat, crit.val = crit.bs, p.value = pValue.bs))
  }
}

