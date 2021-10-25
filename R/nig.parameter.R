#' Generate Parameters for the Normal Inverse Gaussian (NIG) Distribution
#' 
#' Produce four parameters, alpha (tail heavyness), beta (asymmetry), 
#' delta (scale), and mu (location) from the four variables: mean, variance, 
#' kurtosis, and skewness.
#' 
#' @details The parameters are generated with three conditions: 
#' 1) \eqn{3\times kurtosis > 5\times skewness^2}; 
#' 2) \eqn{skewness > 0}, and
#' 3) \eqn{variance > 0}.
#' See \insertCite{Atkinson_1982;textual}{lawstat}, 
#' \insertCite{BarndorffNielsen_Blaesild_1983;textual}{lawstat}, and 
#' \insertCite{Noguchi_Gel_2010;textual}{lawstat}.
#'
#' @param mean mean of the NIG distribution.
#' @param variance variance of the NIG distribution.
#' @param kurtosis excess kurtosis of the NIG distribution.
#' @param skewness skewness of the NIG distribution.
#'
#'
#' @return A list with the following numeric components:
#' \item{alpha}{tail-heavyness parameter of the NIG distribution.}
#' \item{beta}{asymmetry parameter of the NIG distribution.}
#' \item{delta}{scale parameter of the NIG distribution.}
#' \item{mu}{location parameter of the NIG distribution.}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @seealso \code{\link[fBasics]{rnig}}
#' 
#' @keywords distribution
#' 
#' @author Kimihiro Noguchi, Yulia R. Gel
#' 
#' @export
#' @examples
#' library(fBasics)
#' test <- nig.parameter(0, 2, 5, 1)
#' random <- rnig(1000000, alpha = test$alpha, beta = test$beta, 
#'                mu = test$mu, delta = test$delta)
#' mean(random)
#' var(random)
#' kurtosis(random)
#' skewness(random)
#' 
`nig.parameter` <- function(mean = mean,
                            variance = variance,
                            kurtosis = kurtosis,
                            skewness = skewness)
{
    ### stop the code if the parameters do not satisfy the constraints. ###
    if (3 * kurtosis <= 5 * skewness^2) {
        stop("3*kurtosis - 5*skewness^2 must be greater than 0.")
    }
    if (skewness < 0) {
        stop("skewness must be nonnegative.")
    }
    if (variance < 0) {
        stop("variance must be greater than 0.")
    }
    ### parameter calculations ###
    alpha <-
        sqrt(9 * (3 * kurtosis - 4 * skewness ^ 2) / (variance * (3 * kurtosis -
                                                                      5 * skewness ^ 2) ^ 2))
    beta <- sqrt(9 * skewness ^ 2 / (variance * (3 * kurtosis - 5 * skewness ^
                                                     2) ^ 2))
    delta <-
        3 * sqrt(variance * (3 * kurtosis - 5 * skewness ^ 2)) / (3 * kurtosis -
                                                                      4 * skewness ^ 2)
    mu <- mean - delta * beta / sqrt(alpha ^ 2 - beta ^ 2)
    ### display output ###
    return(list(alpha = alpha, beta = beta, delta = delta, mu = mu))
}
