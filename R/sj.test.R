#' Test of Normality -- SJ Test
#' 
#' Perform the robust directed test of normality, which is based on the ratio of the 
#' classical standard deviation \eqn{S} to the robust standard deviation \eqn{J} 
#' (Average Absolute Deviation from the Median, MAAD) of the sample data. 
#' See \insertCite{Gel_etal_2007;textual}{lawstat}.
#' 
#'
#' @param x a numeric vector of data values.
#' @param crit.values a character string specifying how the critical values should be 
#' obtained, i.e., approximated by the \eqn{t}-distribution (default) or empirically.
#' @param N number of Monte Carlo simulations for the empirical critical values.
#'
#'
#' @return A list of class \code{"htest"} with the following components:
#' \item{statistic}{the standardized test statistic.}
#' \item{p.value}{the \eqn{p}-value.}
#' \item{parameter}{the ratio of the classical standard deviation \eqn{S} to
#' the robust standard deviation \eqn{J}.}
#' \item{data.name}{a character string giving the name of the data.}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @seealso \code{\link{rqq}}, \code{\link{rjb.test}}, 
#' \code{\link[tseries]{jarque.bera.test}}
#' 
#' @keywords distribution htest robust
#' 
#' @author Wallace Hui, Yulia R. Gel, Joseph L. Gastwirth, Weiwen Miao
#' 
#' @export
#' @examples
#' data(bias)
#' sj.test(bias)
#' 
`sj.test` <- function(x,
                      crit.values = c("t.approximation", "empirical"),
                      N = 0)
{
    crit.values = match.arg(crit.values)
    if ((crit.values == "empirical") & (N == 0)) {
        stop("number of Monte Carlo simulations N should be provided for the empirical critical values")
    }
    ### SJ Test - New Directional Test
    DNAME = deparse(substitute(x))
    n <- length(x)
    J <- sqrt(pi / 2) * mean(abs(x - median(x)))
    x <- sort(x)
    cw1 <- sd(x) / J
    statistic = sqrt(n) * (cw1 - 1) / sqrt((pi - 3) / 2)
    if (crit.values == "empirical") {
        #### computes empirical critical values for the SJ statistic####
        sj <- double(N)
        for (k in 1:N) {
            e <- rnorm(length(x), mean = 0, sd = sqrt(1))
            J <- sqrt(pi / 2) * mean(abs(e - median(e)))
            sj[k] <- sd(e) / J
        }
        y <- sort(sj)
        if (cw1 >= max(y)) {
            p.value = 0
        } else if (cw1 <= min(y)) {
            p.value = 1
        } else {
            bn <- which(y == min(y[I(y >= cw1)]))
            an <- which(y == max(y[I(y < cw1)]))
            a <- max(y[I(y < cw1)])
            b <- min(y[I(y >= cw1)])
            pa <- (an - 1) / (N - 1)
            pb <- (bn - 1) / (N - 1)
            alpha <- (cw1 - a) / (b - a)
            p.value = 1 - alpha * pb - (1 - alpha) * pa
        }
    } else if (crit.values == "t.approximation") {
        p.value = 1 - pt(statistic, df = (sqrt(n) + 3) / 2)
    }
    METHOD = "Test of Normality - SJ Test"
    ### Display Output ###
    STATISTIC = statistic
    names(STATISTIC) = "Standardized SJ Statistic"
    PARAMETER = cw1
    names(PARAMETER) = "ratio of S to J"
    structure(
        list(
            statistic = STATISTIC,
            parameter = PARAMETER,
            p.value = p.value,
            method = METHOD,
            data.name = DNAME
        ),
        class = "htest"
    )
    
}
