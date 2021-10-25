#' Ranked Version of von Neumann's Ratio Test for Randomness
#' 
#' \insertCite{Bartels_1982;textual}{lawstat} test for randomness that is based 
#' on the ranked version of von Neumann's ratio (RVN). 
#' Users can choose whether to test against two-sided, negative, 
#' or positive correlation. \code{NA}s from the data are omitted.
#'
#'
#' @param y a numeric vector of data values.
#' @param alternative a character string specifying the alternative hypothesis, 
#' must be one of \code{"two.sided"} (default), \code{"negative.correlated"}, or 
#' \code{"positive.correlated"}.
#'
#'
#' @return A list of class \code{"htest"} with the following components:
#' \item{statistic}{the value of the standardized Bartels statistic.}
#' \item{parameter}{RVN ratio.}
#' \item{p.value}{the \eqn{p}-value for the test.}
#' \item{data.name}{a character string giving the names of the data.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @seealso \code{\link{runs.test}}
#' 
#' @keywords distribution htest
#' 
#' @author Kimihiro Noguchi, Wallace Hui, Yulia R. Gel, Joseph L. Gastwirth, Weiwen Miao
#' 
#' @export
#' @examples
#' ## Simulate 100 observations from an autoregressive model of 
#' ## the first order AR(1)
#' y = arima.sim(n = 100, list(ar = c(0.5)))
#' 
#' ## Test y for randomness
#' bartels.test(y)
#' 
#' ## Sample Output
#' ##
#' ##        Bartels Test - Two sided
#' ## data:  y
#' ## Standardized Bartels Statistic -4.4929, RVN Ratio =
#' ## 1.101, p-value = 7.024e-06
#' 
bartels.test <- function(y, alternative = c("two.sided", "positive.correlated", 
                                            "negative.correlated"))
{
    alternative <- match.arg(alternative)
    DNAME = deparse(substitute(y))
    ##Strip NAs
    y <- na.omit(y)
    ### Calculate the rank of the input data ###
    R <- rank(y)
    sum<-0
    T <- length(y)
    for (i in 1:(T-1)) {sum<-sum+(R[i]-R[i+1])^2}
    num <- sum
    sum <- 0
    for (i in 1:T){sum <- sum + (R[i]-(T+1)/2)^2}
    den <- sum
    
    #### Ratio of Von Neumann ####
    RVN <- num/den
    #statistic<-((RVN-2)*sqrt(5*T*(T+1)*(T-1)^2))/sqrt(4*(T-2)*(5*T^2-2*T-9))
    statistic <- (RVN-2)*sqrt(T)/2
    
    ### One-sided or Two-sided Test ###
    ### Users will select the test alternative. Two-sided test is default ###
    if (alternative == "positive.correlated") {
        p.value = pnorm(statistic)
        METHOD = "Bartels Test - Positive Correlated"
    } else if (alternative == "negative.correlated") {
        p.value = 1-pnorm(statistic)
        METHOD = "Bartels Test - Negative Correlated"
    } else {
        p.value = 2*min(pnorm(statistic), 1-pnorm(statistic))
        alternative = "two.sided"
        METHOD = "Bartels Test - Two sided"
    }
    ### Display Output ###
    PARAMETER = RVN
    names(PARAMETER) = "RVN Ratio"
    STATISTIC = statistic
    names(STATISTIC) = "Standardized Bartels Statistic"
    structure(list(statistic = STATISTIC, parameter = PARAMETER, 
                   p.value = p.value, method = METHOD, 
                   data.name = DNAME), class = "htest")
}

