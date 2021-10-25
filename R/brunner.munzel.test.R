#' Brunner--Munzel Test for Stochastic Equality
#' 
#' The Brunner--Munzel test for stochastic equality of two samples, 
#' which is also known as the Generalized Wilcoxon test. 
#' \code{NA}s from the data are omitted.
#' 
#' @details There exist discrepancies with \insertCite{Brunner_Munzel_2000;textual}{lawstat}
#' because there is a typo in the paper. The corrected version is in 
#' \insertCite{Neubert_Brunner_2007;textual}{lawstat}
#' (e.g., compare the estimates for the case study on pain scores). 
#' The current function follows \insertCite{Neubert_Brunner_2007;textual}{lawstat}.
#' 
#' 
#' @param x the numeric vector of data values from the sample 1.
#' @param y the numeric vector of data values from the sample 2.
#' @param alpha significance level, default is 0.05 for 95\% confidence interval.
#' @param alternative a character string specifying the alternative hypothesis, 
#' must be one of \code{"two.sided"} (default), \code{"greater"} or 
#' \code{"less"}. User can specify just the initial letter.
#'
#'
#' @return A list of class \code{"htest"} with the following components:
#' \item{statistic}{the Brunner--Munzel test statistic.}
#' \item{parameter}{the degrees of freedom.}
#' \item{conf.int}{the confidence interval.}
#' \item{p.value}{the \eqn{p}-value of the test.}
#' \item{data.name}{a character string giving the name of the data.}
#' \item{estimate}{an estimate of the effect size, i.e., \eqn{P(X < Y) + 0.5 P(X =Y )}.}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @seealso \code{\link[stats]{wilcox.test}}, \code{\link[stats]{pwilcox}}
#' 
#' @keywords htest nonparametric
#' 
#' @author Wallace Hui, Yulia R. Gel, Joseph L. Gastwirth, Weiwen Miao. 
#' This function was updated with the help of Dr. Ian Fellows.
#' 
#' @export
#' @examples
#' ## Pain score on the third day after surgery for 14 patients under
#' ## the treatment Y and 11 patients under the treatment N
#' ## (see Brunner and Munzel, 2000; Neubert and Brunner, 2007).
#' 
#' Y <- c(1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 4, 1, 1)
#' N <- c(3, 3, 4, 3, 1, 2, 3, 1, 1, 5, 4)
#' 
#' brunner.munzel.test(Y, N)
#' 
#' ##       Brunner-Munzel Test
#' ## data: Y and N
#' ## Brunner-Munzel Test Statistic = 3.1375,  df = 17.683, p-value = 0.005786
#' ## 95 percent confidence interval:
#' ##  0.5952169 0.9827052
#' ## sample estimates:
#' ## P(X<Y)+.5*P(X=Y)
#' ##        0.788961
#' 
brunner.munzel.test <- function (x, y, alternative = c("two.sided", "greater", "less"), 
                                 alpha = 0.05) 
{
    alternative <- match.arg(alternative)
    DNAME = paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    x <- na.omit(x)
    y <- na.omit(y)
    n1 = length(x)
    n2 = length(y)
    r1 = rank(x)
    r2 = rank(y)
    r = rank(c(x, y))
    m1 = mean(r[1:n1])
    m2 = mean(r[n1 + 1:n2])
    pst = (m2 - (n2 + 1)/2)/n1
    v1 = sum((r[1:n1] - r1 - m1 + (n1 + 1)/2)^2)/(n1 - 1)
    v2 = sum((r[n1 + 1:n2] - r2 - m2 + (n2 + 1)/2)^2)/(n2 - 1)
    statistic = n1 * n2 * (m2 - m1)/(n1 + n2)/sqrt(n1 * v1 + 
                                                       n2 * v2)
    dfbm = ((n1 * v1 + n2 * v2)^2)/(((n1 * v1)^2)/(n1 - 1) + 
                                        ((n2 * v2)^2)/(n2 - 1))
    if ((alternative == "greater") | (alternative == "g")) {
        p.value = pt(statistic, dfbm)
    }
    else if ((alternative == "less") | (alternative == "l")) {
        p.value = 1-pt(statistic, dfbm)
    }
    else {
        alternative = "two.sided"
        p.value = 2 * min(pt(abs(statistic), dfbm), (1 - pt(abs(statistic), 
                                                            dfbm)))
    }
    conf.int = c(pst - qt(1 - alpha/2, dfbm) * sqrt(v1/(n1 * 
                                                            n2^2) + v2/(n2 * n1^2)), pst + qt(1 - alpha/2, dfbm) * 
                     sqrt(v1/(n1 * n2^2) + v2/(n2 * n1^2)))
    estimate = pst
    ESTIMATE = pst
    names(ESTIMATE) = "P(X<Y)+.5*P(X=Y)"
    STATISTIC = statistic
    names(STATISTIC) = "Brunner-Munzel Test Statistic"
    PARAMETER = dfbm
    names(PARAMETER) = "df"
    CONF.INT = conf.int
    names(CONF.INT) = c("lower", "upper")
    attr(CONF.INT, "conf.level") = (1 - alpha)
    METHOD = "Brunner-Munzel Test"
    structure(list(estimate = ESTIMATE, conf.int = CONF.INT, 
                   statistic = STATISTIC, parameter = PARAMETER, p.value = p.value, 
                   method = METHOD, data.name = DNAME), class = "htest")
}