#' Runs Test for Randomness
#' 
#' Performs the runs test for randomness \insertCite{Mendenhall_Reinmuth_1982}{lawstat}. 
#' Users can choose whether to plot the 
#' correlation graph or not, and whether to test against two-sided, negative, 
#' or positive correlation. \code{NA}s from the data are omitted.
#' 
#' @details On the graph, observations that are less than the sample median are 
#' represented by red letters "A", and observations that are greater or equal to the 
#' sample median are represented by blue letters "B".
#'
#'
#' @param y a numeric vector of data values.
#' @param plot.it logical. If \code{TRUE}, then the graph will be plotted. 
#' If \code{FALSE} (default), then it is not plotted.
#' @param alternative a character string specifying the alternative hypothesis, 
#' must be one of \code{"two.sided"} (default), \code{"negative.correlated"}, 
#' or \code{"positive.correlated"}.
#'
#'
#' @return A list of class \code{"htest"} with the following components:
#' \item{statistic}{the value of the standardized runs statistic.}
#' \item{p.value}{the \eqn{p}-value for the test.}
#' \item{data.name}{a character string giving the names of the data.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @seealso \code{\link{bartels.test}}
#' 
#' @keywords distribution htest
#' 
#' @author Wallace Hui, Yulia R. Gel, Joseph L. Gastwirth, Weiwen Miao
#' 
#' @export
#' @examples
#' ##Simulate 100 observations from an autoregressive model 
#' ## of the first order (AR(1))
#' y = arima.sim(n = 100, list(ar = c(0.5)))
#' 
#' ##Test y for randomness
#' runs.test(y)
#' 
#' 
`runs.test` <- function(y,
                        plot.it = FALSE,
                        alternative = c("two.sided", "positive.correlated", "negative.correlated"))
{
    alternative <- match.arg(alternative)
    DNAME = deparse(substitute(y))
    ##Strip NAs
    y <- na.omit(y)
    ### Calculate the runs of the data ###
    med <- median(y, na.rm = TRUE)
    for (k in 2:length(y)) {
        if ((y[k] == med) & (y[k - 1] < med)) {
            y[k] = y[k - 1]
        } else if ((y[k] == med) & (y[k - 1] > med)) {
            y[k] = y[k - 1]
        }
    }
    q <- rep(0.05, length(y))
    p <- rep(-0.05, length(y))
    d <- y
    q[I(d < med) | I(d == med)] <- NA
    p[I(d >= med)] <- NA
    
    ## User can select whether to plot Graph or not. Default is not plotted ##
    if (plot.it) {
        plot(
            q,
            type = "p",
            pch = "A",
            col = 'red',
            ylim = c(-0.5, 0.5),
            xlim = c(1, length(y)),
            xlab = "",
            ylab = ""
        )
        points(p, pch = "B", col = 'blue')
        abline(h = 0)
    }
    m <- length(na.omit(q))
    n <- length(na.omit(p))
    R <- 1
    s <- sign(y - med)
    for (k in 1:(length(y) - 1)) {
        if (s[k] != s[k + 1]) {
            R <- R + 1
        }
    }
    E <- 1 + 2 * n * m / (n + m)
    s2 <- (2 * n * m * (2 * n * m - n - m)) / ((n + m) ^ 2 * (n + m - 1))
    statistic <- (R - E) / sqrt(s2)
    ### One-sided or Two-sided Test ###
    ### Users will select the test method. Two sided test is default ###
    if (alternative == "positive.correlated") {
        p.value = pnorm(statistic)
        METHOD = "Runs Test - Positive Correlated"
    } else if (alternative == "negative.correlated") {
        p.value = 1 - pnorm(statistic)
        METHOD = "Runs Test - Negative Correlated"
    } else {
        p.value = 2 * min(pnorm(statistic), 1 - pnorm(statistic))
        alternative = "two.sided"
        METHOD = "Runs Test - Two sided"
    }
    ### Display Output ###
    STATISTIC = statistic
    names(STATISTIC) = "Standardized Runs Statistic"
    structure(
        list(
            statistic = STATISTIC,
            p.value = p.value,
            method = METHOD,
            data.name = DNAME
        ),
        class = "htest"
    )
}
