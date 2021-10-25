#' Test of Symmetry
#' 
#' Perform test for symmetry about an unknown median. Users can choose among the 
#' Cabilio--Masaro test \insertCite{Cabilio_Masaro_1996}{lawstat},
#' the Mira test \insertCite{Mira_1999}{lawstat}, 
#' or the MGG test \insertCite{Miao_etal_2006}{lawstat}; 
#' and between using asymptotic distribution of the respective statistics or 
#' a distribution from \eqn{m}-out-of-\eqn{n} bootstrap 
#' \insertCite{Lyubchich_etal_2016_symmetry}{lawstat}.
#' Additionally to the general distribution asymmetry, the function allows to test 
#' for negative or positive skeweness (see the argument \code{side}). 
#' \code{NA}s from the data are omitted.
#' 
#' @details If the bootstrap option is used (\code{boot = TRUE}), a bootstrap 
#' distribution is obtained for each candidate subsample size \eqn{m}. Then, a heuristic 
#' method \insertCite{Bickel_etal_1997,Bickel_Sakov_2008}{lawstat}
#' is used for the choice of optimal \eqn{m}. Specifically, we use the Wasserstein metric 
#' \insertCite{Ruschendorf_2001}{lawstat} to calculate distances between different 
#' bootstrap distributions and select \eqn{m}, which corresponds to the minimal distance. 
#' See \insertCite{Lyubchich_etal_2016_symmetry;textual}{lawstat} for more details.
#'
#'
#' @param x data to be tested for symmetry.
#' @param option test statistic to be applied. The options include statistic by 
#' \insertCite{Miao_etal_2006;textual}{lawstat} (default),
#' \insertCite{Cabilio_Masaro_1996;textual}{lawstat}, 
#' and \insertCite{Mira_1999;textual}{lawstat}.
#' @param side choice from the three possible alternative hypotheses: 
#' general distribution asymmetry (\code{side = "both"}, default), 
#' left skewness (\code{side = "left"}), or right skewness (\code{side = "right"}).
#' @param boot logical value indicates whether \eqn{m}-out-of-\eqn{n} bootstrap will 
#' be used to obtain critical values (default), 
#' or asymptotic distribution of the chosen statistic.
#' @param B number of bootstrap replications to perform (default is 1000).
#' @param q scalar from 0 to 1 to define a set of possible \eqn{m} for the 
#' \eqn{m}-out-of-\eqn{n} bootstrap. Default \code{q = 8/9}. 
#' Possible \eqn{m} are then set as the values \code{unique(round(n*(q^j))} 
#' greater than 4, where \code{n = length(x)} and \code{j = c(0:20)}.
#'
#'
#' @return A list of class \code{"htest"} with the following components:
#' \item{method}{name of the method.}
#' \item{data.name}{name of the data.}
#' \item{statistic}{value of the test statistic.}
#' \item{p.value}{\eqn{p}-value of the test.}
#' \item{alternative}{alternative hypothesis.}
#' \item{estimate}{bootstrap optimal \eqn{m} (given in the output only if bootstrap 
#' was used, i.e., \code{boot = TRUE}).}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @keywords distribution htest robust
#' 
#' @author Joseph L. Gastwirth, Yulia R. Gel, Wallace Hui, Vyacheslav Lyubchich, 
#' Weiwen Miao, Xingyu Wang (in alphabetical order)
#' 
#' @export
#' @examples
#' data(zuni) #run ?zuni to see the data description
#' symmetry.test(zuni[,"Revenue"], boot = FALSE)
#' 
symmetry.test <- function(x,
                          option = c("MGG", "CM", "M"),
                          side = c("both", "left", "right"),
                          boot = TRUE,
                          B = 1000,
                          q = 8 / 9) 
{
    DNAME <- deparse(substitute(x))
    j <- c(0:20)
    x <- na.omit(x)
    n <- length(x)
    m <- unique(round(n * (q ^ j))[I(round(n * (q ^ j)) > 4)])
    x <- sort(x) - mean(x)
    a <- x[x < mean(x)]
    b <- x[x > mean(x)]
    xx <- c(x, -a, -b)
    option <- match.arg(option)
    side <- match.arg(side)
    if (option == "MGG") {
        METHOD <- "test by Miao, Gel, and Gastwirth (2006)"
        stat.function <-
            function(x)
                sqrt(length(x)) * (mean(x) - median(x)) / (sqrt(pi / 2) * mean(abs(x - median(x))) * sqrt(0.5707963))
    }
    if (option == "M") {
        METHOD <- "test by Mira (1999)"
        stat.function <- function(x) {
            x <- sort(x)
            N <- length(x)
            D <-
                N ^ (1 / 5) * (x[N / 2 + 0.5 * N ^ (4 / 5)] - x[N / 2 - 0.5 * N ^ (4 /
                                                                                       5) + 1])
            g <- mean(x) - 2 / N * sum(x[x <= median(x)])
            S <- 4 * var(x) + D ^ 2 - 4 * D * g
            sqrt(N) * 2 * (mean(x) - median(x)) / sqrt(S)
        }
    }
    if (option == "CM") {
        METHOD <- "test by Cabilio and Masaro (1996)"
        stat.function <-
            function(x)
                sqrt(length(x)) * (mean(x) - median(x)) / (sd(x) * sqrt(0.5707963))
    }
    STATISTIC <- stat.function(x)
    names(STATISTIC) <- "Test statistic"
    if (boot) {
        BootstrapStatistic <- array(NA, c(length(m), B))
        for (i in 1:length(m)) {
            M <-
                sapply(1:B, function(x)
                    sort(sample(
                        xx, size = m[i], replace = TRUE
                    )))
            tmp <- sapply(1:B, function(x)
                length(unique(M[, x])))
            while (any(tmp == 1)) {
                M[, tmp == 1] <-
                    sapply(which(tmp == 1), function(x)
                        sort(sample(
                            xx, size = m[i], replace = TRUE
                        )))
                tmp <- sapply(1:B, function(x)
                    length(unique(M[, x])))
            }
            BootstrapStatistic[i, ] <- sort(apply(M, 2, stat.function))
        }
        Distance <-
            sapply(1:(length(m) - 1), function(x)
                sqrt(sum((BootstrapStatistic[(x + 1), ] - BootstrapStatistic[x, ]) ^ 2
                )))
        di <- which.min(Distance)
        ESTIMATE <- m[di]
        names(ESTIMATE) <- "bootstrap optimal m"
        Tcrit <- sum(STATISTIC < BootstrapStatistic[di, ]) / B
    } else {
        Tcrit <- 1 - pnorm(STATISTIC)
    }
    if (side == "both") {
        ALTERNATIVE <- "the distribution is asymmetric."
        if (Tcrit < 0.5) {
            P.VALUE <- Tcrit * 2
        } else {
            P.VALUE <- 2 * (1 - Tcrit)
        }
    }
    if (side == "left") {
        ALTERNATIVE <- "the distribution is negatively skewed."
        P.VALUE <- 1 - Tcrit
    }
    if (side == "right") {
        ALTERNATIVE <- "the distribution is positively skewed."
        P.VALUE <- Tcrit
    }
    if (boot) {
        METHOD <- paste("m-out-of-n bootstrap symmetry", METHOD)
        structure(
            list(
                method = METHOD,
                data.name = DNAME,
                statistic = STATISTIC,
                p.value = P.VALUE,
                alternative = ALTERNATIVE,
                estimate = ESTIMATE
            ),
            class = "htest"
        )
    } else {
        METHOD <- paste("Symmetry", METHOD)
        structure(
            list(
                method = METHOD,
                data.name = DNAME,
                statistic = STATISTIC,
                p.value = P.VALUE,
                alternative = ALTERNATIVE
            ),
            class = "htest"
        )
    }
}