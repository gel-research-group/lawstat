#' Goodness-of-fit Test Statistics for the Laplace Distribution
#' 
#' Goodness-of-fit test statistics 
#' \code{A2} (Anderson--Darling), 
#' \code{W2} (Cramer--von Mises), 
#' \code{U2} (Watson), 
#' \code{D} (Kolmogorov--Smirnov), and 
#' \code{V} (Kuiper). 
#' By default, \code{NA}s are omitted. For the tables of critical values, see 
#' \insertCite{Stephens_1986;textual}{lawstat} and 
#' \insertCite{Puig_Stephens_2000;textual}{lawstat}.
#'
#' @details The function originally used \code{plaplace} function from R package \code{VGAM}
#' \insertCite{VGAM}{lawstat}, however, to resolve dependencies between packages, 
#' the \code{plaplace} function was copied entirely to the current package under the name \code{VGAM_plaplace}.
#'
#' @param y a numeric vector of data values.
#'
#'
#' @return A list with the following numeric components:
#' \item{A2}{the Anderson--Darling statistic.}
#' \item{W2}{the Cramer--von Mises statistic.}
#' \item{U2}{the Watson statistic.}
#' \item{D}{the Kolmogorov--Smirnov statistic.}
#' \item{V}{the Kuiper statistic.}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @seealso \code{\link[VGAM]{plaplace}}
#' 
#' @keywords distribution
#' 
#' @author Kimihiro Noguchi, Yulia R. Gel
#' 
#' @export
#' @examples
#' ## Differences in flood levels example taken from Puig and Stephens (2000)
#' y <- c(1.96,1.97,3.60,3.80,4.79,5.66,5.76,5.78,6.27,6.30,6.76,7.65,7.84,7.99,8.51,9.18,
#'      10.13,10.24,10.25,10.43,11.45,11.48,11.75,11.81,12.33,12.78,13.06,13.29,13.98,14.18,
#'      14.40,16.22,17.06)
#' laplace.test(y)$D
#' ## [1] 0.9177726
#' ## The critical value at the 0.05 significance level is approximately 0.906.
#' ## Thus, the null hypothesis should be rejected at the 0.05 level.
`laplace.test` <- function(y)
{
    y <- y[!is.na(y)]
    y <- sort(y)
    n <- length(y)
    a <- median(y)
    b <- mean(abs(y-a))
    z <- VGAM_plaplace((y - a)/b)
    ### Anderson-Darling statistic ###
    A2 <- -mean((2 * seq(1:n) - 1) * (log(z) + log(1 - rev(z))))-n
    ### Cramer-von Mises statistic ###
    W2 <- sum((z-(2*seq(1:n)-1)/(2*n))^2)+1/(12*n)
    ### Watson statistic ###
    U2 <- W2-n*(mean(z)-0.5)^2
    ### Kolmogorov statistics (D and V) ###
    D <- sqrt(n)*max(max(seq(1:n)/n-z), max(z-(seq(1:n)-1)/n))
    V <- sqrt(n)*(max(seq(1:n)/n-z)+max(z-(seq(1:n)-1)/n))
    ### display output ###
    list(A2 = A2, W2 = W2, U2 = U2, D = D, V = V)
}

