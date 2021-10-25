#' Coefficient of Dispersion -- a Measure of Relative Variability
#' 
#' Measure of relative inequality (or relative variation) of the data. 
#' Coefficient of dispersion (CD) is the ratio of the mean absolute deviation from 
#' the median (MAAD) to the median of the data. \code{NA}s from the data are omitted.
#' See \insertCite{Gastwirth_1988v1;textual}{lawstat}  and
#' \insertCite{Bonett_Seier_2006;textual}{lawstat}.
#'
#'
#' @param x a numeric vector of data values.
#'
#'
#' @return The coefficient of dispersion.
#' 
#' @references
#' \insertAllCited{}
#' 
#' @seealso \code{\link{gini.index}}, \code{\link{j.maad}}
#' 
#' @keywords homogeneity variability
#' 
#' @author Wallace Hui, Yulia R. Gel, Joseph L. Gastwirth, Weiwen Miao
#' 
#' @export
#' @examples
#' ## The Baker v. Carr Case: one-person-one-vote decision. 
#' ## Measure of Relative Inequality of Population data in 33 districts 
#' ## of the Tennessee Legislature in 1900 and 1972. See 
#' ## popdata (see Gastwirth, 1988).
#' 
#' data(popdata)
#' cd(popdata[,"pop1900"])
#' cd(popdata[,"pop1972"])
#' 
`cd` <- function(x)
{
    x <- na.omit(x)
    x <- sort(x)
    n <- length(x)
    M <- median(x)
    coef <- mean(abs(x - M))/M
    return(coef)
}

