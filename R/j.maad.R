#' MAAD Robust Standard Deviation
#' 
#' Compute average absolute deviation from the sample median, 
#' which is a consistent robust estimate of the population standard deviation 
#' for normally distribution data \insertCite{Gastwirth_1982}{lawstat}. 
#' \code{NA}s from the data are omitted.
#'
#'
#' @param x a numeric vector of data values.
#'
#'
#' @return Robust standard deviation.
#' 
#' @references
#' \insertAllCited{}
#' 
#' @seealso \code{\link{cd}}, \code{\link{gini.index}}, \code{\link{rqq}}, 
#' \code{\link{rjb.test}}, \code{\link{sj.test}}
#' 
#' @keywords homogeneity robust variability
#' 
#' @author Wallace Hui, Yulia R. Gel, Joseph L. Gastwirth, Weiwen Miao
#' 
#' @export
#' @examples
#' ## Sample 100 observations from the standard normal distribution
#' x = rnorm(100)
#' j.maad(x)
#' 
`j.maad` <- function(x) 
{
    x <- na.omit(x)
    ### Robust Standard Deviation J
    J <- sqrt(pi/2)*mean(abs(x-median(x))) 
    return(J)
}

