#' Measures of Relative Variability -- Gini Index
#' 
#' Gini index for measuring relative inequality (or relative variation) of the data 
#' \insertCite{Gini_1912}{lawstat}. \code{NA}s from the data are omitted.
#' 
#' @details See also \insertCite{Gastwirth_1988v1;textual}{lawstat}.
#'
#'
#' @param x the input data.
#'
#'
#' @return A list with the following components:
#' \item{statistic}{the Gini index.}
#' \item{parameter}{the mean difference of the set of numbers.}
#' \item{data.name}{a character string giving the name of the data.}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @seealso \code{\link{cd}}, \code{\link{j.maad}}, \code{\link{lorenz.curve}}
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
#' ## popdata (see Gastwirth (1988)).
#' data(popdata)
#' gini.index(popdata[,"pop1900"])
#' gini.index(popdata[,"pop1972"])
#' 
`gini.index` <- function(x)
{
    DNAME = deparse(substitute(x))
    x <- na.omit(x)
    x = sort(x)
    n = length(x)
    ### Calculate the delta and Gini Index ###
    a <- 0
    for (i in 1:(n-1)) {
        for (j in (i+1):n) {
            a <- a + abs(x[j]- x[i])
        }
    }
    delta <- 2*a/n/(n-1)
    GI <- delta/(2*mean(x))
    METHOD = "Measures of Relative Variability - Gini Index"
    ### Display Output ###
    STATISTIC = GI
    names(STATISTIC) = "Gini Index"
    PARAMETER = delta
    names(PARAMETER) = "delta"
    structure(list(statistic = STATISTIC, parameter = PARAMETER, method = METHOD, 
                   data.name = DNAME), 
              class = "htest")
    
}

