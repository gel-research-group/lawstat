#' The Cochran--Mantel--Haenszel Chi-square Test
#' 
#' The Cochran--Mantel--Haenszel (CMH) procedure tests homogeneity of population 
#' proportions after taking into account other factors. This procedure is widely 
#' used in law cases, for example, on equal employment and discrimination, 
#' and in biological and phamaceutical studies.
#' 
#' @details The test is based on the CMH procedure discussed 
#' by \insertCite{Gastwirth_1984;textual}{lawstat}. The data should be input in an array 
#' of 2 rows \eqn{\times} 2 columns \eqn{\times} \eqn{k} levels. 
#' The output includes the Mantel--Haenszel Estimate, the pooled Odd Ratio, 
#' and the Odd Ratio between the rows and columns at each level. The Chi-square 
#' test of significance tests if there is an interaction or association 
#' between rows and columns.
#' 
#' The null hypothesis is that the pooled Odd Ratio is equal to 1, i.e., 
#' there is no interaction between rows and columns. For more details see 
#' \insertCite{Gastwirth_1984;textual}{lawstat}.
#' 
#' The \code{cmh.test} can be viewed as a subset of 
#' \code{\link[stats]{mantelhaen.test}}, in the sense that \code{cmh.test} is for a 
#' 2 by 2 by \eqn{k} table without continuity correction, whereas 
#' \code{\link[stats]{mantelhaen.test}} allows for a larger table, 
#' and for a 2 by 2 by \eqn{k} table, it has an option of performing continuity correction. 
#' However, in view of \insertCite{Gastwirth_1984;textual}{lawstat}, continuity 
#' correction is not recommended as it tends to overestimate the \eqn{p}-value.
#' 
#'
#' @param x a numeric \eqn{2 \times 2 \times k} array of data values.
#'
#'
#' @return A list of class \code{"htest"} containing the following components:
#' \item{MH.ESTIMATE}{the value of the Cochran--Mantel--Haenszel estimate.}
#' \item{OR}{pooled Odd Ratio of the data.}
#' \item{ORK}{vector of Odd Ratio of each level.}
#' \item{cmh}{the test statistic.}
#' \item{df}{degrees of freedom.}
#' \item{p.value}{the \eqn{p}-value of the test.}
#' \item{method}{type of the performed test.}
#' \item{data.name}{a character string giving the name of the data.}
#' 
#' @references
#' \insertAllCited{}
#' 
#' @seealso \code{\link[stats]{mantelhaen.test}}
#' 
#' @keywords htest homogeneity
#' 
#' @author Min Qin, Wallace W. Hui, Yulia R. Gel, Joseph L. Gastwirth
#' 
#' @export
#' @examples
#' ## Sample Salary Data
#' data(blackhire)
#' cmh.test(blackhire)
#' 
`cmh.test` <- function(x)
{
    pooled = apply(x, 1:2, sum)
    OR = pooled[1,1] * pooled[2,2] / pooled[1,2] / pooled[2,1]
    k = dim(x)[3]
    n11k = x[1,1,]
    n21k = x[2,1,]
    n12k = x[1,2,]
    n22k = x[2,2,]
    ORK = x[1,1,]*x[2,2,]/x[1,2,]/x[2,1,]
    row1sums = n11k + n12k
    row2sums = n21k + n22k
    col1sums = n11k + n21k
    n = apply(x,3,sum)
    
    u11 = row1sums*col1sums/n
    var11 = row1sums*row2sums*col1sums*(n-col1sums)/(n^2)/(n-1)
    num = (sum(n11k - u11))^2
    deno = sum(var11)
    cmh = num/deno
    cmh.p.value = 1 - pchisq(cmh,1)
    
    DNAME = deparse(substitute(x))
    METHOD = "Cochran-Mantel-Haenszel Chi-square Test"
    
    s.diag <- sum(x[1, 1, ] * x[2, 2, ]/n)
    s.offd <- sum(x[1, 2, ] * x[2, 1, ]/n)
    MH.ESTIMATE <- s.diag/s.offd
    
    orkname = paste("Odd Ratio of level", 1:k)
    
    PARAMETER = c(cmh, 1, cmh.p.value, MH.ESTIMATE, OR, ORK)
    names(PARAMETER) = c("CMH statistic", "df", "p-value", 
                         "MH Estimate", "Pooled Odd Ratio", orkname)
    structure(list(parameter = PARAMETER, method = METHOD, data.name = DNAME), 
              class = "htest")
}

